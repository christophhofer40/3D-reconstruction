package christoph;
import java.util.ArrayList;
//import christoph.Quaternion;
import java.util.Arrays;

//fits a parabel on 2 points, each representing the gradient
//the minimum of the parabel is assumed as the optimized position
public class Optimization
{
    public ArrayList<Twopoints> points=new ArrayList<Twopoints>(10);
    
    public Optimization(float xcenter, float ycenter, float zcenter, float[] p1, float p1grad, float[] p2, float p2grad)
    {
        this.points.add(new Twopoints(xcenter,ycenter,zcenter,  p1, p1grad, p2, p2grad));
    }
    
    public void add_two_points(float[] center, float[] p1, float p1grad, float[] p2, float p2grad)
    {
        points.add(new Twopoints(center[0],center[1],center[2],p1, p1grad, p2, p2grad));
    }
    
    @Override
    public String toString()
    {
        return "new point: "+points.get(0).optimum;
    }
    
    public float[] get_shift(int pointsindex)
    {
        return this.points.get(pointsindex).shift;
    }
    
    
    
    class Twopoints
    {
        float tolerance=0.001f;
        float[] p1=new float[3];
        float p1grad;
        float[] p2=new float[3];
        float p2grad;
        float[] center=new float[4];
        float[] directionVector= new float[3]; //p1-p2 normed
        float optimum;
        float[] newPosition=new float[3];
        float shift[] = new float[3];
        
        Twopoints(float xcenter, float ycenter, float zcenter, float[] p1, float p1grad, float[] p2, float p2grad)
        {
            center[0]=xcenter;
            center[1]=ycenter;
            center[2]=zcenter;
            this.p1=p1;
            this.p1grad=p1grad;
            this.p2=p2;
            this.p2grad=p2grad;
            optimum=find_optimum();
            calc_new_position();
            calc_shift();
            //System.out.println(optimum);
            //System.out.println(Arrays.toString(newPosition));
            
        }
        
        boolean get_directionVector()
        {
            float norm=(float) Math.sqrt((p2[0]-p1[0])*(p2[0]-p1[0])+
                    (p2[1]-p1[1])*(p2[1]-p1[1])+(p2[2]-p1[2])*(p2[2]-p1[2]));
            float[] controll= new float[3];
            controll[0]=-(p2[0]-p1[0])/norm;
            controll[1]=-(p2[1]-p1[1])/norm;
            controll[2]=-(p2[2]-p1[2])/norm;
            
            
            norm=(float) Math.sqrt((center[0]-p1[0])*(center[0]-p1[0])+
                    (center[1]-p1[1])*(center[1]-p1[1])+(center[2]-p1[2])*(center[2]-p1[2]));
            directionVector[0]=(p1[0]-center[0])/norm;
            directionVector[1]=(p1[1]-center[1])/norm;
            directionVector[2]=(p1[2]-center[2])/norm;
            if((directionVector[0]-controll[0])>tolerance||(directionVector[1]-controll[1])>tolerance||(directionVector[2]-controll[2])>tolerance)
            {
                System.out.println("Points not colinear!");
                return false;
            }
            else
                return true;
        }
        /*
        float[] get_orthogonalVector(float[] vec)
        {
            return new float[]{vec[1],-vec[0],0};
        }
        
        float[] project_point(float[] point)
        {
            float[] ret;
            Quaternion q= new Quaternion(0,get_orthogonalVector(directionVector));
            ret= q.rotate(point);
            return new float[]{ret[0],ret[1]};
        }*/
        
        
        float CalcParabolaVertex(float x1, float y1, float x2, float y2, float x3, float y3)
        {
            float denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
            float A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
            float B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
            float C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

            return (-B / (2*A));
             
        }
        
        final float find_optimum()
        {
            boolean colinear=get_directionVector();
            //System.out.println(Arrays.toString(directionVector));
            //directionVector=get_orthogonalVector(directionVector);
            float lp1=(float) Math.sqrt((center[0]-p1[0])*(center[0]-p1[0])+(center[1]-p1[1])*(center[1]-p1[1])+(center[2]-p1[2])*(center[2]-p1[2]));
            float lp2= (float) Math.sqrt((center[0]-p2[0])*(center[0]-p2[0])+(center[1]-p2[1])*(center[1]-p2[1])+(center[2]-p2[2])*(center[2]-p2[2]));
            //System.out.println(lp1+" "+lp2);
            return CalcParabolaVertex(0.0f,0.0f,lp1,p1grad,-lp2,p2grad);
            
            
        }
        
        final void calc_new_position()
        {
            newPosition[0]=center[0]+directionVector[0]*optimum;
            newPosition[1]=center[1]+directionVector[1]*optimum;
            newPosition[2]=center[2]+directionVector[2]*optimum;
        }
        
        final void calc_shift()
        {
            shift[0]=directionVector[0]*optimum;
            shift[1]=directionVector[1]*optimum;
            shift[2]=directionVector[2]*optimum;
        }
        
    }
    //test sampler
    public static void main(String[] args)
    {
        Optimization opt= new Optimization(2.0f,2.0f,2.0f,
                new float[]{3.0f,3.0f,3.0f},-1.0f,
                new float[]{1.0f,1.0f,1.0f},-0.05f);
    }
}
