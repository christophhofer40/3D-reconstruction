package christoph;

import ij.*;
import java.util.Arrays;
import java.util.Collections;

public class ImageCalculations
{
   
    public static float[] normalize(float[] pixels,double norm)
    {
            double sum = 0.0f;
            final int area = pixels.length;
            for(int ind = 0; ind < area; ++ind)
            {
                    sum += pixels[ind];	
            }
            double asum = Math.abs(sum); //normalize negative images to -1 * norm
            asum = (float)(norm/asum);
            for(int ind = 0; ind < area; ++ind)
            {
                    pixels[ind]*=asum;	
            }
            return pixels;
    }
    public static float[] normalize_image(float[] im, float[] mask)
    {
        int len=im.length;
        double mean=get_mean(im, mask);
        for(int i=0; i<len;++i)
        {
            if(!Float.isNaN(mask[i]))
                im[i]-=mean;
            else
            	im[i]=Float.NaN;
        }
        double sdt=get_std(im,0.0,mask);
        for(int i=0; i<len;++i)
        {
            if(!Float.isNaN(mask[i]))
                im[i]/=sdt;
        }
        
        return im;
    }
    
    public static float[] sum_arrays(float[] a1, float[] a2 )
    {
        float[] ret=new float[a1.length];
        for(int i=0; i<ret.length;++i)
        {
            if(!Float.isNaN(a1[i]) && !Float.isNaN(a2[i]))
            {
                ret[i]=a1[i]+a2[i];
            }
            else
            	ret[i]=Float.NaN;
        }
        return ret;
    }
    
    static double get_mean(float[] ar,float[] mask)
    {
        double sum=0;
        int counter=0;
        int len=ar.length;
        for(int i=0; i<len;++i)
        {
            if(!Float.isNaN(mask[i]))
            {
                sum+=ar[i];
                counter+=1;
            }
            
        }
        return sum/counter;
    }
    
    public static double get_std(float[] ar,double mean, float[] mask)
    {
        double sum2=0;
        int counter=0;
        int len=ar.length;
        for(int ind=0; ind < len; ++ind)
        {
            if(!Float.isNaN(mask[ind]))
            {
                final double v = ar[ind]-mean; 
                counter+=1;
                sum2 += v*v;
            }
        }
        return(Math.sqrt(sum2/counter));
    }
    
    public static float get_minarg(float[] ar)
    {
        int counter=0;
        int len=ar.length;
        float min=ar[0];
        int minind=0;
        for(int ind=1; ind < len; ++ind)
        {
            if (ar[ind]<min)
            {
                min=ar[ind];
                minind=ind;
            }
        }
        return minind;
    }
    //offset= n --> n-th extrema
    public static int get_maxarg(float[] ar)
    {
        int counter=0;
        int len=ar.length;
        float max=ar[0];
        int maxind=0;
        for(int ind=1; ind < len; ++ind)
        {
            if (ar[ind]>max)
            {
                max=ar[ind];
                maxind=ind;
            }
        }
        return maxind;
    }
    public static int get_maxarg(float[] ar,int offset)
    {
        int counter=0;
        int len=ar.length;
        int ind;
        float[] cp=(float[])ar.clone();
        final int index = Arrays.asList(ar).indexOf(cp[len-offset-1]);
        return index;
    }
    
    public static void plot_image(float[] array, String title)
    {
        int length= (int)Math.sqrt(array.length);
        ImageStack stack= new ImageStack(length,length,1);
        stack.setPixels(array, 1);
        ImagePlus imp=new ImagePlus(title,stack);
        imp.show();
    }
    
    public static double img_match(float[] img1, float[] img2)
    {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double std1 = 0.0;
        double std2 = 0.0;
        double mix = 0.0;
        int counts=0;
        int area = img1.length;

        for(int ind = 0; ind < area; ++ind)
        {
            if(!Float.isNaN(img1[ind]) && !Float.isNaN(img2[ind]))
            {
                counts+=1;
                sum1 += img1[ind];
                sum2 += img2[ind];
            }
        }


        final float avg1 = (float)( sum1 / (  counts) );
        final float avg2 = (float)( sum2 / ( counts ) );

        for(int ind = 0; ind < area; ++ind)
        {
            if(!Float.isNaN(img1[ind]) && !Float.isNaN(img2[ind]))
            {
                final float val1 = (img1[ind]-avg1);
                final float val2 = (img2[ind]-avg2);
                std1 += val1*val1;
                std2 += val2*val2;
                mix += val1*val2; 
            }

        }
       

        //final double sig1 = Math.sqrt(std1/area);
        //final double sig2 = Math.sqrt(std2/area);
        if( (std1 > 0) && (std2 > 0) )
        {
                return -mix/(Math.sqrt(std1*std2));
        }
        else return 0.0; //neutral match to a constant image
    }
    
    public static double img_match2(float[] img1, float[] img2)
    {
        double sum1 = 0.0;
        double mix = 0.0;
        int counts=0;
        int area = img1.length;

        for(int ind = 0; ind < area; ++ind)
        {
            if(!Float.isNaN(img1[ind]) && !Float.isNaN(img2[ind]))
            {
                counts+=1;
                sum1 += (img1[ind]-img2[ind])*(img1[ind]-img2[ind]);
            }
        }


        final float avg1 = (float)( sum1 / (  counts) );
        return avg1;
        
    }
    
    

}
