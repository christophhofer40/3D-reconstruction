package christoph;
import ij.*;
import static ij.measure.Measurements.MEAN;
import ij.process.*;
import ij.plugin.filter.GaussianBlur;
import ij.io.FileSaver;
import ij.io.Opener;
import java.util.ArrayList;
import java.util.Arrays;
//import org.python.core.PyException;
//import org.python.core.PyInteger;
//import org.python.core.PyObject;
//import org.python.util.PythonInterpreter;

public class Fakesimulation
{
    float[] deltalattice;
    ImageStack deltalatticeSt;
    ImageProcessor deltalatticePr;
    FHT deltafht;
    FHT simfht;
    ImageStack simSt;
    ImageProcessor simPr;
    ImagePlus simPl;
    
    ImageStack beamstack3;
    
    public Beam beam;
    GaussianBlur gb= new GaussianBlur();
    public float sigma=0;
    
    int imLength;
    int imWidth;
    public float fov;
        
    public Fakesimulation(float[] delta_graphene,float fieldOfView ,int imWidth)
    {
        this.deltalattice=delta_graphene;
        this.imLength=this.deltalattice.length;
        this.imWidth=imWidth;
        this.deltalatticeSt=new ImageStack(imWidth, imWidth,1);
        this.deltalatticeSt.setPixels(this.deltalattice, 1);
        this.deltalatticePr=this.deltalatticeSt.getProcessor(1);
        this.simSt=new ImageStack(imWidth, imWidth,1);
        this.deltafht=new FHT(this.deltalatticePr);
        this.fov=fieldOfView;
        this.beam=new Beam();
        
    }
    
    public void reset_lattice(float[] delta_graphene)
    {
        this.deltalattice=  delta_graphene;    
        this.deltalatticeSt.setPixels(this.deltalattice, 1);
        this.deltalatticePr=this.deltalatticeSt.getProcessor(1);
        //this.deltalatticePr.resize((int)imWidth/2,(int)imWidth/2,true);
        //System.out.println(((float[])deltalatticePr.getPixels()).length);
        this.deltafht=new FHT(this.deltalatticePr) ;
    }

    
    public float[] simulate_image(float[] delta_graphene)
    {
        
        reset_lattice(delta_graphene);
        // Converts lattice model to a complex Fourier transform 
        
        deltafht.swapQuadrants();
        deltafht.transform();  
        //System.out.println(Arrays.toString((float[])beam.beamfht.getPixels()));       
        beam.beamfht.swapQuadrants();
        beam.beamfht.transform();
        
        simfht=deltafht.multiply(beam.beamfht);
       
        simfht.inverseTransform();
        simfht.swapQuadrants();
        deltafht.inverseTransform();
        deltafht.swapQuadrants();
        beam.beamfht.inverseTransform();
        beam.beamfht.swapQuadrants();
        //simfht.resize((int)imWidth,(int)imWidth,true);
        simSt.setPixels(gaussianblur((float[])simfht.getPixels(),sigma), 1);
        return (float[])simSt.getPixels(1);
                    
    }
    
    float[] gaussianblur(float[] sim,float sigma_)
    {           
        ImageStack tempSt=new ImageStack(imWidth,imWidth,1);
        tempSt.setPixels(sim,1);
        ImageProcessor tempPr=tempSt.getProcessor(1);
        //System.out.println(Reconstructor.get_mean((float[])tempPr.getPixels(),(float[])tempPr.getPixels()));
        gb.blurGaussian(tempPr, sigma_);
        //System.out.println(Reconstructor.get_mean((float[])tempPr.getPixels(),(float[])tempPr.getPixels()));
        float[] ret= new float[imWidth*imWidth];
        ret=(float[]) tempPr.getPixels();
        return ret;
    }
   
    /*
    void optimize_fieldOfView()
    {
        double step=1.0;
        float[] expf=(float[])expSt.getPixels(1);
        float[] impAr= (float[])impSt.getPixels(1);
        double merit= img_match(expf,impAr);
        double newmerit=merit;
        IJ.log("adjusting field of view... ");
        IJ.log("correlation before: "+merit);
        do
        {
            fov+=step;
            simulate_image();
            //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
            newmerit=img_match(expf,(float[])impSt.getPixels(1));
            //System.out.println(newmerit);
            if (newmerit>merit)
            {
                merit=newmerit;
                step*=1.1;
            }
            else
            {
                fov-=step;
                step*=-0.7;
            }
        }
        while(step>=0.01);
        IJ.log("correlation after: "+merit);
        IJ.log("new field of view: "+fov+"\n");
    }
    
    */
    float[] absoluteArray(float[] ar)
    {
        float ret[]= new float[ar.length];
        for(int i=0;i<ar.length;++i)
        {
            ret[i]=(float)Math.abs(ar[i]);
        }
        return ret;
    }
    
    public float[] squareArray(float[] ar)
    {
        float ret[]= new float[ar.length];
        for(int i=0;i<ar.length;++i)
        {
            ret[i]=ar[i]*ar[i];
        }
        return ret;
    }
    
    
    float[] normalize_array(float[] ar)
    {
        double sum=0;
        for(int i=0;i<ar.length;++i)
        {
            sum+=ar[i];
        }
        for(int i=0;i<ar.length;++i)
        {
            ar[i]/=(float)sum;
        }
        return ar;    
    }
  
    
    float zernikePolynom(float x, float y, float[] abarray)
    {
        float val=0;
        val=(float) (-abarray[0]*(x*x+y*y)+
                Math.sqrt(abarray[1]*abarray[1]+abarray[2]*abarray[2])*(x*x+y*y)*Math.cos(2*(Math.atan2(y, x)-Math.atan2(abarray[2], abarray[1])))+
                2.0f/3.0f*Math.sqrt(abarray[3]*abarray[3]+abarray[4]*abarray[4])*4.87*0.001*
                Math.pow(Math.sqrt(x*x+y*y),3)*Math.cos(Math.atan2(y, x)-Math.atan2(abarray[4], abarray[3]))+
                2.0f/3.0f*Math.sqrt(abarray[5]*abarray[5]+abarray[6]*abarray[6])*4.87*0.001*
                Math.pow(Math.sqrt(x*x+y*y),3)*Math.cos(3*Math.atan2(y, x)-Math.atan2(abarray[6], abarray[5])))*(float)(3.1415*4.87*0.001);
        return val;
    }
    
    

    

    public class Beam
    {
        public float[] aberrations = new float[7];
        float[] kernel;
        float[] reKer;
        float[] imKer;
        float[] fhtpx;
        float aperturesize;
        float[] probetail;
        float[] probe;
        
        public ImageStack beamSt;
        ImageProcessor beamPr;
        public FHT beamfht;
        
        ImageStack complexSt;
        ImagePlus complexPl;
        
        Process process=null;
        
        boolean pythonfft=false;
        
        public Beam()
        {
            kernel=new float[imLength];
            reKer=new float[imLength];
            imKer=new float[imLength];
            probe=new float[imLength];
            probetail=new float[imLength];
            probetail[(int)imLength/2+(int)imWidth/2]=2;
            aperturesize=(0.025f/1)*(fov/(4.87f*0.001f));
            probetail=gaussianblur(probetail,aperturesize);
            //IJ.log("aperture size for 60 keV electrons and 25 mrad: "+aperturesize+" px");
            fhtpx=new float[imLength];
            beamSt=new ImageStack(imWidth,imWidth,1);
            
            complexSt=new ImageStack(imWidth, imWidth, 2);
            complexPl= new ImagePlus();
            //this.aberrations[0]=-1.0f;
            initBeam();      
        }
        
        public void initBeam()
        {
            int area=imWidth*imWidth;
            int px=imWidth/2;           
            for(int ind=0;ind<area;++ind)
            {
                int rx=ind%imWidth-px;
                //float rx2=rx/(fov*kersize);
                int ry=(int)(ind/imWidth)-px;
                //float ry2=ry/(fov*kersize);
                this.kernel[ind]=zernikePolynom(rx/fov,ry/fov,this.aberrations);
            }
            update_complex();
            fold_aperture_with_kernel();
            if(pythonfft)
            {
            	try
            	{
					Opener op=new Opener();
					
					complexSt.setPixels(reKer,1);
					complexSt.setPixels(imKer,2);
					complexPl.setStack(complexSt);
					FileSaver fs=new FileSaver(complexPl);
					complexPl.setSlice(1);
					fs.saveAsText(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/complex1.txt");
					complexPl.setSlice(2);
					fs.saveAsText(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/complex2.txt");
					//IJ.saveAsTiff(complexPl,this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/complex.tif");
					process=Runtime.getRuntime().exec("/usr/bin/python3 "+this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/fft.py");
					process.waitFor();
					if(process.exitValue()!=0)
					{
						System.out.println("Warning: Error in python script for simulation");	
					}
                                        ImageProcessor readtxt=op.openTextImage(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/","beamprofile.txt").getProcessor();
					//this.beamPr=IJ.openImage(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/temps/beamprofile.tif").getProcessor();
					this.beamPr=readtxt;
					this.beamSt.setProcessor(this.beamPr,1);
					this.beamfht= new FHT(this.beamPr,false);
            	}
            	catch(Exception e)
				{
					e.printStackTrace();
				}
            	
            }
            
            else
            {
				this.beamSt.setPixels(this.fhtpx, 1);
				//this.beamSt.setPixels(Reconstructor.sum_arrays((float[])this.beamSt.getPixels(1), probetail),1);
				this.beamPr=this.beamSt.getProcessor(1);
				//this.beamPr.resize((int)imWidth/2,(int)imWidth/2,true);
				
				this.beamfht= new FHT(this.beamPr,true);
				this.beamfht.swapQuadrants();
				
				this.beamfht.inverseTransform();
				this.beamfht.setPixels(normalize_array(squareArray((float[])beamfht.getPixels())));
				this.beamfht.swapQuadrants();
			}
			
            
            if(false) //might be included as checkbox
            {
            	beamstack3=new ImageStack(imWidth,imWidth,1);
            	beamstack3.setPixels((float[])beamfht.getPixels(), 1);
                ImagePlus beamp3=WindowManager.getImage("Electronprobe");
        
                if (beamp3==null)
                {

                        beamp3=new ImagePlus("Electronprobe",beamstack3);
                        beamp3.show();
                    
                }
                else
                {
                    beamp3.setStack(beamstack3);
                    beamp3.updateAndRepaintWindow();
                }
            }
    
            //this.beamfht.setPixels(Reconstructor.sum_arrays((float[])beamfht.getPixels(),probetail));
            //beamSt.setPixels((float[])beamfht.getPixels(), 1);
            //ImagePlus temp=WindowManager.getImage("EP1");
            //if (temp==null)
             //   Reconstructor.plot_image((float[])this.beamfht.getPixels(), "EP1");

            

        }
        
        public void update_complex()
        {
            for(int i=0; i<imLength;++i)
            {
                reKer[i]=(float)Math.cos(kernel[i]);
                imKer[i]=(float)Math.sin(kernel[i]);
                fhtpx[i]=reKer[i]-imKer[i];
            }

        }
        
        public void fold_aperture_with_kernel()
        {
            int ix;
            int iy;
            aperturesize=(0.025f/1)*(fov/(4.87f*0.001f));
            for(int i=0;i<imLength;++i)
            {
                ix=i%imWidth;
                iy=i/imWidth;
                if(Math.sqrt((iy-imWidth/2)*(iy-imWidth/2)+(ix-imWidth/2)*(ix-imWidth/2))>=aperturesize)
                {
                    kernel[i]=0;
                    reKer[i]=0;
                    imKer[i]=0;
                    fhtpx[i]=0;
                }
            }
           
        }
        
        
        void normalize()
        {
            float sum=0;
            for(int i=0;i<imLength;++i)
            {
                sum+=kernel[i];
            }
            for(int i=0;i<imLength;++i)
            {
                kernel[i]=kernel[i]/sum;
                reKer[i]/=sum;
                imKer[i]/=sum;
            }
        }
    }
    
    
    
 
}
    
