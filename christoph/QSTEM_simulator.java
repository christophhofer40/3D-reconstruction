package christoph;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.io.FileInfo;
import ij.io.FileOpener;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;


public class QSTEM_simulator
{
    public String path;
    public float modelsize; // in A
    public QscFile qsc;
    public CfgFile cfg= new CfgFile();
    int imsize;
    GaussianBlur gb= new GaussianBlur();
    
    public QSTEM_simulator(float modelsize, int impWidth, String path, String qscname)
    {
        //this.path=this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"QSTEM_files/";
        this.path=path;
        this.modelsize=modelsize;  
        this.qsc=new QscFile(qscname,imsize);
        this.imsize=impWidth;
    }
    
    public class CfgFile
    {
        
        public void write_cfg(ArrayList<float[]> rel_coordlist, boolean mirror)// throws IOException
        {
            try
            {
                File cfg = new File(path+"new.cfg");
                FileWriter cfgw= new FileWriter(cfg,false);
                PrintWriter pw= new PrintWriter(cfgw);
                pw.println("Number of particles = "+rel_coordlist.size());
                pw.println("A = 1.0 Angstrom (basic length-scale)");
                pw.println("H0(1,1) = "+modelsize+" A");
                pw.println("H0(1,2) = 0 A");
                pw.println("H0(1,3) = 0 A");
                pw.println("H0(2,1) = 0 A");
                pw.println("H0(2,2) = "+modelsize+" A");
                pw.println("H0(2,3) = 0 A");
                pw.println("H0(3,1) = 0 A");
                pw.println("H0(3,2) = 0 A");
                pw.println("H0(3,3) = "+modelsize+" A");
                pw.println(".NO_VELOCITY.");
                pw.println("entry_count = 3");
                pw.println("12");
                pw.println("C");
                if (mirror)
                    for (float[] xyz : rel_coordlist)
                        pw.println(xyz[0]+ " "+(1.0f-xyz[1])+ " "+xyz[2]);
                else
                    for (float[] xyz : rel_coordlist)
                        pw.println(xyz[0]+ " "+xyz[1]+ " "+xyz[2]);

                 pw.close();
            }
            catch(Exception e)
            {
                e.printStackTrace();
            }
        }
    }
    public class QscFile
    {
        String name;
        float[] aberrations= new float[7]; //[C10,C12a,C12b,C21a,C21b,C23a,C23b]
        File qscfile;
        int imsize;
        
        
        QscFile(String name, int imsize)
        {
            this.name=name;
            this.qscfile=new File(path+name);
            this.imsize=imsize;
        }
        
        
        void change_aberrations(float[] new_abarray) 
        {
           
            File fn = null;
            try {
                fn= new File(path+"copy.qsc");
                FileWriter fw = new FileWriter(fn,false);
                PrintWriter pw= new PrintWriter(fw);
                FileReader fr = new FileReader(qscfile);
                final Scanner s= new Scanner(fr).useDelimiter(" ");
                while(s.hasNextLine())
                {
                    String line = s.nextLine();
                    if(line.startsWith("defocus"))
                        pw.println("defocus: "+new_abarray[0]);
                    else if(line.startsWith("astigmatism: "))
                        pw.println("astigmatism: "+Math.sqrt(new_abarray[2]*new_abarray[2]+new_abarray[1]* new_abarray[1]));
                    else if(line.startsWith("astigmatism angle: "))
                        pw.println("astigmatism angle: "+Math.atan2(new_abarray[2], new_abarray[1])*180/3.1415);
                    else if(line.startsWith("a_31:"))
                        pw.println("a_31: "+Math.sqrt(new_abarray[4]*new_abarray[4]+new_abarray[3]* new_abarray[3]));
                    else if(line.startsWith("phi_31:"))
                        pw.println("phi_31: "+Math.atan2(new_abarray[4], new_abarray[3])*180/3.1415);
                    else if(line.startsWith("a_33:"))
                        pw.println("a_33: "+Math.sqrt(new_abarray[6]*new_abarray[6]+new_abarray[5]* new_abarray[5]));
                    else if(line.startsWith("phi_33:"))
                        pw.println("phi_33: "+Math.atan2(new_abarray[6], new_abarray[5])*180/3.1415);
                    else
                        pw.println(line);
                    
                }  
                pw.close();
            } catch (IOException ex) {
                Logger.getLogger(QSTEM_simulator.class.getName()).log(Level.SEVERE, null, ex);
            } finally 
            {
                try {
                    
                    Files.move(fn.toPath(),qscfile.toPath(),StandardCopyOption.ATOMIC_MOVE);
                } catch (IOException ex) {
                    Logger.getLogger(QSTEM_simulator.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            
        }
        
        
    }
    

    public float[] simulate_image(float sigma_) //throws IOException
    {
        try
        {
            Process process = new ProcessBuilder(this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/stem3d",path+qsc.name).start();
            InputStream is = process.getInputStream();
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) 
            {
                //if (line.startsWith("DEBUG")||line.startsWith("Debug")||line.startsWith("ERROR")||line.startsWith("Error"))
                System.out.println(line);
            }
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
        

        float[] arr= Read_img.read_img("plugins/qstem_images/detector1.img");
        FloatProcessor fp= new FloatProcessor(imsize,imsize,arr);
        fp.rotate(-90);
        gb.blurGaussian(fp, sigma_);
        return (float[])fp.getPixels();
        
    }

    
    
    /*void open_img() 4.895309E-13, 4.063621E-13
    {
        FileInfo fi= new FileInfo();
        fi.directory=this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"qstem_images/";
        fi.fileName="detector1.img";
        fi.width=512;
        fi.fileType = FileInfo.GRAY32_FLOAT;
        fi.intelByteOrder=true;
        fi.nImages=3;
        fi.height=512;
        fi.offset=82;
        fi.gapBetweenImages=0;
        System.out.println("try opening "+fi.directory+fi.fileName);
        FileOpener fo = new FileOpener(fi); 
        ImagePlus imp = fo.open(true); 
        if (imp==null)
            System.out.println("null simulation");
        imp.show();
    }*/
    
  
    
    
}
