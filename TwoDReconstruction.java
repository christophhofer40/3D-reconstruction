/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


/**
 *
 * @author christoph
 */

import ij.*;
import ij.gui.*;
import ij.plugin.*;
import ij.ImageStack.*;
import ij.WindowManager.*;
import ij.io.FileInfo;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;
import java.io.FileReader;
import java.io.File;
import java.io.FileWriter;          
import java.io.*;
import java.util.Collections;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Locale;
import java.nio.file.*; 
import java.nio.file.Files; 
//import christoph.Master;
import christoph.Fakesimulation;
import christoph.QSTEM_simulator;
import christoph.Quaternion;
import christoph.Shear;
import christoph.ImageCalculations;
import christoph.Optimization;
//import christoph.JythonBEAM;
import ij.plugin.filter.MaximumFinder;
import java.awt.Color;
import java.awt.Frame;
import java.net.InetAddress;
import javax.swing.SwingConstants;


public class TwoDReconstruction implements PlugIn 
{
    
    static ImageStack stack;
    static ImagePlus imp;
    static ImagePlus masterp;
    static ImageStack masterSt;
    static int impWidth;
    static int impHeight;
    static int stackNb;
    static int viewcount_one=0;
    static String impT= "title";
    static ImageStack stack_model;
    static ImageStack stack_sim;
    static ImagePlus stackp;
    static ImagePlus stacksimp;
    static ImageStack modelSt;
    static ImagePlus model;
    static String modelT= "modeltitle";
    static ImageStack mergedSt;
    static ImagePlus merged;
    static ImageStack diffSt;
    static ImagePlus diffp;
    static ImageStack beamstack;
    static ImagePlus beamplus;
    
    static ArrayList<Atoms> masters=new ArrayList<Atoms>(10);
    static ArrayList<Bond> bonds=new ArrayList<Bond>(3000);
    float merit;
    
    
    static float[] masterarray;
    static float[] simarray;

    //float fov=3f;
    float stepsize=1;
    float threshold=30;
    int cropradius=20;
    float smoothweight=2f;
    float smoothscale=0.5f;
    float factor=0;
    String präfix;
    int slice=0;
   
    int globalcounter=0;
    int automatic_rounds=0;
    int iterations=10;
    int iterations2=1000;
    int totalrounds=1;
   

    int best_view=0;
    float solid_diameter=6;
    long starttime=System.currentTimeMillis();

    
    FileInfo fi;
    String path;
    
    String topFileStr;
    File topFile;  
    String angFileStr;
    File angFile;
    File datafile;
    static Random random = new Random(System.currentTimeMillis());
    
    //static ArrayList<Ring> rings=new ArrayList<Ring>();
    
    static boolean reset=true;
    static boolean initmodel=false;
    static boolean move_atoms=false;
    static boolean tune_quaternion=false;
    static boolean fix_tune_quaternion=false;
    static boolean fixtilt_tune=false;
    static boolean badd_remove_atoms=false;
    static boolean tune_skewmatrix=false;
    static boolean smooth=false;
    static boolean smooth2=false;
    static boolean opt_fov=false;
    static boolean xyz=false;
    static boolean orthmove;
    static boolean rdmmove;
    static boolean bluropt;
    static boolean abopt;
    static boolean dialog;
    static boolean switch_=false;
    static boolean show_sims=true;
    static boolean show_moves=true;
    static boolean qstem=false;
    static boolean automatic=false;
    static boolean edgeatoms=false;
    static boolean analyze=false;
    static boolean optfov=false;
    static boolean adjust_intensities=false;
    static boolean adjust_global_intensities=false;
    static boolean uncorr_views=true;
    static boolean read_beamparam=false;
    static boolean save_topfile=false;
    
    
    @Override
    public void run(String arg) 
    {
        Locale.setDefault(Locale.UK);
        while (DoDialog())
        {
        if(validateInput()!=0)
            {continue;}
            IJ.resetEscape(); 
            run_slice();
        }
    }
    
    
    public boolean DoDialog()
    {
        int[] idArray = WindowManager.getIDList();
        int idlen = 0;
        if(idArray != null)
        {	idlen = idArray.length;}
        String[] imptitleArray = new String[idlen+1];
        imptitleArray[idlen] = "<current>";
        for (int i = 0; i < idlen; ++i)
        {	
            String title = WindowManager.getImage(idArray[i]).getTitle();
            imptitleArray[i] = title;
        }
        
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("2D-Reconstructor");
        gd.setSmartRecording(true);


        gd.addChoice("raw image (GRAY32)", imptitleArray, impT);
        //gd.addChoice("2 bit model", imptitleArray, modelT);
		
        gd.addMessage("Initialization");
        gd.addCheckbox("reset model", false);
        gd.addCheckbox("read beamparamters", false);
        
        gd.addCheckbox("\t init from imagetitle 'model.tif'", false);
        
        gd.addMessage("Move atoms");
        gd.addCheckbox("fine tune atoms", rdmmove);
        gd.addNumericField("\t number of iterations", iterations, 0);
        //gd.addCheckbox("add/remove atoms", badd_remove_atoms);
        //gd.addNumericField("solid diameter (px)", solid_diameter, 1);
        //gd.addNumericField("\t number of iterations", iterations2, 0);
        //gd.addCheckbox("smooth rings", smooth2);
        //gd.addNumericField("\tweight of the smooth (between 0 and 1 recommended)", smoothscale, 2);
        gd.addMessage("Adjust beam parameters");
        gd.addCheckbox("match Gaussian blur", bluropt);
        gd.addCheckbox("match aberrations", abopt);
        gd.addCheckbox("match fov", optfov);
        
        gd.addMessage("Adjust intensities");
        gd.addCheckbox("adjust atom intensities", adjust_intensities);
        gd.addCheckbox("adjust global atom intensities", adjust_global_intensities);
        //gd.addCheckbox("uncorrelated_views", uncorr_views);
        //gd.addCheckbox("write xyz", xyz);
        gd.addNumericField("\t mean stepsize (sigma of normaldistribution in pixel)", stepsize, 2);
        //gd.addCheckbox("\t use qstem?", qstem);
        //gd.addCheckbox("if not: show simulations? (slower, but cooler)", show_sims);
        //gd.addCheckbox("show all moves? (slower, but cooler)", show_moves);
       
        gd.addNumericField(" totalrounds (-1...forever)", totalrounds, 0); 
        gd.addNumericField("slice to optimize (0...all slices)", slice, 0);
        gd.addCheckbox("save topology file", save_topfile);
        gd.showDialog();

        impT=gd.getNextChoice();
        //modelT=gd.getNextChoice();
        
        reset=gd.getNextBoolean();
        read_beamparam=gd.getNextBoolean();
        initmodel=gd.getNextBoolean();
        rdmmove=gd.getNextBoolean();
        iterations=(int)gd.getNextNumber();
        //badd_remove_atoms=gd.getNextBoolean();
        //solid_diameter=(float)gd.getNextNumber();
        //iterations2=(int)gd.getNextNumber();
        //smooth2=gd.getNextBoolean();
        //smoothscale=(float)gd.getNextNumber();       
        bluropt=gd.getNextBoolean();
        abopt=gd.getNextBoolean();
        optfov=gd.getNextBoolean();
        adjust_intensities=gd.getNextBoolean();
        adjust_global_intensities=gd.getNextBoolean();
        //uncorr_views=gd.getNextBoolean();
        
        stepsize=(float) gd.getNextNumber();
        //qstem=gd.getNextBoolean();
        //show_sims=gd.getNextBoolean();
        //show_moves=gd.getNextBoolean();      
        
        totalrounds=(int)gd.getNextNumber();
        slice=(int)gd.getNextNumber();
        save_topfile=gd.getNextBoolean();
        if(gd.wasCanceled())
        {
            return false;
        }
        return true;
    }

    int validateInput()
    {
        if(impT.equals("<current>"))
        {	
            imp = IJ.getImage();
            if(imp != null)
            {	impT = imp.getTitle();}	
        }
        
        else
        { 
            imp  = WindowManager.getImage(impT);
            //model = WindowManager.getImage(modelT);
        }
        impT=impT.substring(0, impT.length()-4);
        fi = imp.getOriginalFileInfo();
        path = fi.directory;
        topFileStr = path+impT+".top";
        topFile= new File(topFileStr);
        imp=imp.duplicate();
        //model=model.duplicate();
        stack=imp.getStack();
        stackNb=stack.getSize();
        impWidth=stack.getWidth();
        impHeight=stack.getHeight();
        //modelSt=model.getStack();
        
        
        
        if(imp == null)
        {
            IJ.log("Error: could not find " + impT);
            return 2; 
        }
        return 0;

    }
    void run_slice()
    {
        IJ.resetEscape();
        if(reset)
        {
            if (initmodel)
                init_from_model();
            else
            {
                try 
                {
                    read_top(topFile);
                    IJ.log("atom structure initialized");

                } catch (IOException ex) 
                {
                    Logger.getLogger(TwoDReconstruction.class.getName()).log(Level.SEVERE, null, ex);
                    IJ.log("something went wrong with reading the topfile!");
                }
                /*if (master.views.size()==1)
                {
                    for (Master.Atom atom: master.atoms)
                    {
                        atom.viewcounts=2;
                    }
                    master.seen_atoms=master.atoms.size();
                }*/
            }
            File beamfile=new File(path+"beamparameters.txt");
            if (!beamfile.exists())
                 write_beam_param();
            try {
                read_beam_param(beamfile);
                IJ.log("aberrations setted from "+path+"beamparameters.txt");
                
            } catch (FileNotFoundException ex) {
                Logger.getLogger(TwoDReconstruction.class.getName()).log(Level.SEVERE, null, ex);
            }
            show_views();
            //change_fov(fov);
        }   
        else
        {
        	if(read_beamparam)
        		try {
					read_beam_param(new File(path+"beamparameters.txt"));
					for (int i=0;i<masters.size();++i)
                    {
						beamstack.setPixels((float[])masters.get(i).sim.beam.beamfht.getPixels(),i+1);
                    }
					IJ.log("aberrations setted from "+path+"beamparameters.txt");
					
				} catch (FileNotFoundException ex) {
					Logger.getLogger(TwoDReconstruction.class.getName()).log(Level.SEVERE, null, ex);
				}

        }
        check_imagestatus();
        int it=0;
        update();
        IJ.log("derviation before: "+merit );
 
        do
        {
            it+=1;
            for(int k=0;k< iterations;k++)
            {
                if(rdmmove)
                random_move(); 
                if(IJ.escapePressed())
                    break;
            }
            for(int k=0;k< iterations2;k++)
            {
                if(badd_remove_atoms)
                    add_and_remove_atoms();
            }
                
            update_merit();
            
            if(abopt)
                optimize_aberrations2();
            
            if(bluropt)
                optimize_gaussianblur();
            if(optfov)
            {
            	for(int i=0;i<stackNb;++i)
                	optimize_fov(i);
                write_beam_param();
            }
            if(adjust_intensities)
            	adjust_atomintensities();
            if(adjust_global_intensities)
            	adjust_global_atomintensities(true);
            
            
            update_merit();
            write_top();
            //save_tifdata2(path);
            IJ.log("deviation after: "+merit );
            //check_imagestatus();
            stackp.updateAndRepaintWindow();  
            stacksimp.updateAndRepaintWindow();
            diffp.updateAndRepaintWindow();
            if(save_topfile)
                try 
                {
                    save_top();
                }
                catch (Exception ex) {
					Logger.getLogger(TwoDReconstruction.class.getName()).log(Level.SEVERE, null, ex);
				}
            if(IJ.escapePressed())
                break;
        }while(it!=totalrounds); 
    }
    

    void init_from_model()
    {
        ImagePlus modelp = WindowManager.getImage("model.tif");
        ImageStack modelst=modelp.getImageStack();
        int width=modelp.getWidth();    
        for (int i=1;i<=modelst.getSize();++i)
        {
            int id=0;
            masters.add(new Atoms(i-1));
            float[] pxarray=(float[])modelst.getPixels(i);
            
            for (int j=0;j<pxarray.length;++j)
            {
                float val=pxarray[j];
                if (val==0)
                    continue;
                float x=(int)(j%width);
                float y=(int)(j/width);
                masters.get(i-1).atoms.add(new Atom(x,y,id,i-1));
                id+=1;
            }
            if (id==0)
                masters.get(i-1).atoms.add(new Atom(width/2,width/2,id,i-1));

        }
        
    }

    
    void read_top(File topFile)
        throws FileNotFoundException, IOException 
    {   
        masters.clear();
        if (topFile == null || !topFile.canRead()) 
        {
            throw new IllegalArgumentException("file not readable: " + topFile);
        }
        FileReader topReader=new FileReader(topFile);
        final Scanner s= new Scanner(topFile).useDelimiter("\t");
        int chapter=-1;
        boolean bcorrelated=false;
        ArrayList<Integer> znumber=new ArrayList<Integer>();
        while(s.hasNextLine())
        {    
            String line = s.nextLine();
            Scanner ls = new Scanner(line);
            if(line.startsWith("#") || line.equals(""))
            {
                continue;
            }
            
            if(line.startsWith("UNCORRELATED"))
            {
                uncorr_views=true;
                continue;
            }
            
            if(line.startsWith("VIEW"))
            {
                chapter+=1;
                masters.add(new Atoms(chapter));
                continue;
            }
            if(line.startsWith("BACKGROUND"))
            {
            	ls.next();
            	masters.get(chapter).backgroundlevel=ls.nextFloat();
            }

            if(chapter==-1)
            {

                //Scanner ls = new Scanner(line);
                Scanner ls2 = new Scanner(line);
                if (line.startsWith("ATOM"))
                {
                    int element=6;
                    for(int i=0; i<6;++i)
                    {
                        ls.next();
                        ls2.next();
                    }
                    if (ls.hasNextInt())
                    {
                        element=ls.nextInt();
                        if (element!=6)
                        {
                        	bcorrelated=true;
                        	System.out.println("Foreigner element with atomnumber "+element);
                        }
                        	
                    }
                    
                    znumber.add(element);                       
                }
                if(line.startsWith("BOND"))
                {
                    //masters.add(new Atoms(chapter));
                    continue;
                }
                continue;
            }
            if(line.startsWith("ATOM"))
            {
                //Scanner ls = new Scanner(line);
                ls.next();
                int id=ls.nextInt();
                float x=ls.nextFloat();
                float y=ls.nextFloat(); 
                masters.get(chapter).atoms.add(new Atom(x,y,id,chapter));
                
                if(id>-1 && bcorrelated)
                        masters.get(chapter).atoms.get(masters.get(chapter).atoms.size()-1).element=znumber.get(id);
                else
                    masters.get(chapter).atoms.get(masters.get(chapter).atoms.size()-1).element=6;
                int elem=masters.get(chapter).atoms.get(masters.get(chapter).atoms.size()-1).element;
                while (ls.hasNext())
                {
                        String kwarg=ls.next();
                        if(kwarg.startsWith("#intensity"))
                        {
                                float intensity=ls.nextFloat();
                                //System.out.println(intensity);
                                masters.get(chapter).atoms.get(masters.get(chapter).atoms.size()-1).intensity_factor=intensity;
                        }
                }
            }
        }     
    }
    
    void save_top() throws FileNotFoundException, IOException 
    {
        FileReader topReader=new FileReader(topFile);
        String src = path+impT+".top";
        File destFile=new File(path+impT+"_cp.top");
        if(destFile.exists()) 
        {
            destFile.delete();
        }
        Path temp=Files.move(Paths.get(src), Paths.get(path+impT+"_cp.top"));
        File topFile=new File(path+impT+"_cp.top");
        final Scanner s= new Scanner(topFile).useDelimiter("\t");
        File file = new File(src);
        FileWriter fw = new FileWriter(file, false);
        PrintWriter pw = new PrintWriter(fw);
        while(s.hasNextLine())
        {    
            String line = s.nextLine();
            //Scanner ls = new Scanner(line);
            
            
            
            if(line.startsWith("VIEW"))
            {
                break;
            }
            pw.println(line);

            

        }
            
        for(int v = 0; v<stackNb; ++v)
        {
            if( (masters.get(v).atoms.size() == 0))
            {	continue;}
            pw.println("\n\nVIEW\t" + v);
            pw.println("QUATERNION\t0\t0");
            pw.println("SKEWMATRIX\t1\t0\t0\t1");
            pw.println();
            pw.println();

            for(Atom atom : masters.get(v).atoms)
            {
                if (!Float.isNaN(atom.x)&&!Float.isNaN(atom.y))
                    pw.println("ATOM\t"+atom.id+"\t"+atom.x+"\t"+atom.y+
                            "\t#intensity:\t"+(atom.intensity_factor));//intensity factor != integrated intensity
            }	
        }
            
    }

    void read_beam_param(File f) throws FileNotFoundException
    {
           
        if (f == null || !f.canRead()) 
        {
            throw new IllegalArgumentException("file not readable: " + f);
        }
        final Scanner s= new Scanner(f).useDelimiter("\t");  
        int slice=0;
        while(s.hasNextLine())
        {
            String line = s.nextLine();
            Scanner ls = new Scanner(line);
            if (line.startsWith("#"))
                continue;
            if(line.startsWith("FieldOfView"))
            {
                ls.next();
                masters.get(slice).fov=ls.nextFloat();
                continue;
            }    
            if(line.startsWith("VIEW"))
            {
                ls.next();
                slice=ls.nextInt();
                continue;
            }
            if(line.startsWith("Sigma"))
            {
                ls.next();
                float sig=ls.nextFloat(); 
                masters.get(slice).sim.sigma=sig;                          
                continue;
            }
			if(line.startsWith("Aberrations"))
            {
				ls.next();
				for(int k=0;k<7;++k)
				{
					float ab=0;
					if (ls.hasNext())
					{
						ab=ls.nextFloat();
					}
					masters.get(slice).sim.beam.aberrations[k]=ab;  			
				}
				
			}
			
            
        }
        for(int sl=0;sl<masters.size();sl++)
        {
        	masters.get(sl).sim.fov=masters.get(sl).fov;
        	masters.get(sl).sim.beam.initBeam();
        }
       
        
            
    }
    
    void write_beam_param()
    {
        try
        {
            File file = new File(path+"beamparameters.txt");
            FileWriter fw = new FileWriter(file, false);
            PrintWriter pw = new PrintWriter(fw);
            pw.println("# reads out aberrations and sigma for gaussblur (for ImageJ plugin)");
            pw.println("# Slicenumber	C10	C12a	C12b	C21a	C21b	C23a	C23b");
            for( int i=0; i<masters.size();i++)
            {
            	pw.println("VIEW:\t"+i);
                pw.println("Aberrations:\t"+masters.get(i).sim.beam.aberrations[0]+"\t"
                                +masters.get(i).sim.beam.aberrations[1]+"\t"
                                +masters.get(i).sim.beam.aberrations[2]+"\t"
                                +masters.get(i).sim.beam.aberrations[3]+"\t"
                                +masters.get(i).sim.beam.aberrations[4]+"\t"
                                +masters.get(i).sim.beam.aberrations[5]+"\t"
                                +masters.get(i).sim.beam.aberrations[6]+"\t");
                pw.println("Sigma:\t"+ +masters.get(i).sim.sigma);
                pw.println("FieldOfView:\t"+masters.get(i).fov);
                beamstack.setPixels((float[])masters.get(i).sim.beam.beamfht.getPixels(),i+1);
            }
            pw.close();
            IJ.log("rewrote " + file.getName());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
        beamplus.updateAndRepaintWindow();
    }

    /*
    void rotate_and_translate_every_atom()
    {
        
        for(Master.View view : master.views)
            for(Master.Atom2D atom2d : view.atoms2D)
            {
                atom2d.update_rotation();
            }
        translate_atoms();
    }
    */

    void write_top()
    {
        try
        {
            File topfile ;
            topfile= new File(path+imp.getShortTitle() + "_new.top");

            FileWriter fw = new FileWriter(topfile, false);
            PrintWriter pw = new PrintWriter(fw);

            pw.println("#creator: " + getClass().getSimpleName() + " source: " + imp.getTitle());
            pw.println("#2D optimized Top file");
            pw.println("#only 2D positions in the views optimized");
            pw.println("#copy and paste to the original top file");
            pw.println("MASTER");
            
            for (Atom atom: masters.get(0).atoms)
            {
            	if (atom.id<0)
            		continue;
            	pw.println("ATOM\t"+atom.id+"\t1\t"+atom.x+"\t"+atom.y+"\t0\t"+
                        atom.element);
            }
            pw.print("\n\n");
            for(int v = 0; v<stackNb; ++v)
            {
                if( (masters.get(v).atoms.size() == 0))
                {	continue;}
                pw.println("\n\nVIEW\t" + v);
                pw.println("BACKGROUND\t"+masters.get(v).backgroundlevel);
                pw.println();
                pw.println();

                for(Atom atom : masters.get(v).atoms)
                {
                    if (!Float.isNaN(atom.x)&&!Float.isNaN(atom.y))
                        pw.println("ATOM\t"+atom.id+"\t"+atom.x+"\t"+atom.y+
                                "\t#intensity:\t"+(atom.intensity_factor));//intensity factor != integrated intensity
                }	
            }
            pw.close();
            IJ.log("rewrote " + topfile.getName());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    
    }
    void check_imagestatus()
    {
        stackp=WindowManager.getImage("model_of_views");
        stacksimp=WindowManager.getImage("simulation_of_views");
        diffp=WindowManager.getImage("difference");
        beamplus=WindowManager.getImage("beamprobe");
        if(stackp==null)
        {
            stackp=new ImagePlus("model_of_views",stack_model);
            stackp.setStack(stack_model);
            stackp.show();
        }
        else
        	stackp.setStack(stack_model);
        
        if(stacksimp==null)
        {
            stacksimp=new ImagePlus("simulation_of_views",stack_sim);
            stacksimp.setStack(stack_sim);
            stacksimp.show();
        }
        else
        	stacksimp.setStack(stack_sim);
        
        if(diffp==null)
        {
            diffp=new ImagePlus("difference",diffSt);
            diffp.setStack(diffSt);
            diffp.show();
        }
        else
        	diffp.setStack(diffSt);
        
        if(beamplus==null)
        {
            beamplus=new ImagePlus("beamprobe",beamstack);
            beamplus.setStack(beamstack);
            beamplus.show();
        }
        else
        	beamplus.setStack(beamstack);
        
        
        //merged.show();
        
    }
   
    void show_views()
    { 
        stack_model=new ImageStack(impWidth,impHeight,stackNb);
        stack_sim=new ImageStack(impWidth,impHeight,stackNb);
        diffSt=new ImageStack(impWidth,impHeight,stackNb);
        beamstack=new ImageStack(impWidth,impHeight,stackNb);
        //update();
        for(int i=0;i<masters.size();++i)
        {
        	beamstack.setPixels((float[])masters.get(i).sim.beam.beamfht.getPixels(),i+1);
            masters.get(i).update_pixels();
			
        }

        check_imagestatus();
        stackp.updateAndRepaintWindow();
        stacksimp.updateAndRepaintWindow();
        diffp.updateAndRepaintWindow();
        beamplus.updateAndRepaintWindow();
        starttime=System.currentTimeMillis();        
    }
    
    void update_merit()
    {
        merit=0;
        for (Atoms atoms:masters)
        {
            merit+=atoms.difference;
        }
        merit=merit/masters.size();
    }
    
    void update()
    {
    	int counter=0;
        for(Atoms a:masters)
        {
        	counter+=1;
            if (slice!=counter && slice>0)
                continue;
            a.update_pixels();
            
        }
        update_merit();
    }
    
    void doProcess()
    {
       
    }
    
    boolean in_triangle(float[] p1, float[] p2, float[] p3, float[] m)
    {
        float[] v1={p2[0]-p1[0],p2[1]-p1[1]};
        float[] v2={p3[0]-p2[0],p3[1]-p2[1]};
        float[] v3={p1[0]-p3[0],p1[1]-p3[1]};
        
        float a = (p2[0]-p1[0])*(m[1]-p1[1])-(p2[1]-p1[1])*(m[0]-p1[0]);
        float b = (p3[0]-p2[0])*(m[1]-p2[1])-(p3[1]-p2[1])*(m[0]-p2[0]);
        float c = (p1[0]-p3[0])*(m[1]-p3[1])-(p1[1]-p3[1])*(m[0]-p3[0]);
        if((a>=0 && b>=0 && c>=0)||(a<=0 && b<=0 && c<=0))
            return true;
        else
            return false;
    }
  
    void random_move()
    {
        float[] shift = new float[2];
        int counter=0;
        for (Atoms atoms:masters)
        {
            counter+=1;
            if (slice!=counter && slice>0)
                continue;
            for(Atom atom :atoms.atoms)
            {
                if(IJ.escapePressed())
                    break;
                if(Float.isNaN(atoms.expIm[((int)(atom.x+0.5)+(int)(atom.y+0.5)*impWidth)]))
            	{
            		continue;
            	}
                float norm=0;
                do
                {
                    shift[0]=random.nextFloat()-0.5f;
                    shift[1]=random.nextFloat()-0.5f;
                    norm=shift[0]*shift[0]+shift[1]*shift[1];
                }while(norm>0.25f || norm==0.0f);
                shift[0]=shift[0]/norm;
                shift[1]=shift[1]/norm;
                optimize_atom(shift,atom);
            }
            
        }
    }
    
    void add_and_remove_atoms()
    {
        int counter=0;
        for (Atoms atoms:masters)
        {
            float olddiff=atoms.difference;
            counter+=1;
            if (slice!=counter && slice!=0)
                continue;
            int offset=random.nextInt(5);
            if (offset==5)
                offset+=random.nextInt(5);
            int maxind=ImageCalculations.get_maxarg(atoms.diffarray,offset);
            
            maxind=(int)(random.nextFloat()*imp.getWidth()*imp.getWidth());
            int nancounter=0;
            while (Float.isNaN(atoms.expIm[maxind])&&nancounter<100)
            {
                maxind=(int)(random.nextFloat()*imp.getWidth()*imp.getWidth());
                nancounter+=1;
            }
                
            int id=atoms.atoms.size();
            int x=(int)(maxind%imp.getWidth());
            int y=(int)(maxind/imp.getWidth());
            Atom atom=new Atom(x,y,id,counter-1);        
            int replace=check_dist(atom); 
            float rx=0;
            float ry=0;
            Atom replaceat=null;
            if (replace>=0)
            {
                replaceat=atoms.atoms.get(replace);
                rx=replaceat.x;
                ry=replaceat.y;            
                replaceat.x=x;
                replaceat.y=y;
            }
            else
            {
                atoms.atoms.add(atom);
            }
            atoms.update_pixels();
            
            if (olddiff<=atoms.difference)
            {
                if (replace>=0)
                {
                    replaceat.x=rx;
                    replaceat.y=ry;
                }
                else
                {
                    atoms.atoms.remove(id);
                }
                    
                atoms.update_pixels();
            }
            else
                System.out.println("new diff:\t"+atoms.difference);
        } 
    }
    
    
    //>=0: new atom replaces this id
    //-1: new atom added in model
    int check_dist(Atom atom)
    {
        int smallestid=get_closest_atom(atom);
        if (masters.size()==0)
        {
            System.out.println("Master is empty?");
        }
        float mindist=get_dist(atom, masters.get(atom.view).atoms.get(smallestid));
        if (mindist<solid_diameter*solid_diameter)
            return smallestid;
        else
            return -1;
    }
    
    int get_closest_atom(Atom atom)
    {
        int closest=0;
        float mindist=-1;
        for (Atom atomi: masters.get(atom.view).atoms)
        {
            if (atom==atomi)
                continue;
            float dist= (atomi.x-atom.x)*(atomi.x-atom.x)+(atomi.y-atom.y)*(atomi.y-atom.y);
            if(dist<mindist || mindist==-1)
            {
                mindist=dist;
                closest=atomi.id;           
            }
            
        }
        return closest;
    }
    float get_dist(Atom atomi, Atom atom)
    {
        float dist= (atomi.x-atom.x)*(atomi.x-atom.x)+(atomi.y-atom.y)*(atomi.y-atom.y);
        return dist;
    }
    
    
    void optimize_atom(float[] vector, Atom atom)
    {
        float[] vector2=atom.force.clone();
        if(Float.isNaN(vector[0]))
        {
            vector[0]=vector2[0];
            vector[1]=vector2[1];
        }
        float oldx=atom.x;
        float oldy=atom.y;
        float olddiff=masters.get(atom.view).difference; 
        float inorm = 1/(float)Math.sqrt(vector[0]*vector[0]+
        vector[1]*vector[1]);
        float[] shift= new float[2];
        float length;
        do
        {
            length = (float) (random.nextGaussian()*stepsize);
        }while (length==0);
        shift[0]=(vector[0]*inorm)*length;
        shift[1]=(vector[1]*inorm)*length;
        
        float[] p1=new float[3];
        float[] p2= new float[3];
        float p1grad;
        float p2grad;
        float norm;
        //System.out.println("pos before: "+ atom.x+ " "+ atom.y+ " "+ atom.z);
        atom.shift_coordinates( shift);
        //System.out.println("after:  "+atom.sum_distances());
        p1grad=(float)(masters.get(atom.view).difference-olddiff);
        if (p1grad<0)
        {         
            System.out.println("new merit1: "+merit);
            push_pixels();
            atom.force[0]=shift[0];
            atom.force[1]=shift[1];
            return;
        }
        p1[0]=atom.x;
        p1[1]=atom.y;
        
        shift[0]=-2*shift[0];
        shift[1]=-2*shift[1];
        atom.shift_coordinates(shift);    
        p2grad=(float)(masters.get(atom.view).difference-olddiff);
        if (p2grad<0)
        {
            System.out.println("new merit2: "+merit);
            push_pixels();
            atom.force[0]=0.5f*shift[0];
            atom.force[1]=0.5f*shift[1];
            return;
        }
        p2[0]=atom.x;
        p2[1]=atom.y;       

        atom.set_coordinates(oldx,oldy);
        Optimization opt = new Optimization(atom.x,atom.y,0,p1,p1grad,p2,p2grad);
        shift=opt.get_shift(0);
        norm=shift[0]*shift[0]+shift[1]*shift[1];
        //System.out.println("p1: "+Arrays.toString(p1));
        //System.out.println("p2: "+Arrays.toString(p2));
        //System.out.println(Arrays.toString(shift));
        
        atom.shift_coordinates( shift);
        //System.out.println("pos after: "+ atom.x+ " "+ atom.y+ " "+ atom.z);
        if((masters.get(atom.view).difference>olddiff))
        {
            System.out.println("bad merit : "+merit);
            atom.set_coordinates( oldx,oldy);
            //System.out.println("correlation after reset: "+master.correlation);
            //System.out.println("energy after reset:      "+master.energy);
            //master.merit=oldmerit;
            //master.correlation=oldcorr;
            //master.energy=oldenergy;

        }
        else
        {
            System.out.println("new merit : "+merit);

        }
                
        


        atom.force[0]=vector2[0];
        atom.force[1]=vector2[1];
        atom.force[2]=vector2[2];
        push_pixels();
        //System.out.println("before: "+ olddev);
        //System.out.println("after: "+ master.deviation);
    }
    
  
    
    void optimize_aberrations2()//only one round
    {
        //master.update_deviations();
        IJ.log("optimizing aberrations...");
        float[] step = new float[7];
        double match;
        ArrayList<Boolean> ind=new ArrayList<Boolean>();
        
        for(int i=0;i<7;++i)
        {
            ind.add(false);
        }
        int counter=0;
        for (Atoms atoms:masters)
        {
            counter+=1;
            if (slice!=counter && slice>0)
            continue;
            for(int i=0;i<7;++i)
            {
                ind.set(i,false);
            }
            for (int k=0 ; k<3;++k)
            {
                step[k]= ((float)random.nextGaussian());
            }
            for (int k=3 ; k<7;++k)
            {
                step[k]=((float)random.nextGaussian()) * 50.0f;
            }
            int s=0;
            float olddif;
            do
            {

                for(int i=ind.size()-1;i>=0;--i)
                {
                    if(ind.get(i)==true)
                        continue;

                    for(int r=0;r<2;++r)
                    {
                        olddif=atoms.difference;
                        atoms.sim.beam.aberrations[i]+=step[i];
                        atoms.sim.beam.initBeam();
                        atoms.simarray=ImageCalculations.normalize_image(atoms.simulate_image(qstem),atoms.expIm);
                        //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                        atoms.calc_difference();                        //System.out.println(newmerit);
                        if(show_moves)
                            push_pixels();
                        if (atoms.difference<olddif)
                        {
                            ind.set(i, true);
                            break;
                        }
                        else
                        {
                            atoms.sim.beam.aberrations[i]-=step[i];
                            atoms.sim.beam.initBeam();
                            atoms.simarray=ImageCalculations.normalize_image(atoms.simulate_image(qstem),atoms.expIm);
                            //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                            atoms.calc_difference();     
                            step[i]*=-1;
                        }
                    }
                    step[i]*=0.5;
                        
                }
                
                s+=1;

            }while ( s<=1);
        }
        for (Atoms atoms:masters)
        {
            atoms.sim.beam.initBeam();
            atoms.simarray=ImageCalculations.normalize_image(atoms.simulate_image(qstem),atoms.expIm);
            atoms.update_pixels();
        }
        update_merit();
        write_beam_param();
    }
    
    void optimize_fov(int slice)
    {
        double oldcorr=merit;
        System.out.println("corr with fov old fov: "+ oldcorr);
        float oldfov=masters.get(slice).fov;
        int sign=random.nextBoolean() ? 1:-1;
        float step=random.nextFloat()*sign;    
        for(int i=0; i<8;i++)
        {
        	change_fov(oldfov+step,slice);  	
            if(merit<oldcorr)
            {
                System.out.println("new fov: "+(masters.get(slice).fov+step));
                break;
            }
            else
            {
            	change_fov(oldfov,slice);
            	step*=-0.7;
            }
        }
        
    }
    
    public void change_fov(float fov,int slice)
    {
        masters.get(slice).fov=fov;
        masters.get(slice).sim.fov=fov;
        masters.get(slice).sim.beam.initBeam();
        update();
        
    }
    
    public void adjust_atomintensities()
    {
    	int counter=0;
    	for (Atoms atoms:masters)
        {
            counter+=1;
            if (slice!=counter && slice>0)
                continue;
            float olddiff=atoms.difference;
            float mean=0;
            for(Atom atom :atoms.atoms)
            {
            	float empiricfactor=(float) Math.pow(atom.element,1.6);
            	float intensity=empiricfactor*atom.intensity_factor;
            	mean+=intensity;
            }
            mean/=atoms.atoms.size();
            
            for(Atom atom :atoms.atoms)
            {
            	//if (atom.id!=95)
            		//continue;
            	if(Float.isNaN(atoms.expIm[((int)(atom.x+0.5)+(int)(atom.y+0.5)*impWidth)]))
            	{
            		continue;
            	}
            		
				int sign=random.nextBoolean() ? 1:-1;
				float step=(random.nextFloat()*mean/10f+mean/200)*sign;		
				float empiricfactor=(float) Math.pow(atom.element,1.6);
                if(IJ.escapePressed())
                    break;
                for (int it=0;it<2;++it)
                {
                	float oldint=atom.intensity_factor;
					atom.intensity_factor=(atom.intensity_factor*empiricfactor+step)/empiricfactor;
					//System.out.println( atom.intensity_factor);
					atoms.update_pixels();
					if (atoms.difference<olddiff)
					{
						olddiff=atoms.difference;
						System.out.println("adjusting atom intensities:\tnew difference:\t"+olddiff);
	
					}
					else
					{
						atom.intensity_factor=oldint;
						step*=-1;
					}
					
				}
            }
            
        }
    	update();
    }
    
    public void adjust_global_atomintensities(boolean background)//True: also adjust background
    {
    	int counter=0;
    	for (Atoms atoms:masters)
        {
            counter+=1;
            if (slice!=counter && slice>0)
                continue;
            float olddiff=atoms.difference;
            float step=0.05f;
            float stepbg=0.002f;
            while(Math.abs(step)>0.000001)
            {
                if(IJ.escapePressed())
                    break;
                for(Atom atom:atoms.atoms)
                	atom.intensity_factor*=(step+1);
                //System.out.println( atom.intensity_factor);
                atoms.update_pixels();
                if (atoms.difference<olddiff)
                {
                    olddiff=atoms.difference;
                    step*=1.1;
                    System.out.println("adjusting global intensities:\tnew difference:\t"+olddiff);


                }
                else
                {
                    for(Atom atom:atoms.atoms)
                            atom.intensity_factor/=(step+1);
                    step*=-0.7;
                }
                
                //adjust background
                if (background)
                {
                    float oldstepbg=stepbg;
                    float oldbglevel=atoms.backgroundlevel;
                    atoms.backgroundlevel+=stepbg;
                    if (atoms.backgroundlevel<0)
                            atoms.backgroundlevel=0;
                    atoms.update_pixels();
                    if (atoms.difference<olddiff)
                    {
                            olddiff=atoms.difference;
                            stepbg*=1.1;
                            System.out.println("adjusting backgroundlevel:\tnew difference:\t"+olddiff);
                    }
                    else
                    {
                            atoms.backgroundlevel=oldbglevel;
                            stepbg*=-0.7;
                    }
                }
                
            }
            
        }
    	update();
    }


    void push_pixels()
    {
        
    }
    
    void optimize_gaussianblur()
    {
        //master.update_deviations();
        IJ.log("optimizing aberrations...");
        float step=(random.nextFloat()-0.5f)*2;
        double match;
        int counter=0;
        for (Atoms atoms:masters)
        {
        	counter+=1;
        	if (slice!=counter && slice>0)
                continue;
            
            int s=0;
            float olddif;
            do
            {
                olddif=atoms.difference;
                atoms.sim.sigma+=step;
                atoms.simarray=ImageCalculations.normalize_image(atoms.simulate_image(qstem),atoms.expIm);
                //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                atoms.calc_difference();                        //System.out.println(newmerit);
                if(show_moves)
                        push_pixels();
                if (atoms.difference<olddif)
                        break;

                else
                {
                        atoms.sim.sigma-=step;
                        atoms.simarray=ImageCalculations.normalize_image(atoms.simulate_image(qstem),atoms.expIm);
                        //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                        atoms.calc_difference();     
                        step*=-0.7;
                }                        
                
                
                s+=1;

            }while ( s<=3);
        }

        update();
        write_beam_param();
    }
  
    
    /*void smooth_structure2()
    {
        float[] vector;
        float x,y;
        int counter;
        for(Atoms matoms:masters)
        {
            matoms.update_rings();
            for(Atom atom: matoms.atoms)
            {
               // System.out.println("\n\nbefore: "+atom.toString());

                vector=new float[3];
                x=0;
                y=0;
                counter=0;
                for(Master.Ring ring: atom.rings)
                {
                    //System.out.println(ring.toString());
                    vector=modify_length(get_vector(ring.center,new float[] {atom.x,atom.y,0}),ring.radius);
                    //System.out.println(Arrays.toString(vector)+ "... center - projection");
                    x+=ring.center[0]+(vector[0]);
                    y+=ring.center[1]+(vector[1]);
                    counter+=1;
                    //System.out.println(ring.center[0]+(vector[0])+" "+ring.center[1]+(vector[1])+" "+ring.center[2]+(vector[2])+" [xyz] projection");

                }
                if (counter==0)
                    continue;
                x=x/counter;
                y=y/counter;
                //System.out.println(x+" "+y+" "+z+" [xyz] targetpoint");
                //System.out.println(Arrays.toString(new float[]{(x-atom.x)*scale,(y-atom.y)*scale,(z-atom.z)*scale}));
                atom.shift_coordinates(new float[]{(x-atom.x)*smoothscale,(y-atom.y)*smoothscale});
                //System.out.println("after: "+atom.toString());
                //System.out.println("deviation after: "+master.deviation);
            }
            
        }
    }*/
    void save_tifdata2(String path)
    {
        boolean success;
        success = (new File(path+"models")).mkdir();
        success = (new File(path+"sims")).mkdir();
        success = (new File(path+"diffs")).mkdir();
                            
        IJ.saveAsTiff(stackp, path+"models/model_of_views"+globalcounter+".tif");
        IJ.saveAsTiff(stacksimp, path+"sims/simulation_of_views"+globalcounter+".tif");
        IJ.saveAsTiff(diffp, path+"diffs/difference"+globalcounter+".tif");
        //IJ.saveAsTiff(merged, path+"sum_of_models.tif");

        String log= IJ.getLog();
        File file = new File(path+"log.txt");
        try 
        {
        FileWriter fw = new FileWriter(file, false);
        PrintWriter pw = new PrintWriter(fw);
        pw.println(log);
        pw.close();
        } catch (IOException ex) 
        {
            System.out.println("could not save "+path+"log.txt");
            Logger.getLogger(TwoDReconstruction.class.getName()).log(Level.SEVERE, null, ex);
        }       
        globalcounter+=1;
        
    }
    class Atom
    {
        float x;
        float y;
        int id;
        int view;
        float[] force=new float[]{Float.NaN,Float.NaN,Float.NaN}; 
        int element=6;
        float intensity_factor=1;
        //public ArrayList<Ring> rings =  new ArrayList<Ring>(3);
        
        Atom(float x, float y, int id, int view)
        {
            this.x=x;
            this.y=y;
            this.id=id;
            this.view=view;
        }
        
        void shift_coordinates(float[] shift)
        {
            this.set_coordinates(this.x+shift[0], this.y+shift[1]);
        }
        
        void set_coordinates(float x, float y)
        {
            this.x=x;
            this.y=y;
            masters.get(view).update_pixels();
            update_merit();
        }
    }
    
    class Atoms
    {
        ArrayList<Atom> atoms=new ArrayList<Atom>(2200);
        //ArrayList<Ring> rings=new ArrayList<Ring>(1000);
        int view;
        float[] imarray;
        float[] imarray_model;
        float[] simarray;
        float[] diffarray;
        float[] expIm;
        float difference=0;
        float backgroundlevel=0;
        float fov=3f;
        public Fakesimulation sim;
        
        Atoms(int view)
        {
            this.view=view;
            imarray=new float[impWidth*impHeight];
            imarray_model=new float[impWidth*impHeight];
            simarray=new float[impWidth*impHeight];
            diffarray=new float[impWidth*impHeight];
            sim=new Fakesimulation(imarray,fov, impWidth);
            expIm=ImageCalculations.normalize_image((float[])stack.getPixels(view+1),(float[])stack.getPixels(view+1));
        }
        
        public void update_pixels()
        {
            Arrays.fill(imarray, backgroundlevel);
            Arrays.fill(imarray_model, 0);
            Arrays.fill(simarray, 0.0f);
            float px;
            float py;
            if(!atoms.isEmpty())
            {
                for(Atom atom : atoms )  
                {
                    float intensity=(float) Math.pow(atom.element,1.6);
                    intensity*=atom.intensity_factor;//should compensate intensity fluctuations
                    px= atom.x;
                    py= atom.y;
                    //System.out.println(px+"\t"+py);
					//if(Float.isNaN(expIm[(int)(px+py*impWidth+0.5f)]))//This code makes troubles with edges
					//	continue;

                    if((int)px>0.5 && (int)(px+1)<impWidth  && ((int)py>0.5 && (int)(py+1)<impWidth ))
                    {
                        imarray[(int)px+(int)py*impWidth ]+=(1-(px-(int)px))*(1-(py-(int)py))*intensity;
                        imarray[(int)px+1+(int)py*impWidth ]+=(px-(int)px)*(1-(py-(int)py))*intensity;
                        imarray[(int)px+1+((int)py+1)*impWidth ]+=(px-(int)px)*(py-(int)py)*intensity;
                        imarray[(int)px+((int)py+1)*impWidth ]+=(1-(px-(int)px))*(py-(int)py)*intensity;
                    }
                    imarray_model[(int)(px+0.5)+(int)(py+0.5)*impWidth]+=intensity;
                    
                }
                
                simarray=ImageCalculations.normalize_image(this.simulate_image(qstem),expIm);
                calc_difference();
                stack_model.setPixels(imarray_model, view+1);
                stack_sim.setPixels(simarray,view+1); 
                diffSt.setPixels(diffarray, view+1);
                //System.out.println(Arrays.toString(view.diffarray));
        

            }   
            if((System.currentTimeMillis()-starttime)/1000.0>1)
            {
                check_imagestatus();
                //System.out.println("Update");
                stackp.updateAndRepaintWindow();  
                stacksimp.updateAndRepaintWindow();
                diffp.updateAndRepaintWindow();
                starttime=System.currentTimeMillis();
            }

            
        }
        public void calc_difference()
        {
            Arrays.fill(diffarray, 0.0f);
            this.difference=0;
            for (int i=0;i<impWidth *impWidth ;++i)
            {
                if(Float.isNaN(expIm[i]))
                    continue;
                diffarray[i]=(simarray[i]-expIm[i])*(simarray[i]-expIm[i]);
                
                this.difference+=diffarray[i];
            }
            this.difference=this.difference/(impWidth*impWidth);
            //difference=(float)ImageCalculations.img_match2(simarray,expIm);
        }
        
        
        public float[] simulate_image(boolean qstem)
        {
            /*if(qstem)
            {
                //update_relative_coordlist();
                //qsim.cfg.write_cfg(relative_coordlist,true);
                //return qsim.simulate_image();
            }
            else*/
            {
                return sim.simulate_image(imarray);
            }
        }
        
        /*public void update_rings()
        {
            for (Ring ring : this.rings)
            {
                ring.update();
            }
        }*/
    }
    /*
    class Ring
    {
        public ArrayList<Atom> atoms;
        public float[] center= new float[3];
        public int length;
        public float radius;
        int view;

        public Ring(ArrayList<Atom> atom, int view)
        {
            this.view=view;
            this.atoms=(ArrayList<Atom>) atoms;
            for(Atom at : atoms)
            {
                at.rings.add(this);
            }
            this.length= atoms.size();
            masters.get(view).rings.add(this); 
            this.update();
            //System.out.println(toString());
        }

        void calc_center()
        {
            float x=0;
            float y=0;
            float z=0;
            int counter=0;
            //System.out.println("\nnew ring");
            for (Atom atom: this.atoms)
            {
                x+=atom.x;
                y+=atom.y;
                counter+=1;
                //System.out.println(atom.toString());
            }
            center[0]=x/counter;
            center[1]=y/counter;

        }

        void calc_radius()
        {
            float dx=0;
            float dy=0;
            float dz=0;
            int counter=0;
            radius=0;
            for (Atom atom: this.atoms)
            {
                dx=(atom.x-center[0])*(atom.x-center[0]);
                dy=(atom.y-center[1])*(atom.y-center[1]);
                counter+=1;
                radius+=(float)Math.sqrt(dx+dy);
            }
            radius=radius/counter;

        }

        void update()
        {
            calc_center();
            calc_radius();
        }

        void get_neighbors()
        {

        }

        @Override
        public String toString()
        {
            return "Ring has radius "+radius+" and center "+Arrays.toString(center);
        }

    }*/

    public class Bond
    {
        public Atom a1,a2;
        public Bond( Atom a, Atom b )
        {
                a1 = a;
                a2 = b;
                bonds.add(this);
        }
    }
}



    
   


