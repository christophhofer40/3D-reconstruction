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
import christoph.Master;
import christoph.Fakesimulation;
import christoph.QSTEM_simulator;
import christoph.Quaternion;
import christoph.Shear;
import christoph.ImageCalculations;
import christoph.Optimization;
import java.awt.Color;
import java.awt.Frame;
import java.net.InetAddress;
import javax.swing.SwingConstants;


public class Reconstructor implements PlugIn 
{
    
    static Master master;
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
    
    static float[] masterarray;
    static float[] simarray;
    ArrayList<Double> meritcurve = new ArrayList<Double>();

    float fov=8;
    float stepsize=1;
    float threshold=30;
    int cropradius=20;
    float smoothweight=2f;
    float smoothscale=0.5f;
    float factor=0;
    String präfix;
   
    int globalcounter=0;
    int automatic_rounds=0;
    int iterations=10;
    int totalrounds=1;
    float energyinc=0.05f;
    
    
   

    int best_view=0;

    
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
    static boolean resetz=false;

    static boolean move_atoms=false;
    static boolean tune_quaternion=false;
    static boolean fix_tune_quaternion=false;
    static boolean fixtilt_tune=false;
    static boolean tune_skewmatrix=false;
    static boolean bpull=false;
    static boolean bpull2=false;
    static boolean smooth=false;
    static boolean smooth2=false;
    static boolean opt_fov=false;
    static boolean opt_fov2=false;
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
    static boolean bcompress=false;
    static boolean opt_trans=false;
    //static boolean bdeltaintensity=false;
    
    
    @Override
    public void run(String arg) 
    {
        Locale.setDefault(Locale.UK);
        
        while (DoDialog())
        {
        if(validateInput()!=0)
            {continue;}
            IJ.resetEscape(); 
            markToolbar("energymixing: "+factor);
            run_slice();
        }
    }
    
    public void markToolbar(String comment)
    {
        String hostname = null;
        try
        {	
            hostname = InetAddress.getLocalHost().getHostName();


            Frame[] activeframes = Frame.getFrames();
            activeframes[0].setTitle(System.getProperty("user.name")+"@"+ hostname+
                    (comment!=null?" # "+comment:""));
        }
        catch(Exception e)
        {
                hostname = "<none>";
                e.printStackTrace();
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
        
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("3D-Reconstructor");
        gd.setSmartRecording(true);


        gd.addChoice("raw image (GRAY32)", imptitleArray, impT);
        //gd.addChoice("2 bit model", imptitleArray, modelT);
        
        //gd.addNumericField("field of view [nm]",fov,2);
        gd.addNumericField("energy - deviation weighting", factor, 9);    
        
        gd.addMessage("Initialization");
        gd.addCheckbox("reset model", false);
        gd.addCheckbox("reset z", false);
        
        gd.addMessage("Move atoms");
        gd.addCheckbox("z-direction", move_atoms);
        gd.addCheckbox("random directions", rdmmove);
        gd.addCheckbox("orthogonal movement", orthmove);
        
        //gd.addCheckbox("optimize edgeatoms", edgeatoms);
        //gd.addNumericField("threshold for selecting atoms [px]", threshold, 0);
        //gd.addNumericField("radius for calculating local deviation [px]", cropradius, 0);
        gd.addNumericField("\t number of finetune iterations", iterations, 0);

        //gd.addCheckbox("pull back bad atom to center", bpull);
        //gd.addCheckbox("pull back bad atom to seen position", bpull2);
        
        gd.addMessage("Relax model geometrically");
        gd.addCheckbox("smooth z-coorinates", smooth);
        gd.addNumericField("\tweight of atom in the smooth (>=1 recommended)", smoothweight, 2);
        gd.addCheckbox("smooth rings", smooth2);
        gd.addNumericField("\tweight of the smooth (between 0 and 1 recommended)", smoothscale, 2);
        
        gd.addMessage("Adjust view parameters");
        gd.addCheckbox("optimize field of view (with energyminimization)", opt_fov);
        gd.addCheckbox("optimize field of view (with image correlation)", opt_fov2);
        gd.addCheckbox("tune azimuthal angle", fix_tune_quaternion);
        gd.addCheckbox("tune tilt (varies the tilt angle of the microscope)", fixtilt_tune);
        gd.addCheckbox("compensate distortions", tune_skewmatrix);
        
        gd.addMessage("Adjust beam parameters");
        gd.addCheckbox("adjust gaussianblur", bluropt);
        gd.addCheckbox("find aberrations", abopt);
        gd.addCheckbox("compress/stretch model", bcompress);
        gd.addCheckbox("optimize translation (recommended when 2D positions are not nicely defined)", opt_trans);
        //gd.addCheckbox("adjust delta intensity", bdeltaintensity); //doesnt make difference with same atoms
        //gd.addCheckbox("write xyz", xyz);
        gd.addCheckbox("switch to correlations as merit", switch_);
        gd.addNumericField("\t mean stepsize (sigma of normaldistribution in pixel)", stepsize, 2);
        //gd.addCheckbox("\t use qstem?", qstem);
        //gd.addCheckbox("if not: show simulations? (slower, but cooler)", show_sims);
        //gd.addCheckbox("show all moves? (slower, but cooler)", show_moves);
       
        gd.addNumericField(" totalrounds (-1...forever)", totalrounds, 0);
        //gd.addCheckbox("increase energy contribution automaticly?", automatic);
        //gd.addNumericField("\tby how much?", energyinc, 3);
        //gd.addNumericField("\t how many rounds? ", automatic_rounds, 0);
        
        gd.addCheckbox("analyze energy and image contribution", analyze);
        
        gd.showDialog();

        impT=gd.getNextChoice();
        //modelT=gd.getNextChoice();
        
        //fov = (float)gd.getNextNumber();
        factor=(float)gd.getNextNumber();
        
        reset=gd.getNextBoolean();
        resetz=gd.getNextBoolean();
        
        move_atoms=gd.getNextBoolean();
        rdmmove=gd.getNextBoolean();
        orthmove=gd.getNextBoolean();
        //edgeatoms=gd.getNextBoolean();
        
        //threshold=(float)gd.getNextNumber();
        //cropradius=(int)gd.getNextNumber();
        iterations=(int)gd.getNextNumber();

        //bpull=gd.getNextBoolean();
        //bpull2=gd.getNextBoolean();
        smooth=gd.getNextBoolean();
        smoothweight=(float)gd.getNextNumber();
        smooth2=gd.getNextBoolean();
        smoothscale=(float)gd.getNextNumber();
        
        opt_fov=gd.getNextBoolean();
        opt_fov2=gd.getNextBoolean();
        fix_tune_quaternion=gd.getNextBoolean();
        fixtilt_tune=gd.getNextBoolean();
        tune_skewmatrix=gd.getNextBoolean();    
        
        bluropt=gd.getNextBoolean();
        abopt=gd.getNextBoolean();
        bcompress=gd.getNextBoolean();
        opt_trans=gd.getNextBoolean();
        //bdeltaintensity=gd.getNextBoolean();
        
        //xyz=gd.getNextBoolean();
        switch_=gd.getNextBoolean();
        stepsize=(float) gd.getNextNumber();
        //qstem=gd.getNextBoolean();
        //show_sims=gd.getNextBoolean();
        //show_moves=gd.getNextBoolean();      
        
        totalrounds=(int)gd.getNextNumber();
        //automatic=gd.getNextBoolean(); 
        //energyinc=(float)gd.getNextNumber();
        //automatic_rounds=(int)gd.getNextNumber();

        analyze=gd.getNextBoolean(); 
        
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
        angFileStr=path+"angles.txt";
        topFile= new File(topFileStr);
        angFile=new File(angFileStr);
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
        
        boolean change=false;
        int nochange=0; //if overshoots threshold: next round or end
        
        if(reset || resetz)
        {
            globalcounter=0;
            try 
            {
                read_top(topFile);
                IJ.log("atom structure initialized");
                //master.update(true, true, false);

            } catch (IOException ex) 
            {
                Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
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
            
            master.init_lammps(path);
            try {
                read_beam_param(new File(path+"beamparameters.txt"));
                IJ.log("aberrations setted from "+path+"beamparameters.txt");
                
            } catch (FileNotFoundException ex) {
                Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
            }
            datafile=new File(path+"data.txt");
            
            master.factor=factor;
            
            show_views();
            master.update(true, true, true);
            append_datafile(datafile,true,true,"#resetted structure\t#energy scalefactor = "+factor);
        }
        
        if (master.qstem!=qstem)
        {
            master.qstem=qstem;
            change=true;
        }
        int it=0;
        if(master.fov!=fov)
        {
            master.change_fov(fov);
            write_beam_param();
            change=true;
            push_pixels();
        }
        if(master.factor!=factor)
        {
            master.factor=factor;
            append_datafile(datafile,true,true,"#new scalefactor = "+factor);
            change=true;
        }
        if(change)
        {
            master.update(false,true,true);
            push_pixels();
        }
        if(automatic)
            append_datafile(datafile,true,true,"starting automatic increase of energy contribution for "+totalrounds +" rounds each with starting factor "+factor);

        
        do
        {
            if (master.factor>1 || master.factor<0)
            {
                IJ.log("factor bigger than 1 or smaller than 0");
                break;
            }
            

            best_view=find_bestView()-1;
            //System.out.println(best_view);
            double oldmerit=master.merit;
            IJ.log("\n\nstarting new round... merit before = \t\b"+oldmerit);
            IJ.log("energy - deviation scalefactor: "+factor);
            IJ.log("merit before = "+master.merit+"\n");
            IJ.log("energy = " +master.energy);
            IJ.log("image deviation² = "+master.correlation);

            if (analyze)
                analyze_energy_deviation_contribution();
            
            if(move_atoms)
            {
                IJ.log("fine tuning atoms...");
                for(int k=0; k<iterations;k++)
                {
                    for(Master.View view:master.views)
                    {
                        if (!view.atoms2D.isEmpty())
                            doProcess();
                        break;
                        
                    }
                    if(IJ.escapePressed())
                         break;
                    master.update(true, true, true);
                    meritcurve.add(master.energy);
                    meritcurve.add(master.correlation);
                    meritcurve.add((double)master.mdeviation);
            
                }
            }
                    
            if(orthmove)
            {
                for(int k=0; k<iterations;k++)
                {
                    projection_move();
                    
                    if(IJ.escapePressed())
                        break;
                    master.update(true, true, true);
                    meritcurve.add(master.energy);
                    meritcurve.add(master.correlation);
                    meritcurve.add((double)master.mdeviation);
            
                }
            }
            if(rdmmove)
            {
                for(int k=0; k<iterations;k++)
                {

                    random_move();
                    
                    if(IJ.escapePressed())
                        break;
                    
                    master.update(true, true, true);
                    meritcurve.add(master.energy);
                    meritcurve.add(master.correlation);
                    meritcurve.add((double)master.mdeviation);
            
                }
            }
            if(edgeatoms)
            {
                for(int k=0; k<iterations;k++)
                {

                    optimize_edgeatoms();
                    
                    if(IJ.escapePressed())
                        break;
                    
                    master.update(true, true, true);
                    meritcurve.add(master.energy);
                    meritcurve.add(master.correlation);
            
                }
            }
            
            
            master.write_xyz(); 
            //master.update(true, true, true);
            push_pixels();
            write_top("");
            if(IJ.escapePressed())
                break;

            if(move_atoms||orthmove||rdmmove) 
            {
                if(master.merit!=oldmerit)
                {
                    append_datafile(datafile,true,false,"");
                    IJ.log("distance-deviation after finetune:\t\b" +master.mdeviation+" pixel");           
                    IJ.log("new merit  = "+master.merit+"\n");
                    IJ.log("image deviation² = "+master.mdeviation);
                }
                else
                    IJ.log("no change after fineoptimize");
                if(master.merit>oldmerit&&switch_)
                {
                    IJ.log("merit bigger again");
                    IJ.log("maybe rounding error");
                    nochange+=1;
                }
                else
                    nochange=0;
            }
            
            if(bpull)
            {
                IJ.log("pulling back bad atom");
                pull_back_to_circumcenter();
                IJ.log("new merit (might be worse: "+master.merit);
                
            }
            
            if (bpull2)
            {
                IJ.log("pulling back bad atom to seen pos");
                pull_back_to_seenPosition();
                master.update(true, true, true);
                IJ.log("new merit (might be worse: "+master.merit);
            }

            if (smooth)
            {
                smooth_structure();
                IJ.log("distance-deviation after smoothing z coordinates:\t\t" +master.mdeviation+" pixel");           
                IJ.log(" new merit  = "+master.merit+"\n");
                append_datafile(datafile,true,true,"#z-smooth");
            }

            if (smooth2)
            {
                smooth_structure2();
                IJ.log("distance-deviation after smoothing rings:\t\t" +master.mdeviation+" pixel");           
                IJ.log(" new merit  = "+master.merit+"\n");
                append_datafile(datafile,true,true,"#ringsmooth");
            }
            
            
            if(opt_fov)
            {
                oldmerit=master.merit;
                IJ.log("optimizing field of View");
                IJ.log("energy before: "+master.energy);
                optimize_fieldOfView();  
                if(oldmerit!=master.merit)
                {
                    IJ.log("new field of view: "+fov);
                    IJ.log("energy after :" +master.energy);
                    IJ.log("new merit  = "+master.merit+"\n");
                    append_datafile(datafile,true,true,"#field of view optimized");
                    
                }
                else
                    IJ.log("no change after optimizing field of view");
            }
            
            if(opt_fov2)
            {
                oldmerit=master.merit;
                IJ.log("optimizing field of View");
                IJ.log("correlation before: "+master.correlation);
                optimize_fieldOfView2();  
                if(oldmerit!=master.merit)
                {
                    IJ.log("new field of view: "+fov);
                    IJ.log("correlation after: "+master.correlation);
                    IJ.log("new merit  = "+master.merit+"\n");
                    append_datafile(datafile,true,true,"#field of view optimized");
                    
                }
                else
                    IJ.log("no change after optimizing field of view");
            }
            
            if(opt_trans)
            {
                IJ.log("optimizing translation");
                IJ.log("correlation before: "+master.correlation);
                optimize_translation();
                IJ.log("correlation after: "+master.correlation);
            }
            
            if(tune_quaternion)
            {
                //optimize_quaternions();
            }
            if(fix_tune_quaternion)
            {
                fix_tilt_tune3();
                
                IJ.log("\n\ndistance-deviation after "+master.mdeviation);
                for(Master.View view:master.views)
                    IJ.log("new inplaneangle of view "+view.slicenumber+": " +view.inplaneangle+"\n");
                IJ.log("new merit  = "+master.merit+"\n");
            }
            if (fixtilt_tune)
            {
                tune_fixtilt();
                IJ.log("\n\ndistance-deviation after "+master.mdeviation);
                for(Master.View view:master.views)
                    IJ.log("fixed angle of view "+view.slicenumber+": " +view.first_angle+"\n");
                IJ.log("new merit  = "+master.merit+"\n");
            }
            if(tune_skewmatrix)
            {
                optimize_skewmatrix2();
                IJ.log("skewmatrix opimized... new deviation: "+master.mdeviation);
                IJ.log("new matrix: "+master.views.get(0).skewmatrix.getStr()+"\n");
                IJ.log("new merit  = "+master.merit+"\n");
            }
            if(bluropt)
            {
                optimize_gaussianblur();
                IJ.log("new blur (in sigma): "+master.views.get(0).sim.sigma);
                IJ.log(" new merit  = "+master.merit+"\n");
                append_datafile(datafile,true,true,"#gaussianblur optimized");
                
            }
            if(abopt)
            {
                optimize_aberrations2();
                IJ.log(" new merit  = "+master.merit+"\n");
                IJ.log("optimized aberrations:");
                for(Master.View view:master.views)
                    IJ.log("view "+view.slicenumber+" [C1,C12a,C12b,C21a,C21b,C23a,C23b] = " + Arrays.toString((float[]) view.sim.beam.aberrations)+"\n");
                append_datafile(datafile,true,true,"#aberrations optimized");

            }
            if(bcompress)
            {
                IJ.log("stretching/compressing model...");
                modify_compression();
                IJ.log(" new merit  = "+master.merit+"\n");
            }
            
            /*if(bdeltaintensity)
            {
                IJ.log("adjusting delta intensity...");
                modify_deltaintensity();
                IJ.log(" new merit  = "+master.merit+"\n");
            }*/
            //if (xyz)
            //    write_xyz();
            
            //master.update(true, true, true);
            master.write_xyz();
            //save_tifdata2(path);
            push_pixels();
            write_top("");
            it++;
            if (nochange>=3)
            {
                it=totalrounds;
                IJ.log("\nstructure not improving");
            }
                
            if(it==totalrounds && automatic)
            {
                save_tifdata2(path,String.valueOf(factor));
                master.write_xyz(String.valueOf(factor));
                write_top(String.valueOf(factor));
                globalcounter+=1;
                IJ.log("top and xyz with energyfactur "+factor+ " has counter "+globalcounter);
                
                factor+=energyinc;
                master.factor=factor;
                master.update(true, true, true);
                append_datafile(datafile,true,true,"#automatic increase of energyfactor, new energy scalefactor = "+factor);
                it=0;
                
                IJ.log("increasing energy contribution to "+factor);
                markToolbar("energymixing: "+factor);
                if(factor>1 || globalcounter==automatic_rounds)
                {                  
                    IJ.log("automatic energy increase over");
                    break;
                }
            }
            
            
            if(IJ.escapePressed())
                break;
            
            
        }        
        while (it!=totalrounds);

        append_datafile(datafile,true,false,"");
    }
    
    void append_datafile(File file, boolean append, boolean adddata,String comment)
    {
        if(adddata)
        {
            meritcurve.add(master.energy);
            meritcurve.add(master.correlation);
            meritcurve.add((double)master.mdeviation);
            System.out.println((double)master.mdeviation);
        }
        
        try
        {
            FileWriter fw = new FileWriter(file, append);
            PrintWriter pw = new PrintWriter(fw);
            for( int i=0; i<meritcurve.size();i+=3)
            {
                pw.println(meritcurve.get(i)+"\t"+meritcurve.get(i+1)+"\t"+meritcurve.get(i+2)+"\t"+comment);
            }
            pw.close();
            IJ.log("rewrote " + file.getName());
            meritcurve.clear();
            
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
    
    void read_top(File topFile)
        throws FileNotFoundException, IOException 
    {   
        if (topFile == null || !topFile.canRead()) 
        {
            throw new IllegalArgumentException("file not readable: " + topFile);
        }
        FileReader topReader=new FileReader(topFile);
        final Scanner s= new Scanner(topFile).useDelimiter("\t");
        int chapter=-1;
        int counter=0;
        int viewcounts = -1001;
        int atomid = -1001;
        float x = Float.NaN;
        float y = Float.NaN;
        float z = Float.NaN;
        int a1;
        int a2;
        boolean add;
        //float[] quaternionVector=new float[4];
        float fixedangle;
        float inplaneangle;
        int lammpsid=1;
        while(s.hasNextLine())
        {    
            String line = s.nextLine();
            if(line.startsWith("#") || line.equals(""))
            {
                continue;
            }
            if(line.startsWith("VIEW"))
            {
                chapter+=1;
                counter=0;
                master.views.add(master.new View(chapter+1, (float[]) stack.getPixels(chapter+1)));
                if (chapter==0)
                    master.atomnumber=master.atoms.size();
                continue;
            }
            
            Scanner ls = new Scanner(line);
            if(line.startsWith("MASTER"))
            {
                ls.next();
                int atomnumber=ls.nextInt();
                master = new Master(atomnumber, impWidth, fov);
                master.set_xyzPath(path);               
                continue;

            }
            if(line.startsWith("LABEL"))
            {
                ls.next();
                master.views.get(chapter).label=ls.next();
                continue;
            }
            
            if(line.startsWith("STRETCH"))
            {
                ls.next();
                master.views.get(chapter).stretch_factor=ls.nextFloat();
                continue;
            }
            
            
            if(line.startsWith("TRANSLATION"))
            {
                ls.next();
                master.views.get(chapter).translation[0]=ls.nextFloat();
                master.views.get(chapter).translation[1]=ls.nextFloat();
                //master.views.get(chapter).translation[2]=ls.nextFloat();
                continue;
            }
            if(line.startsWith("SKEWMATRIX"))
            {
                ls.next();
                master.views.get(chapter).skewmatrix.map[0][0]=ls.nextFloat();
                master.views.get(chapter).skewmatrix.map[0][1]=ls.nextFloat();
                master.views.get(chapter).skewmatrix.map[1][0]=ls.nextFloat();
                master.views.get(chapter).skewmatrix.map[1][1]=ls.nextFloat();
                continue;
            }
            if(line.startsWith("QUATERNION"))
            {
                ls.next();
                /*quaternionVector[0]=ls.nextFloat();
                quaternionVector[1]=ls.nextFloat();
                quaternionVector[2]=ls.nextFloat();
                quaternionVector[3]=ls.nextFloat();
                master.views.get(chapter).quaternion=
                        new Quaternion(quaternionVector[0],quaternionVector[1]
                                ,quaternionVector[2],quaternionVector[3]);
                */
                fixedangle=ls.nextFloat();
                inplaneangle=ls.nextFloat();
                x = (float) Math.cos(inplaneangle);
                y = (float) Math.sin(inplaneangle);           
                float[] vec= new float[]{x,y,0.0f};
                master.views.get(chapter).quaternion= new Quaternion (fixedangle, vec);
                master.views.get(chapter).inplaneangle=inplaneangle*180/3.1415f;
                master.views.get(chapter).first_angle=fixedangle;   
                master.views.get(chapter).get_directionVector();
                
                continue;
            }
            if(line.startsWith("BACKGROUND"))
            {
                ls.next();
                master.views.get(chapter).backgroundlevel=ls.nextFloat();
                continue;
            }

            if(line.startsWith("ATOM"))
            {
                if(chapter==-1)
                {
                    int element=6;
                    try
                    {
                        ls.next();
                        atomid=ls.nextInt();  
                        viewcounts=ls.nextInt();
                        x=ls.nextFloat();
                        y=ls.nextFloat();
                        z=ls.nextFloat();
                        if (resetz)
                            z=0f;
                        
                        if(Float.isNaN(((float[])stack.getPixels(1))[(int)(x+0.5)+(int)(y+0.5)*impWidth]))
                        {
                            System.out.println("Atom at the edge not included");
                            viewcounts-=1;
                            
                        }
                        
                        if (ls.hasNextInt())
                        {                      
                            element=ls.nextInt();
                            System.out.println("Foreign atomnumer "+element);
                        }        
                                
                        System.out.println("reading atom: id: " + atomid + " viewcounts: " + viewcounts + " x: " + x + " y: " + y + " z: " +z);
                    }
                    catch(Exception e)
                    {
                        e.printStackTrace();
                        System.out.println("Something went wrong with reading topfile at line\n "+line);
                        System.out.println("id: " + atomid + " viewcounts: " + viewcounts + " x: " + x + " y: " + y + " z: " +z);
                        return;
                    }
                    master.atoms.add(master.new Atom(x,y,z,atomid,viewcounts));
                    if(element!=6)
                        master.atoms.get(atomid).element=element;
                    if(viewcounts>=stack.size() || stack.size()==1)
                    {
                        if (viewcounts>0)
                        {
                            master.seen_atoms+=1;
                            master.atoms.get(atomid).lammpsid=lammpsid;
                            master.atoms.get(atomid).viewcounts=stack.size();
                            lammpsid+=1;
                        }
                    }  
                    
                    continue;
                  
                }
                

            
                ls.next();
                atomid=ls.nextInt();
                x=ls.nextFloat();
                y=ls.nextFloat();
                
                //if(master.atoms.get(atomid).viewcounts<stackNb)
                  //  continue;
                if(atomid!=-1 && master.atoms.get(atomid).viewcounts==stack.size() )
                {  
                    Master.Atom2D atom2d=master.new Atom2D(x,y,master.atoms.get(atomid),
                                    master.views.get(chapter));
                    master.atoms.get(atomid).exp_atoms2D.add(atom2d);
                    master.views.get(chapter).atoms2D.add(atom2d);
                    if(atomid==0) 
                    {
                        master.views.get(chapter).atom0index=counter;
                    }
                    else
                        counter+=1;
                    while (ls.hasNext())
                    {
                        String kwarg=ls.next();
                        if(kwarg.startsWith("#intensity"))
                        {
                            //int elem=master.atoms.get(atomid).element;
                            float intensity=ls.nextFloat();
                            //System.out.println(intensity);
                            master.views.get(chapter).atoms2D.get( master.views.get(chapter).atoms2D.size()-1).intensity_factor=intensity;
                        }
                    }
                }
                
                //counter+=1;
                if(atomid==-1 ||  master.atoms.get(atomid).viewcounts<stack.size())
                {
                    
                    if(stack.size()==1)
                    {
                        z=0;
                        /*if(Float.isNaN(((float[])stack.getPixels(1))[(int)(x+0.5)+(int)(y+0.5)*impWidth]))
                        {
                            System.out.println("Atom at the edge not included");
                            continue;
                        }*/
                        try
                        {
                            atomid=master.atoms.size(); 
                            viewcounts=stack.size();
                            System.out.println("reading fakeatom: id: " + atomid + " viewcounts: " + viewcounts + " x: " + x + " y: " + y + " z: " +z);
                        }
                        catch(Exception e)
                        {
                            e.printStackTrace();
                            System.out.println("Something went wrong with reading topfile at line\n "+line);
                            System.out.println("id: " + atomid + " viewcounts: " + viewcounts + " x: " + x + " y: " + y + " z: " +z);
                            return;
                        }
                        master.atoms.add(master.new Atom(x,y,z,atomid,viewcounts));
                        if(viewcounts>1)
                        {
                            master.seen_atoms+=1;
                            master.atoms.get(atomid).lammpsid=lammpsid;
                            lammpsid+=1;
                        }  
                        Master.Atom2D atom2d=master.new Atom2D(x,y,master.atoms.get(atomid),
                                    master.views.get(chapter));
                        master.atoms.get(atomid).exp_atoms2D.add(atom2d);
                        master.views.get(chapter).atoms2D.add(atom2d);
                        
                        while (ls.hasNext())
                        {
                            String kwarg=ls.next();
                            if(kwarg.startsWith("#intensity"))
                            {
                                float intensity=ls.nextFloat();
                                //System.out.println(intensity);
                                master.views.get(chapter).atoms2D.get( master.views.get(chapter).atoms2D.size()-1).intensity_factor=intensity;
                            }
                        }

                        continue;
                    }
                    else
                    {
                        master.views.get(chapter).oneatoms.add(master.new Minus1Atom((int)x,(int)y));
                        while (ls.hasNext())
                        {
                            String kwarg=ls.next();                        
                            if(kwarg.startsWith("#intensity"))
                            {
                                float intensity=ls.nextFloat();
                                //System.out.println(intensity);
                                master.views.get(chapter).oneatoms.get( master.views.get(chapter).oneatoms.size()-1).intensity_factor=intensity;
                            }
                        }
                        
                        continue;
                    }
                }
                
            }
            if(line.startsWith("BOND"))
            {
                ls.next();
                a1= ls.nextInt();
                a2= ls.nextInt();
                if(master.atoms.get(a1).viewcounts>1 &&master.atoms.get(a2).viewcounts>1)
                {
                    master.atoms.get(a1).neighbors.add(master.atoms.get(a2));
                    master.atoms.get(a2).neighbors.add(master.atoms.get(a1));
                    master.new Bond(master.atoms.get(a1),master.atoms.get(a2));
                }
            }
            if(line.startsWith("RING"))
            {
                ls.next();
                ArrayList<Master.Atom> ratoms =new ArrayList<Master.Atom>();        
                add=true;
                while(ls.hasNextInt())
                {
                    a1= ls.nextInt();
                    ratoms.add( master.atoms.get(a1));
                    if(master.atoms.get(a1).viewcounts<2)
                        add=false;
                }
                if (add)
                    master.new Ring(ratoms); 
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
                fov=ls.nextFloat();
                master.change_fov(fov);
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
                master.views.get(slice).sim.sigma=sig;                          
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
                        master.views.get(slice).
                                sim.beam.aberrations[k]=ab;  			
                }

            }
			
            
        }
        for(int sl=0;sl<master.views.size();sl++)
        {
        	master.views.get(sl).sim.fov=master.views.get(sl).sim.fov;
        	master.views.get(sl).sim.beam.initBeam();
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
            for( int i=0; i<master.views.size();i++)
            {
            	pw.println("VIEW:\t"+i);
                pw.println("Aberrations:\t"+master.views.get(i).sim.beam.aberrations[0]+"\t"
                                +master.views.get(i).sim.beam.aberrations[1]+"\t"
                                +master.views.get(i).sim.beam.aberrations[2]+"\t"
                                +master.views.get(i).sim.beam.aberrations[3]+"\t"
                                +master.views.get(i).sim.beam.aberrations[4]+"\t"
                                +master.views.get(i).sim.beam.aberrations[5]+"\t"
                                +master.views.get(i).sim.beam.aberrations[6]+"\t");
                pw.println("Sigma:\t"+ +master.views.get(i).sim.sigma);
                pw.println("FieldOfView:\t"+master.views.get(i).sim.fov);
                //beamstack.setPixels((float[])master.views.get(i).sim.beam.beamfht.getPixels(),i+1);
            }
            pw.close();
            IJ.log("rewrote " + file.getName());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
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
    
    
    
    void translate_atoms()
    {
        for(Master.View view:master.views)
        {
            view.update_translation();
        }
    }*/
    
    

    int find_bestView()
    {
        float bestdev=-1;
        int best_view=-1;
        int deviation=0;
        for(Master.View view : master.views)
        {
            if(view.atoms2D.size()!=0)
            {
                deviation=0;
                for(Master.Atom2D atom2d: view.atoms2D)
                {
                    deviation+=atom2d.ldeviation;
                }
                if(deviation < bestdev || bestdev==-1)
                {
                    bestdev=deviation;
                    best_view=view.slicenumber-1;
                }
            }
        }
        //IJ.log("best view = "+best_view);        
        return best_view;
    }

    
    void show_views()
    {
        /*stackp=WindowManager.getImage("model_of_views.tif");
        if (stackp!=null)
            stackp.close();
        stacksimp=WindowManager.getImage("simulation_of_views.tif");
        if (stacksimp!=null)
            stacksimp.close();
        merged=WindowManager.getImage("sum_of_models.tif");
        if (merged!=null)
            merged.close();
        diffp=WindowManager.getImage("difference.tif");
        if (diffp!=null)
            diffp.close();
        */
        stackp=WindowManager.getImage("model_of_views.tif");
        stacksimp=WindowManager.getImage("simulation_of_views.tif");
        //merged=WindowManager.getImage("sum_of_models.tif");
        diffp=WindowManager.getImage("difference.tif");
        
        stack_model=new ImageStack(impWidth,impHeight,stackNb);
        stack_sim=new ImageStack(impWidth,impHeight,stackNb);
        //mergedSt=new ImageStack(impWidth,impHeight,stackNb);
        diffSt=new ImageStack(impWidth,impHeight,stackNb);
        //for (Master.Atom atom:master.atoms)
         //   atom.update2D();
        master.update(true,true,true);
        master.relative_coordlist=new float[master.seen_atoms][3];
        for(Master.View view: master.views)
        {
            view.update_pixels(true);           
            stack_model.setPixels(view.imarray, view.slicenumber);
            stack_sim.setPixels(view.simarray, view.slicenumber); 
            diffSt.setPixels(view.diffarray, view.slicenumber);
            view.update_correlation();
            //System.out.println(Arrays.toString(view.diffarray));
        }
        if(stackp==null)
            stackp=new ImagePlus("model_of_views",stack_model);
        else
            stackp.setStack(stack_model);
        if(stacksimp==null)
            stacksimp=new ImagePlus("simulation_of_views",stack_sim);
        else
            stacksimp.setStack(stack_sim);
        //add_models();
        //if(merged==null)
        //    merged=new ImagePlus("sum of models",mergedSt);
        //else
        //    merged.setStack(mergedSt);
        if(diffp==null)
            diffp=new ImagePlus("difference",diffSt);
        else
            diffp.setStack(diffSt);
        stackp.show();
        stacksimp.show();
        //merged.show();
        diffp.show();
        
    }
    
    void doProcess()
    {
        //master.update_deviations();
        do
        {
            //System.out.println(best_view);
            best_view+=1;
            best_view=best_view%stackNb;
        }
        while(master.views.get(best_view).atoms2D.isEmpty());
        //float[] axisVector=master.views.get(best_view).quaternion.getAxisVector();  
        float[] shift = new float[3];
        master.views.get(best_view).get_directionVector();
        shift=master.views.get(best_view).directionVector; 

        float length=(float) (random.nextGaussian()+stepsize);
    
        
        /*for(int a=1;a<master.atomnumber;a++)
        {
            if (master.atoms.get(a).viewcounts>1)
                optimize_atom(a,shift);
            
            if(IJ.escapePressed())
                    {	break;}
        }*/
        
        //to do: change optimize moreatoms to make it compatible with moving more atoms
        /*if (switch_&&factor==0.0f)
            optimize_moreatoms(get_random_atoms(),shift);
        
        else*/
        float inorm = 1/(float)Math.sqrt(shift[0]*shift[0]+
        shift[1]*shift[1]+shift[2]*shift[2]);
        
        //int sucess=0;
        shift[0]=(shift[0]*inorm)*length;
        shift[1]=(shift[1]*inorm)*length;
        shift[2]=(shift[2]*inorm)*length;
        for(int a=1;a<master.atomnumber;a++)
            {
                if (master.atoms.get(a).viewcounts>1)
                    optimize_atom2(a,shift.clone(),false);

                if(IJ.escapePressed())
                        {	break;}
            }
        
        for(Master.View view:master.views)
        {
            view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
            view.update_correlation();
        }
        push_pixels();
    }
    
    int[] get_random_atoms()
    {
        ArrayList<Integer> rdm = new ArrayList<Integer>();
        ArrayList<Integer> atomlist = new ArrayList<Integer>();
        atomlist.clear();
        int start=1;
        boolean add;
        for(int i=start;i<master.atoms.size();++i)
        {
            rdm.add(i);
        }
        Collections.shuffle(rdm);
        int t=0; 
        while (master.atoms.get(rdm.get(t)).exp_atoms2D.size()<2)
        {
            t+=1;
        }
        atomlist.add(rdm.get(t));
        for(int i=t+1;i<master.atoms.size()-start;++i)
        {
            Master.Atom a = master.atoms.get(rdm.get(i));
            add=true;
            if (a.exp_atoms2D.size()<2)
            {
                add=false;
                continue;
            }
            for(int k : atomlist)
            {
                Master.Atom other= master.atoms.get(k);
                
                for (int j=0;j<a.exp_atoms2D.size();j++)
                {
                    
                    //System.out.println(threshold);
                    //System.out.println(other.exp_atoms2D.get(j).distance2D(a.exp_atoms2D.get(j)));
                    //System.out.println();
                    if(other.exp_atoms2D.get(j).distance2D(a.exp_atoms2D.get(j))<threshold)
                    {
                        //System.out.println(other.exp_atoms2D.get(j).distance2D(a.exp_atoms2D.get(j)));
                        add=false;
                        break;
                    }
                }
                if(add==false)
                    break;
            }
            if (add==true)
            {
                atomlist.add(rdm.get(i));
                //System.out.println(rdm.get(i));
            }   
        }
        int[] ret=new int[atomlist.size()];
        //System.out.println(atomlist.size()+" atoms added");
        for(int i=0; i<atomlist.size();++i)
        {
            ret[i]=atomlist.get(i);
        }
        return ret;
    }
    
    
    
    //to do: make it compatible with energies
    void optimize_moreatoms(int[] atomnumbers,float[] vector)
    {
        diffp.deleteRoi();
        ImagePlus temp=diffp.duplicate();
        int x;
        int y;
        int counter=0;
        int count=atomnumbers.length;
        double newdev;
        double firstenergy;
        double[] tempenergies= new double[atomnumbers.length];
        float[] shift= new float[3];
        boolean sign=random.nextBoolean();
        float inorm = 1/(float)Math.sqrt(vector[0]*vector[0]+
        vector[1]*vector[1]+vector[2]*vector[2]);
        float length=(float) (random.nextGaussian()*stepsize);

        if (!sign)
        {
            length=-length;
        }
        shift[0]=(vector[0]*inorm)*length;
        shift[1]=(vector[1]*inorm)*length;
        shift[2]=(vector[2]*inorm)*length;
        
        for(int i : atomnumbers)
        {
            Master.Atom atom=master.atoms.get(i);
            atom.sumOfStats=0;
            firstenergy=master.energy;
            for(Master.View view:master.views)
            {
                Master.Atom2D atom2d= atom.exp_atoms2D.get(view.slicenumber-1);
                x=(int)(atom2d.expx+0.5);
                y=(int)(atom2d.expy+0.5);
                temp.setPosition(view.slicenumber);
                temp.setRoi(x-cropradius, y-cropradius, 2*cropradius+1, 2*cropradius+1);          
                atom.sumOfStats+=temp.getStatistics().mean;
                temp.deleteRoi();
                
            }
            atom.shift_coordinate(shift,false,false);
            tempenergies[counter]=firstenergy-master.energy;
            //atom.sumOfStats=atom.sumOfStats*(1-master.factor)+(firstenergy-master.energy)*master.factor;
            counter+=1;
        }
        master.update(true, true, true);
        diffp.deleteRoi();
        temp=diffp.duplicate();
        //System.out.println();
        counter=0;
        for(int i : atomnumbers)
        {
            Master.Atom atom=master.atoms.get(i);
            newdev=0;
            for(Master.View view:master.views)
            {
                Master.Atom2D atom2d= atom.exp_atoms2D.get(view.slicenumber-1);
                x=(int)(atom2d.expx+0.5);
                y=(int)(atom2d.expy+0.5);
                temp.setPosition(view.slicenumber);
                temp.setRoi(x-cropradius, y-cropradius, 2*cropradius+1, 2*cropradius+1);          
                newdev+=temp.getStatistics().mean;    
                temp.deleteRoi();
                
            }

            
            if(newdev>=atom.sumOfStats)
            {
                atom.shift_coordinate(new float[]{-shift[0],-shift[1],-shift[2]},false, false);
                count-=1;
            }
               
        }
        
        master.update(true, true, true);

        push_pixels();
        System.out.println(count +" atoms moved");
        System.out.println("new corr after fineoptimize = "+master.correlation);
    }
    
    /*
    void optimize_atom(int atomnumber,float[] vector)
    {
        
        float olddev=master.deviation;
        double oldcorr=master.correlation;
        float olddist=master.atoms.get(atomnumber).sum_distances();
        float inorm = 1/(float)Math.sqrt(vector[0]*vector[0]+
        vector[1]*vector[1]+vector[2]*vector[2]);
        float[] shift= new float[3];
        float length = 3;
        boolean sign=random.nextBoolean();
        if(switch_)
            length=(float) (random.nextGaussian()/4+stepsize);
        if (!sign)
        {
            length=-length;
        }
        do
        {
            shift[0]=(vector[0]*inorm)*length;
            shift[1]=(vector[1]*inorm)*length;
            shift[2]=-(vector[2]*inorm)*length;
            master.shift_coordinate(atomnumber, shift);

            //show_and_wait();
            if(!switch_)
            {
                if(master.deviation<olddev)
                {
                    length*=1.1;
                    olddev=master.deviation;
                    
                }
                else
                {
                    shift[0]=-shift[0];
                    shift[1]=-shift[1];
                    shift[2]=-shift[2];
                    master.shift_coordinate(atomnumber, shift);
                    length*=-0.7;
                }
            }
            else
            {
                for(Master.View view:master.views)
                {
                    view.update_model();
                    view.update_correlation();
                }

                //System.out.println(oldcorr);
                //System.out.println(master.correlation);
                //System.out.println();
                if(show_moves)
                    push_pixels();
                            

                
                if(master.correlation<oldcorr && master.atoms.get(atomnumber).sum_distances()*1.5<=olddist)
                {
                    
                    oldcorr=master.correlation;
                    System.out.println(master.correlation);
                    System.out.println(olddist);
                    System.out.println(master.atoms.get(atomnumber).sum_distances());
                    System.out.println();
                }
                else
                {
                    shift[0]=-shift[0];
                    shift[1]=-shift[1];
                    shift[2]=-shift[2];
                    master.shift_coordinate(atomnumber, shift);
                    
                }
                break;
                
            }
        
        }while(Math.abs(length)>=0.1);
        
        //System.out.println("before: "+ olddev);
        //System.out.println("after: "+ master.deviation);
        push_pixels();

        

    }*/
    
    
    void optimize_atom2(int atomnumber,float[] shift, boolean energy) //breaks after 2 successed rounds
    {
        //master.update_deviations();
        float olddev=master.mdeviation;
        double oldcorr;
        
        int sucess=0;
        float length=1;
        do
        {
            oldcorr=master.merit;         
            shift[0]*=length;
            shift[1]*=length;
            shift[2]*=length;
            master.shift_coordinate(atomnumber, shift,(switch_&&factor<1),switch_);
            if(show_moves)
                    push_pixels();
            //show_and_wait();
            if(!switch_)
            {
                //System.out.println("before: "+olddev+ " \tafter: "+master.deviation);
                if(master.mdeviation<olddev)
                {
                    length*=1.1;
                    olddev=master.mdeviation;
                    sucess+=1;
                    if (sucess>=2)
                        break;
                    
                }
                else
                {
                    shift[0]=-shift[0];
                    shift[1]=-shift[1];
                    shift[2]=-shift[2];
                    master.shift_coordinate(atomnumber, shift,switch_,energy);
                    length*=-0.7;
                    
                }
            }
            else
            {
                //System.out.println(oldcorr);
                //System.out.println(master.correlation);
                //System.out.println();
                if(master.merit<oldcorr)
                {
                    
                    System.out.println("new merit = "+master.merit);
                    break;
                }
                else
                {
                    shift[0]=-shift[0];
                    shift[1]=-shift[1];
                    shift[2]=-shift[2];
                    master.shift_coordinate(atomnumber, shift,switch_&&factor<1,switch_);
                    length*=-1;
                    sucess+=1;
                    if (sucess==2)
                        break;
                    
                    
                }
                //break;
                
            }
        
        }while(Math.abs(length)>=0.1);
        
        //System.out.println("before: "+ olddev);
        //System.out.println("after: "+ master.deviation);
        push_pixels();

        

    }
    
    //finds best solution of three colinear points
    //if constraint==-1: no constraint
    void optimize_atom3(int atomnumber,float[] vector2, float constraint) 
    {
        Master.Atom atom=master.atoms.get(atomnumber);
        float[] vector=atom.force.clone();
        if(Float.isNaN(vector[0]))
        {
            vector[0]=vector2[0];
            vector[1]=vector2[1];
            vector[2]=vector2[2];
        }
        float oldx=atom.x;
        float oldy=atom.y;
        float oldz=atom.z;
        double oldmerit=master.merit; 
        float inorm = 1/(float)Math.sqrt(vector[0]*vector[0]+
        vector[1]*vector[1]+vector[2]*vector[2]);
        float[] shift= new float[3];
        int count=0;
        float length=(float) (random.nextGaussian()*stepsize);
        if(constraint!=-1)
            
            do
            {
                constraint=atom.sum_distances();
                //System.out.println("before: "+constraint);
                do
                {
                    length = (float) (random.nextGaussian()*stepsize);
                }while (length==0);
                //System.out.println("correlation before: "+master.correlation);
                //System.out.println("energy before:      "+master.energy); 
                if(count!=0)
                {
                    do
                    {
                        vector[0]=random.nextFloat()-0.5f;
                        vector[1]=random.nextFloat()-0.5f;
                        vector[2]=random.nextFloat()-0.5f;
                        inorm=1/(vector[0]*vector[0]+vector[1]*vector[1]+vector[2]*vector[2]);
                    }while(1/inorm>0.25f || inorm==Float.NaN);
                }
                shift[0]=(vector[0]*inorm)*length;
                shift[1]=(vector[1]*inorm)*length;
                shift[2]=(vector[2]*inorm)*length;   
                count+=1;
                
            }while(atom.sum_newdistances(shift)>constraint && count<20);
        shift[0]=(vector[0]*inorm)*length;
        shift[1]=(vector[1]*inorm)*length;
        shift[2]=(vector[2]*inorm)*length;  
        if(count>=20)
        {
            System.out.println("good enough optimized");
            return;
        }
        float[] p1=new float[3];
        float[] p2= new float[3];
        float p1grad;
        float p2grad;
        float norm;
        //System.out.println("pos before: "+ atom.x+ " "+ atom.y+ " "+ atom.z);
        master.shift_coordinate(atomnumber, shift,(true&&factor<1),true);
        //System.out.println("after:  "+atom.sum_distances());
        p1grad=(float)(master.merit-oldmerit);
        if (p1grad<0)
        {         
            System.out.println("new merit1: "+master.merit);
            push_pixels();
            atom.force[0]=shift[0];
            atom.force[1]=shift[1];
            atom.force[2]=shift[2];
            return;
        }
        p1[0]=atom.x;
        p1[1]=atom.y;
        p1[2]=atom.z;
        
        shift[0]=-2*shift[0];
        shift[1]=-2*shift[1];
        shift[2]=-2*shift[2];
        master.shift_coordinate(atomnumber, shift,true&&factor<1,true);        
        p2grad=(float)(master.merit-oldmerit);
        if (p2grad<0)
        {
            System.out.println("new merit2: "+master.merit);
            push_pixels();
            atom.force[0]=0.5f*shift[0];
            atom.force[1]=0.5f*shift[1];
            atom.force[2]=0.5f*shift[2];
            return;
        }
        p2[0]=atom.x;
        p2[1]=atom.y;
        p2[2]=atom.z;
        

        master.set_coordinates(atomnumber, oldx,oldy,oldz,false,false);
        Optimization opt = new Optimization(atom.x,atom.y,atom.z,p1,p1grad,p2,p2grad);
        shift=opt.get_shift(0);
        norm=shift[0]*shift[0]+shift[1]*shift[1]+shift[2]*shift[2];
        //System.out.println("p1: "+Arrays.toString(p1));
        //System.out.println("p2: "+Arrays.toString(p2));
        //System.out.println(Arrays.toString(shift));
        
        if(norm>2*Math.pow(impWidth/(fov*10000),2)&& norm<Math.pow(impWidth/(fov*10),2))
        {
            master.shift_coordinate(atomnumber, shift,true&&factor<1,true);
            //System.out.println("pos after: "+ atom.x+ " "+ atom.y+ " "+ atom.z);
            if(master.merit>oldmerit)
            {
                System.out.println("bad merit : "+master.merit);
                master.set_coordinates(atomnumber, oldx,oldy,oldz,true&&factor<1,true);
                //System.out.println("correlation after reset: "+master.correlation);
                //System.out.println("energy after reset:      "+master.energy);
                //master.merit=oldmerit;
                //master.correlation=oldcorr;
                //master.energy=oldenergy;
                
            }
            else
            {
                System.out.println("new merit : "+master.merit);
                
            }
                
        }
        else
        {
            if(norm>=Math.pow(impWidth/(fov*10),2))
                System.out.println("shift too long : " + Math.sqrt(norm));  
            else
                System.out.println("shift too short: " + Math.sqrt(norm));
            master.set_coordinates(atomnumber, oldx,oldy,oldz,true&&factor<1,true);
            //master.merit=oldmerit;
            //master.correlation=oldcorr;
            //master.energy=oldenergy;
        }

        atom.force[0]=vector2[0];
        atom.force[1]=vector2[1];
        atom.force[2]=vector2[2];
        push_pixels();
        //System.out.println("before: "+ olddev);
        //System.out.println("after: "+ master.deviation);
        

    }
    
    float find_max_dist()
    {
        float maxdist=-1;
        float dist=-1;
        for(Master.Atom atom : master.atoms)
        {
            dist=atom.sum_distances();
            if(dist>maxdist)
                maxdist=dist;
        }
        return maxdist;
    }
    
    void optimize_edgeatoms()
    {
        
    }
    
    void projection_move()
    {
        for(int a=1;a<master.atomnumber;a++)
        {
            if (master.atoms.get(a).viewcounts<stackNb)
                continue;
            float[] shift= new float[3];
            int badview=-1;
            float deviation=-1;
            for(Master.Atom2D atom2d :master.atoms.get(a).exp_atoms2D)
            {
                if(atom2d.ldeviation>deviation)
                {
                    //System.out.println(atom2d.deviation);
                    deviation=atom2d.ldeviation;
                    badview=atom2d.view.slicenumber-1;
                }
            }
            Quaternion conj = master.views.get(badview).quaternion.conjugate();
            Master.Atom2D bad=master.atoms.get(a).exp_atoms2D.get(badview);
            shift[0]=bad.expx-bad.realx;
            shift[1]=bad.expy-bad.realy;
            shift[2]=0;
            master.atoms.get(a).prefered_dir=shift;
            shift=conj.rotate(shift);
            //if(!switch_)
                optimize_atom2(a,shift,false);
            if(IJ.escapePressed())
                    {	break;}
          
        }
        //if(switch_)
            //optimize_moreatoms_proj(get_random_atoms());
    }
    
    /*
    void optimize_moreatoms_proj(int[] atomnumbers)
    {
    
        diffp.deleteRoi();
        ImagePlus temp=diffp.duplicate();
        boolean tswitch=switch_;
        switch_=false; // not simulating when shifting atom
        show_sims=true; //update after shifting all the atoms
        int x;
        int y;
        int count=atomnumbers.length;
        double newdev;
        
        
        for(int i : atomnumbers)
        {
            Master.Atom atom=master.atoms.get(i);
            float[] vector= new float[3];
            vector=atom.prefered_dir;
            boolean sign=random.nextBoolean();
            float inorm = 1/(float)Math.sqrt(vector[0]*vector[0]+
            vector[1]*vector[1]+vector[2]*vector[2]);
            float length=(float) (random.nextGaussian()/4+stepsize);

            if (!sign)
            {
                length=-length;
            }
            atom.prefered_dir[0]=(vector[0]*inorm)*length;
            atom.prefered_dir[1]=(vector[1]*inorm)*length;
            atom.prefered_dir[2]=(vector[2]*inorm)*length;
            //System.out.println(Arrays.toString(atom.prefered_dir));
 
            atom.sumOfStats=0;
            
            for(Master.View view:master.views)
            {
                Master.Atom2D atom2d= atom.exp_atoms2D.get(view.slicenumber-1);
                x=(int)(atom2d.expx+0.5);
                y=(int)(atom2d.expy+0.5);
                temp.setPosition(view.slicenumber);
                temp.setRoi(x-cropradius, y-cropradius, 2*cropradius+1, 2*cropradius+1);          
                atom.sumOfStats+=temp.getStatistics().mean;
                temp.deleteRoi();
                
            }
            atom.shift_coordinate(atom.prefered_dir);
        }
        for(Master.View view:master.views)
        {
            view.update_pixels();
        }
        diffp.deleteRoi();
        temp=diffp.duplicate();
        //System.out.println();
        int k=0;
        for(int i : atomnumbers)
        {
            k+=1;
            Master.Atom atom=master.atoms.get(i);
            newdev=0;
            atom.prefered_dir[0]=-atom.prefered_dir[0];
            atom.prefered_dir[1]=-atom.prefered_dir[1];
            atom.prefered_dir[2]=-atom.prefered_dir[2];
           // System.out.println(Arrays.toString(atom.prefered_dir));
            for(Master.View view:master.views)
            {
                Master.Atom2D atom2d= atom.exp_atoms2D.get(view.slicenumber-1);
                x=(int)(atom2d.expx+0.5);
                y=(int)(atom2d.expy+0.5);
                temp.setPosition(view.slicenumber);
                temp.setRoi(x-cropradius, y-cropradius, 2*cropradius+1, 2*cropradius+1);          
                newdev+=temp.getStatistics().mean;    
                temp.deleteRoi();
                
            }

            
            if(newdev>=atom.sumOfStats)
            {
                atom.shift_coordinate(atom.prefered_dir);
                count-=1;
            }
               
        }
        
        for(Master.Atom atom:master.atoms)
        {
            atom.update2D();
        }
        for(Master.View view:master.views)
        {
            view.update_model();
            view.update_correlation();
        }
        switch_=true;
        push_pixels();
        switch_=tswitch;
        System.out.println(count +" atoms moved");
        System.out.println("new corr after dineoptimize = "+master.correlation);
    
    }
    */
    
    void random_move()
    {
        float[] shift = new float[3];
        float constraint;
        if(factor!=-1)
            constraint=-1;
        else
        {
            constraint=find_max_dist();
            System.out.println("maxdist = "+constraint);
        }
        //if(!switch_|| factor!=0.0f)
        {
            for(int a=1;a<master.atomnumber;a++)
            {
                int r=0;
                if (master.atoms.get(a).viewcounts>1)
                {
                    do
                    {
                        shift=get_randomVector();
                        if(rdmmove)
                            if(switch_)
                                optimize_atom3(a,shift,constraint);
                            else
                                optimize_atom2(a,shift,false);

                        r++;

                    }
                    while(r<1);

                }
                if(IJ.escapePressed())
                        {	break;}
            }
        }
        //else
          //  optimize_moreatoms(get_random_atoms(),shift);
        
        
    }
    
    float[] get_randomVector() //length==1
    {

        float[] shift = new float[3];
        float norm=0;
        do
        {
            shift[0]=random.nextFloat()-0.5f;
            shift[1]=random.nextFloat()-0.5f;
            shift[2]=random.nextFloat()-0.5f;
            norm=shift[0]*shift[0]+shift[1]*shift[1]+shift[2]*shift[2];
        }while(norm>0.25f || norm==0.0f);
        shift[0]=shift[0]/norm;
        shift[1]=shift[1]/norm;
        shift[2]=shift[2]/norm;
        return shift;
        
    }
    void pull_back_to_circumcenter()
    {
        Master.Atom badatom=master.atoms.get(0);
        float longest=-1;
        for(Master.Atom atom:master.atoms)
        {
            float avbl=atom.get_avgBondlength();
            if(avbl>longest)
            {
                longest=avbl;
                badatom=atom;
            }
        }
        System.out.println("bad atom has coordinates in view 1: " +badatom.exp_atoms2D.get(0).realx+"\t"+badatom.exp_atoms2D.get(0).realy);
        float[] dir=badatom.get_centerDirection(0.4f);
        if (dir[0]!=Float.NaN)
        {
            badatom.shift_coordinate(dir,true,true);
        }
                
    }
    
    void pull_back_to_seenPosition()
    {
        Master.Atom badatom=null; 
        float largestdev=0;
        for(Master.Atom atom:master.atoms)
        {
            float locdev=0;
            for(Master.Atom2D at: atom.exp_atoms2D)
            {
                locdev+=at.ldeviation;        
                if(locdev>largestdev)
                {
                    largestdev=locdev;
                    badatom=atom;
                }
            }
        }
        System.out.println("bad atom has coordinates in view 1: " +badatom.exp_atoms2D.get(0).realx+"\t"+badatom.exp_atoms2D.get(0).realy);
        float[] shift=new float[3];
        shift=get_randomVector();
        if(badatom!=null)
            optimize_atom2(badatom.id,shift,true);
        
        
    
    }
    
    void write_contributions(int atomnumber, float[] shift)
    {
        FileWriter fw = null ;
        File energyfile =new File(path+"energy_image_contribution.txt");
        try 
        {
            float ang2px=impWidth/(fov*10);
            fw = new FileWriter(energyfile, true);
            PrintWriter pw = new PrintWriter(fw);
            if (atomnumber==0)
            {
                pw.println("#analyzing energy and image contribution of an optimized model");
                pw.println("#each atom is shifted in different and their deviation of the energy and image correlation is documented");
                pw.println("#id\tshift\tdEnergy [eV/Atom]\tdImage");               
            }
            double oldenergy=master.energy;
            double oldcorr= master.correlation;
            
            Master.Atom atom=master.atoms.get(atomnumber);
            if(atom.viewcounts<=1)
                return;

            float oldx =atom.x;
            float oldy =atom.y;
            float oldz =atom.z;
            //float length=(float) (random.nextGaussian()*stepsize);  
            float length=random.nextFloat()*stepsize*3;
            shift[0]*=length;
            shift[1]*=length;
            shift[2]*=length;
            atom.shift_coordinate(shift,true,true);
            pw.print(atomnumber+"\t");
            pw.print(length/ang2px+"\t");
            pw.print((master.energy-oldenergy)+"\t");
            pw.println((master.correlation-oldcorr));
            if(show_moves)
                push_pixels();
            atom.set_coordinates(oldx,oldy,oldz,true,true);

            

        } catch (IOException ex) {
            Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    
    void push_pixels()
    {
        for(Master.View view : master.views)
        {
            if (!view.atoms2D.isEmpty())
            {
                view.update_pixels(switch_);
                stack_model.setPixels(view.imarray, view.slicenumber);
                stack_sim.setPixels(view.simarray, view.slicenumber);
                diffSt.setPixels(view.diffarray, view.slicenumber);
            }
        }
        stackp.updateAndRepaintWindow();
        stacksimp.updateAndRepaintWindow();
        diffp.updateAndRepaintWindow();
        //add_models();
        //merged.updateAndRepaintWindow();
        
    }
    void save_tifdata(String path)
    {
        IJ.saveAsTiff(stackp, path+"model_of_views.tif");
        IJ.saveAsTiff(stacksimp, path+"simulation_of_views.tif");
        IJ.saveAsTiff(diffp, path+"difference.tif");
        //IJ.saveAsTiff(merged, path+"sum_of_models.tif");

        String log= IJ.getLog();
        File file = new File(path+"log.txt");
        master.write_xyz();
        try 
        {
        FileWriter fw = new FileWriter(file, false);
        PrintWriter pw = new PrintWriter(fw);
        pw.println(log);
        pw.close();
        } catch (IOException ex) 
        {
            System.out.println("could not save "+path+"log.txt");
            Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
        }   
    }
    
    void save_tifdata2(String path,String appendix)
    {
        boolean success;
        success = (new File(path+"models")).mkdir();
        success = (new File(path+"sims")).mkdir();
        success = (new File(path+"diffs")).mkdir();
        success = (new File(path+"masters")).mkdir();
                            
        IJ.saveAsTiff(stackp, path+"models/model_of_views"+appendix+".tif");
        IJ.saveAsTiff(stacksimp, path+"sims/simulation_of_views"+appendix+".tif");
        IJ.saveAsTiff(diffp, path+"diffs/difference"+appendix+".tif");
        //IJ.saveAsTiff(merged, path+"sum_of_models.tif");

        String log= IJ.getLog();
        File file = new File(path+"log.txt");
        master.write_xyz(appendix);
        try 
        {
        FileWriter fw = new FileWriter(file, false);
        PrintWriter pw = new PrintWriter(fw);
        pw.println(log);
        pw.close();
        } catch (IOException ex) 
        {
            System.out.println("could not save "+path+"log.txt");
            Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
        }       
        globalcounter+=1;
        
    }
    
    void smooth_structure()
    {
        for(Master.Atom atom: master.atoms)
        {
            atom.temp_z=atom.z*smoothweight;                   
            for(Master.Atom neighbor: atom.neighbors)
            {
                atom.temp_z+=neighbor.z;
               // System.out.println(neighbor.z);
            }
            atom.temp_z=atom.temp_z/(atom.neighbors.size()+smoothweight);
            //System.out.println(atom.neighbors.size());
            //System.out.println();
        }
        float[] shift;
        for(Master.Atom atom: master.atoms)
        {
            shift=new float[]{0.0f,0.0f,atom.temp_z-atom.z};
            //System.out.println(Arrays.toString(shift));
            atom.shift_coordinate(shift, false, false);
        }
        master.norm_z();
        master.update(true,true,true);
        
        push_pixels();
    }
    
    void smooth_structure2()
    {
        float[] vector;
        float x,y,z;
        int counter;
        master.update_rings();
        for(Master.Atom atom: master.atoms)
        {
            //if(atom.id>5)
               // break;
            if(atom.viewcounts<2)
                continue;
           // System.out.println("\n\nbefore: "+atom.toString());
            
            vector=new float[3];
            x=0;
            y=0;
            z=0;
            counter=0;
            for(Master.Ring ring: atom.rings)
            {
                //System.out.println(ring.toString());
                vector=modify_length(get_vector(ring.center,new float[] {atom.x,atom.y,atom.z}),ring.radius);
                //System.out.println(Arrays.toString(vector)+ "... center - projection");
                x+=ring.center[0]+(vector[0]);
                y+=ring.center[1]+(vector[1]);
                z+=ring.center[2]+(vector[2]);
                counter+=1;
                //System.out.println(ring.center[0]+(vector[0])+" "+ring.center[1]+(vector[1])+" "+ring.center[2]+(vector[2])+" [xyz] projection");
                
            }
            if (counter==0)
                continue;
            x=x/counter;
            y=y/counter;
            z=z/counter;
            //System.out.println(x+" "+y+" "+z+" [xyz] targetpoint");
            //System.out.println(Arrays.toString(new float[]{(x-atom.x)*scale,(y-atom.y)*scale,(z-atom.z)*scale}));
            atom.shift_coordinate(new float[]{(x-atom.x)*smoothscale,(y-atom.y)*smoothscale,(z-atom.z)*smoothscale}, false, false);
            //System.out.println("after: "+atom.toString());
            //System.out.println("deviation after: "+master.deviation);
        }
        master.norm_z();
        master.update(true, true, true);
        push_pixels();
    }
    
    //creates vector from p1 to p2
    float[] get_vector(float[] p1, float[] p2)
    {
        return new float[]{p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]};
    }
    //returns vector with certain length
    float[] modify_length(float[] vec, float len)
    {
        float l=0;
        for(float i : vec)
        {
            l+=i*i;
        }          
        l=(float)Math.sqrt(l);
        return new float[]{len*vec[0]/l,len/l*vec[1],len/l*vec[2]};
    }
    
    
    void add_models()
    {
        int area=modelSt.getHeight()*modelSt.getWidth();
        
        float[] array1=new float[area];
        float[] array2=new float[area];
        
       
        for(int s=1;s<=modelSt.getSize();++s)
        {
            array1 = ((float[])modelSt.getPixels(s));
            for(int i=0;i<area;i++)
            {
                if (array1[i]!=0)
                 array1[i]=2;
            }
            array2 = ((float[])stack_model.getPixels(s));
            float[] sum=new float[area];
            for(int i=0;i<area;i++)
            {
                sum[i]=array1[i]+array2[i];
            }
            
            mergedSt.setPixels(sum, s);
        }
    }
    

 
        
    
   void optimize_skewmatrix2() //only one round
    {
        System.out.println("deviation before optimizing skewmatrix "+master.mdeviation);
        System.out.println("merit before optimizing skewmatrix "+master.merit);

        float m11;
        float m12;
        float m21;
        float m22;
        
        float olddev_=master.mdeviation;  
        double oldcorr_=master.correlation;

        for(Master.View view:master.views)
        {
            if (view.slicenumber!=1)
                continue;
            float[] steps = new float[4];        
            for(int i=0; i < 4; ++i)
            {
                steps[i] = 0.01f * (random.nextBoolean()?1f:-1f);
            }
            float it=0;
            do
            {
                for(int k=0;k<4;++k)
                { 
                    int k1 = k%2;
                    int k2 = k/2;
                    
                    float m_old=view.skewmatrix.map[k1][k2];     
                    view.skewmatrix.map[k1][k2]+=steps[k];
                    view.skewmatrix.normalize();
                    System.out.println(view.skewmatrix.map[k1][k2]);
                    System.out.println(m_old);

                    master.update(true,switch_,false);
                    if(!switch_)
                    {
                        if(master.mdeviation<olddev_)
                        {
                            System.out.println("olddev: "+olddev_);
                            System.out.println("newdev: "+master.mdeviation);
                            System.out.println();
                            steps[k]*=1.4f;
                            olddev_=master.mdeviation;
                            //break;
                            it+=1;
                        }
                        else
                        {
                            view.skewmatrix.map[k1][k2]=m_old;
                            master.update(true,false,false);
                            //System.out.println("SM: "+view.skewmatrix.map[k1][k2]);
                            System.out.println("dev:\t"+master.mdeviation);
                            steps[k]*=-0.6f;
                        }
                    }
                    else
                    {
                        if(master.correlation<oldcorr_)
                        {
                            //System.out.println("olddev: "+olddev_);
                            //System.out.println("newdev: "+master.deviation);
                            //System.out.println();
                            steps[k]*=1.4f;
                            oldcorr_=master.correlation;
                        }
                        else
                        {
                            //view.skewmatrix.map=m_old.clone();
                            master.update(true,switch_,false);
                            steps[k]*=-0.6f;
                        }
                    }
                }
                it++;
                master.update(true,true,false);

            }
            while (it<10);
        }
        push_pixels();
        System.out.println("deviation after optimizing skewmatrix "+master.mdeviation);

        System.out.println("merit after optimizing skewmatrix "+master.merit);

        
    }
   
    void optimize_translation()
    {
        float[] steps = new float[2];        
        

        for (Master.View view:master.views)
        {
            for(int i=0; i < 2; ++i)
            {
                steps[i] = 1f * (random.nextBoolean()?1f:-1f);
            }
            do
            {
                for(int k=0;k<2;k++)
                {
                    double match=view.lcorrelation;
                    float oldoffset=view.offset[k];
                    view.offset[k]+=steps[k];
                    master.update(true,true,false);
                    if(show_moves)
                        push_pixels();
                    if (view.lcorrelation<match)
                    {
                        steps[k]*=1.1;
                        System.out.println("new merit\t"+master.merit);

                    }
                    else
                    {
                        view.offset[k]=oldoffset;
                        master.update(true,true,false);     
                        steps[k]*=-0.7;
                    }
                }
            }while(Math.abs(steps[0])>0.1 || Math.abs(steps[1])>0.1);
        }
    }
    /*
    void optimize_aberrations()
    {
        float[] step = new float[7];
        double match;
        float bestmatch=-1;
        double max;
        ArrayList<Integer> index=new ArrayList<Integer>();
        for(int i=0;i<7;++i)
        {
            index.add(i);
        }
        
        for (Master.View view:master.views)
        {
            for (int k=0 ; k<3;++k)
            {
                step[k]=2.0f;
            }
            for (int k=3 ; k<7;++k)
            {
                step[k]=200.0f;
            }
            int s=0;
            do
            {
                s+=1;
                for(int i:index)
                {        

                    for(int r=0;r<3;++r)
                    {
                        match=view.lcorrelation;
                        view.sim.beam.aberrations[i]+=step[i];
                        view.sim.beam.initBeam();
                        view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
                        //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                        view.update_correlation();                        //System.out.println(newmerit);
                        if(show_moves)
                            push_pixels();
                        if (view.lcorrelation<match)
                        {
                            step[i]*=1.1;
                        }
                        else
                        {
                            view.sim.beam.aberrations[i]-=step[i];
                            view.sim.beam.initBeam();
                            step[i]*=-0.7;
                        }
                    }
                }
                max=Math.abs(step[0]);
                for(int i =1; i<step.length-4;++i )
                {
                    if (Math.abs(step[i])>max)
                        max=Math.abs(step[i]);
                }
                for(int i =step.length-4; i<step.length;++i )
                {
                    if (Math.abs(step[i]*0.01)>max)
                        max=Math.abs(step[i]*0.01);
                }

            }while (max>0.1 && s<=4);
        }
        for(Master.View view:master.views)
        {
            view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
            view.update_correlation();
        }
        push_pixels();
        write_beam_param();
    }*/
    
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
        
        for (Master.View view:master.views)
        {
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
            do
            {

                for(int i=ind.size()-1;i>=0;--i)
                {
                    if(ind.get(i)==true)
                        continue;

                    for(int r=0;r<2;++r)
                    {
                        match=view.lcorrelation;
                        view.sim.beam.aberrations[i]+=step[i];
                        view.sim.beam.initBeam();
                        view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
                        //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                        view.update_correlation();                        //System.out.println(newmerit);
                        if(show_moves)
                            push_pixels();
                        if (view.lcorrelation<match)
                        {
                            ind.set(i, true);
                            break;
                        }
                        else
                        {
                            view.sim.beam.aberrations[i]-=step[i];
                            view.sim.beam.initBeam();
                            view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
                            //System.out.println(Arrays.toString((float[])impSt.getPixels(1)));
                            view.update_correlation();     
                            step[i]*=-1;
                        }
                    }
                    step[i]*=0.5;
                }
                
                s+=1;

            }while ( s<=1);
        }
        for(Master.View view:master.views)
        {
            view.sim.beam.initBeam();
            view.simarray=ImageCalculations.normalize_image(view.simulate_image(qstem),view.expIm);
            view.update_correlation();
        }
        push_pixels();
        write_beam_param();

    }

    void optimize_gaussianblur()
    {
        //master.update_deviations();
        double oldcorr=master.correlation;
        System.out.println("corr before: "+ oldcorr);
        float oldblur=0;
        int s=0;
        boolean change=false;
        for(Master.View view:master.views)
        {
            if (!view.atoms2D.isEmpty())
            {
                oldblur=view.sim.sigma;
                break;
            }
        }
        float step=0.2f ;//* (random.nextBoolean()?1f:-1f);
        do
        {
            for(Master.View view:master.views)
            {
                view.sim.sigma+=step;
            }
            master.update(false, true, false);
            if (master.correlation<oldcorr)
            {
                change=true;
                break;
            }            
            else
            {
                for(Master.View view:master.views)
                {
                    view.sim.sigma=oldblur;  
                }
                step*=-1;
            }
            s++;
        }while(s<2);
        System.out.println("sigma after: "+master.views.get(0).sim.sigma);
        if (change==false)
        {
            System.out.println("no change");
            master.update(false, true, false);
        }
        push_pixels();
        write_beam_param();
        System.out.println("corr after : "+ master.correlation);
    }
    
    
    
    //minimizes energy for proper analysis of the structure
    //only one round
    void optimize_fieldOfView()
    {
        master.update(false,false,true);
        float oldfov=fov;
        double firstenergy= master.energy;
        float step=1f * (random.nextBoolean()?1f:-1f);;
        do
        {
            fov+=step;
            master.change_fov(fov);
            master.update(false,false,true);
            //master.atoms.get(0).shift_coordinate(new float[]{0,0,0}, true, true);
            //System.out.println(master.energy);
            if(master.energy<firstenergy)
                break;
            else
            {
                fov=oldfov;
                master.change_fov(fov);
                step*=-0.8;             
            }
        }
        while(Math.abs(step)>=0.08);
        master.change_fov(fov);
        master.update(false,false,true);
        write_beam_param();
    }
    
    //used image correlation for optimizing
    //only one round
    void optimize_fieldOfView2()
    {
        master.update(false,true,false);
        float oldfov=master.fov;
        fov=master.fov;
        double oldcorr=master.correlation;
        float step=(float) (random.nextGaussian())/2f;;
        int iterator=0;
        do
        {
            fov+=step;
            master.change_fov(fov);
            //master.atoms.get(0).shift_coordinate(new float[]{0,0,0}, true, true);
            //System.out.println(master.energy);
            if(master.correlation<oldcorr)
                break;
            else
            {
                fov=oldfov;
                master.change_fov(fov);
                step*=-0.8;             
            }
            iterator+=1;
        }
        while(iterator<=3);
        
        write_beam_param();
    }
    
   
    
    void fix_tilt_tune3() //only one round
    {
        float x;
        float y;
        float step;
        for(Master.View view:master.views)
        {
            step=10.0f * (random.nextBoolean()?1f:-1f);
            do
            {
                if(!view.atoms2D.isEmpty())
                {
                    float angle=view.first_angle;
                    float[] vector= new float[]{0f,0f,0f};

                    master.update(true,false,false);
                    if (switch_)
                    {
                        view.update_model();
                        view.update_correlation();
                    }

                    float olddev_=master.mdeviation; 
                    double oldcorr_=master.correlation; 
                    //System.out.println("dev before "+master.deviation);

                    view.inplaneangle+=step;
                    x = (float) Math.cos(view.inplaneangle/180f*3.1415f);
                    y = (float) Math.sin(view.inplaneangle/180f*3.1415f);
                    vector[0] =x;
                    vector[1] =y;
                    //System.out.println(x+" "+y);
                    view.quaternion.angleToQuaternion(angle, vector);
                    //view.quaternion.normalize();


                    master.update(true,false,false);
                    if (switch_)
                    {
                        view.update_model();
                        view.update_correlation();
                    }

                    //show_and_wait();
                    //System.out.println(ang);
                    //System.out.println(master.deviation);
                    //System.out.println("olddev: "+olddev_);
                    //System.out.println("dev: "+master.deviation);
                    //System.out.println();

                    if(!switch_)
                    {
                        if(master.mdeviation>=olddev_)
                        {
                            //System.out.println("view: "+view.slicenumber);
                            //System.out.println("olddev: "+olddev_);
                            //System.out.println("newdev: "+master.deviation);
                            //System.out.println("best angle: "+ang);
                            //System.out.println();
                            view.inplaneangle-=step;   
                            x = (float) Math.cos(view.inplaneangle/180f*3.1415f);
                            y = (float) Math.sin(view.inplaneangle/180f*3.1415f);
                            vector= new float[]{0f,0f,0f};
                            vector[0] =x;
                            vector[1] =y;
                            view.quaternion.angleToQuaternion(angle, vector);
                            step*=-0.8;
                        }
                        else
                            break;
                    }       

                    else
                    {
                        if(master.correlation>=oldcorr_)
                        {
                            view.inplaneangle-=step; 
                            x = (float) Math.cos(view.inplaneangle/180f*3.1415f);
                            y = (float) Math.sin(view.inplaneangle/180f*3.1415f);
                            vector= new float[]{0f,0f,0f};
                            vector[0] =x;
                            vector[1] =y;
                            view.quaternion.angleToQuaternion(angle, vector);                           
                            step*=-0.8;

                        } 
                       else
                           break;
                    }
                }
            }while(Math.abs(step)>0.1);
            
            
        }
        master.update(true,true,false);
     
        push_pixels();
        
        System.out.println("merit after "+master.merit);
        for (Master.View view: master.views)
        {
            System.out.println("angle after "+view.inplaneangle+ " of view "+view.slicenumber);
        }
    }
    
    //tunes the tilt-angle from the microscope data
    //usually not necessary
    //don't use it too early
    void tune_fixtilt() //only one round
    {
        float x;
        float y;
        float step;
        for(Master.View view:master.views)
        {
            step=1.0f/180*3.1415f * (random.nextBoolean()?1f:-1f);
            do
            {
                if(!view.atoms2D.isEmpty())
                {
                    float[] vector= new float[3];

                    //master.update(true,false,false);

                    float olddev_=master.mdeviation; 
                    double oldcorr_=master.correlation; 
                    //System.out.println("dev before "+master.deviation);

                    view.first_angle+=step;
                    x = (float) Math.cos(view.inplaneangle/180f*3.1415f);
                    y = (float) Math.sin(view.inplaneangle/180f*3.1415f);
                    vector[0] =x;
                    vector[1] =y;
                    //System.out.println(x+" "+y);
                    view.quaternion.angleToQuaternion(view.first_angle, vector);
                    //view.quaternion.normalize();


                    master.update(true,false,false);
                    
                    //show_and_wait();
                    //System.out.println(ang);
                    //System.out.println(master.deviation);
                    //System.out.println("olddev: "+olddev_);
                    //System.out.println("dev: "+master.deviation);
                    //System.out.println();

                    if(!switch_)
                    {
                        if(master.mdeviation>=olddev_)
                        {
                            //System.out.println("view: "+view.slicenumber);
                            //System.out.println("olddev: "+olddev_);
                            //System.out.println("newdev: "+master.deviation);
                            //System.out.println("best angle: "+ang);
                            //System.out.println();
                            view.first_angle-=step;   
                            view.quaternion.angleToQuaternion(view.first_angle, vector);
                            step*=-0.8;
                            master.update(true,false,false);
                        }
                        else
                            break;
                    }       

                    else
                    {
                        view.update_model();
                        view.update_correlation();
                
                        if(master.correlation>=oldcorr_)
                        {
                            view.first_angle-=step; 
                            view.quaternion.angleToQuaternion(view.first_angle, vector);                           
                            step*=-0.8;

                        } 
                        else
                        {
                           System.out.println("new fixed angle");
                           break;
                        }
                    }
                }
            }while(Math.abs(step)>0.05);
            
            
        }
        master.update(true,true,false);
     
        push_pixels();
        System.out.println("dev after "+master.mdeviation);
        for (Master.View view: master.views)
        {
            System.out.println("(fixed) angle after "+view.first_angle+ " of view "+view.slicenumber);
        }
    }
    
   
    
    void show_and_wait()
    {
        push_pixels();
        try {
        Thread.sleep(2000);
        } catch(InterruptedException ex) {
            Thread.currentThread().interrupt();
        }
    }
    
    void modify_compression()
    {
        float olddev= master.mdeviation;
        for (Master.View view : master.views)
        {
            float step=0.01f * (random.nextBoolean()?1f:-1f);
            do
            {
                view.stretch_factor+=step;
                /*for(Master.Atom2D at:view.atoms2D)
                {
                    at.update2D();
                }*/
                master.update(true,false,false);
                if (master.mdeviation<olddev)
                {
                    olddev=master.mdeviation;
                    System.out.println("model stretched/compressed in view "+view.slicenumber+"!");
                    System.out.println("new stretchfactor:\t"+view.stretch_factor);
                    continue;
                }
                else
                {
                    view.stretch_factor-=step;
                    step*=-0.8;
                }
            }
            while (Math.abs(step)>0.0001);
                
        }
        master.update(true, true, false);
    }
    
    /*void modify_deltaintensity()
    {
        float olddev= master.mdeviation;
        for (Master.View view : master.views)
        {
            float step=1f * (random.nextBoolean()?1f:-1f);
            do
            {
                view.deltaintensity+=step;
                for(Master.Atom2D at:view.atoms2D)
                {
                    at.update2D();
                }
                master.update(false,true,false);
                System.out.println("merit before: "+olddev);
                System.out.println("merit after: "+master.mdeviation);
                if (master.mdeviation<olddev)
                {
                    olddev=master.mdeviation;
                    System.out.println("deltaintensity matched in view "+view.slicenumber+"!");
                    System.out.println("new deltaintensity multiplicator:\t"+view.deltaintensity);
                    continue;
                }
                else
                {
                    view.deltaintensity-=step;
                    step*=-0.8;
                }
            }
            while (Math.abs(step)>0.001);
                
        }
        master.update(false, true, false);
    }*/
    
    void drawBonds(ImagePlus bild)
    {
        Overlay overlay = new Overlay();
        for(Master.Bond bond:master.bonds)
        {
            
        }
        
        //overlay.crop(impSlice,impSlice);
        bild.setOverlay(overlay);
    }
    
    
    void write_top(String appendix)
    {
        try
        {
            File topfile ;
            topfile= new File(path+imp.getShortTitle() + "_new"+appendix+".top");
            FileWriter fw = new FileWriter(topfile, false);
            PrintWriter pw = new PrintWriter(fw);

            pw.println("#creator: " + getClass().getSimpleName() + " source: " + imp.getTitle());
            pw.println("#MASTER are the ATOMs in the joined topology with their id, view counts and approximate x,y,z coordinates\n" +
                               "#BONDs connect the two ATOMs with id1 and id2\n" +	
                               "#Every VIEW lists a LABEL, a TRANSLATION, a SKEWMATRIX,\n" +
                               "#a rotation QUATERNION as well as all ATOMS with their id, x and y in image pixels and z=0.0");
            pw.println("#Energy\t"+master.energy);
            pw.println("#Imageerror\t"+master.correlation);
            pw.println("#Distanceerror\t"+master.mdeviation);
            
            pw.println("\n\nMASTER\t" + master.atoms.size());
            for(Master.Atom atom: master.atoms)
            {
                if (atom.element==6)
                    pw.println("ATOM\t"+atom.id +"\t"+atom.viewcounts +"\t"+atom.x+"\t"+atom.y+"\t"+atom.z);
                else
                    pw.println("ATOM\t"+atom.id +"\t"+atom.viewcounts +"\t"+atom.x+"\t"+atom.y+"\t"+atom.z+"\t"+atom.element);
            }
            pw.print("\n\n");
            for(Master.Bond b: master.bonds)
            {
                pw.println("BOND\t"+b.a1.id +"\t"+b.a2.id );
            }
            pw.print("\n\n");
            for(Master.Ring r: master.rings)
            {
                pw.print("RING");
                for(Master.Atom a : r.atoms)
                    pw.print("\t"+a.id);
                pw.println();
            }
            pw.print("\n\n");
            for(int v = 0; v<stackNb; ++v)
            {
                if( (master.views.get(v).atoms2D.size() == 0))
                {	continue;}
                pw.println("\n\nVIEW\t" + v);
                pw.println("#view of slice:" + (v+1) + " in " + imp.getTitle());
                pw.println("LABEL\t"+master.views.get(v).label);
                pw.println("TRANSLATION\t"+master.views.get(v).translation[0]+
                        "\t"+master.views.get(v).translation[1]+
                        "\t"+master.views.get(v).translation[2]);
                pw.println("SKEWMATRIX\t"+master.views.get(v).skewmatrix.map[0][0]+"\t"+
                        master.views.get(v).skewmatrix.map[0][1]+"\t"+
                        master.views.get(v).skewmatrix.map[1][0]+"\t"+
                        master.views.get(v).skewmatrix.map[1][1]+"\t");
                pw.println("QUATERNION\t"+master.views.get(v).first_angle+"\t"+master.views.get(v).inplaneangle/180*3.1415f);
                pw.println("STRETCH\t"+master.views.get(v).stretch_factor);
                pw.println("BACKGROUND\t"+master.views.get(v).backgroundlevel);

                for(Master.Atom2D atom2d : master.views.get(v).atoms2D)
                {
                    pw.println("ATOM\t"+atom2d.atom.id+"\t"+atom2d.expx+"\t"+atom2d.expy+"\t#intensity\t"+(atom2d.intensity_factor));
                    pw.println("SEEN\t"+atom2d.atom.id+"\t"+atom2d.realx+"\t"+atom2d.realy+"\t#intensity\t"+(atom2d.intensity_factor));
                }

                for(Master.Minus1Atom atom : master.views.get(v).oneatoms)
                {
                    pw.println("ATOM\t-1\t"+atom.x+"\t"+atom.y+"\t#intensity\t"+(atom.intensity_factor));
                }
            }
            pw.close();
            IJ.log("rewrote " + topfile.getName());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
        save_tifdata(path);
    } 
    
     
    void print_background()
    {
        ImageStack temp = new ImageStack(impWidth, impHeight, stackNb);
        for(Master.View view:master.views)
        {
            float[] temparray= new float[impWidth*impHeight];
            for(Master.Minus1Atom atom1 : view.oneatoms)
            {
                temparray[atom1.x+atom1.y*impWidth]+=1;
            }
            temp.setPixels(temparray, view.slicenumber);           
        }
        ImagePlus tempp = new ImagePlus("background atoms", temp);
        tempp.show();
    }
    
    //shifts one atom and documents the differential of the deviation and energy
    void analyze_energy_deviation_contribution()
    {
        FileWriter fw = null ;
        File energyfile =new File(path+"energy_image_contribution.txt");
        try {
            float ang2px=impWidth/(fov*10);
            fw = new FileWriter(energyfile, true);
            PrintWriter pw = new PrintWriter(fw);
            pw.println("#analyzing energy and image contribution of an optimized model");
            pw.println("#selected atom are shifted in each direction (x,y and z) and their deviation of the energy and image correlation is saved");
            pw.println("#id\tshift\tdE\tdC");
            double oldenergy=master.energy;
            double oldcorr= master.correlation;
            for(int atomnumber=0;atomnumber<master.atoms.size();atomnumber++)
            {
                if(atomnumber!=8 && atomnumber!=602&&atomnumber!=898)
                    continue;
                Master.Atom atom=master.atoms.get(atomnumber);
                if(atom.viewcounts<=1)
                    continue;
                System.out.println("analyzing atom "+atomnumber);
                
                float[] shift= new float[3];
                float oldx =atom.x;
                float oldy =atom.y;
                float oldz =atom.z;
                float end=5*ang2px;
                for(float i=0.01f*ang2px;i<=end;i+=0.01f*ang2px)
                {                                  
                    System.out.println("radius [A]: "+i/ang2px);
                    for(float j=0f;j<3f;j+=0.5f)
                    {
                        if(j==Math.round(j))
                            shift[(int)j]=1*i;
                        else
                            shift[(int)j]=-1*i;
                        
                        atom.shift_coordinate(shift,true,true);
                        pw.print(atomnumber+"\t");
                        pw.print(shift[0]/ang2px+"\t"+shift[1]/ang2px+"\t"+shift[2]/ang2px+"\t");
                        pw.print((master.energy-oldenergy+"\t"));
                        pw.print((master.correlation-oldcorr));
                        pw.println();
                        shift[(int)j]=0;
                        if(show_moves)
                            push_pixels();
                        atom.set_coordinates(oldx,oldy,oldz,true,true);
                    }
                     if(IJ.escapePressed())
                        break;
                    
                }
                if(IJ.escapePressed())
                    break;
            }

        } catch (IOException ex) {
            Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fw.close();
            } catch (IOException ex) {
                Logger.getLogger(Reconstructor.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }

  
}



    
   


