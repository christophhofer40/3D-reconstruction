package christoph;
import ij.*;
import ij.ImageStack.*;
import ij.WindowManager.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.File;
import java.io.FileWriter;
import java.io.*;
import christoph.run_lammps.Run_lammps;
import christoph.Quaternion;
import christoph.Shear;
import christoph.Fakesimulation;
import christoph.QSTEM_simulator;
import christoph.ImageCalculations;
import java.util.logging.Level;
import java.util.logging.Logger;


public class Master
{
    public int atomnumber=0;
    public int seen_atoms;
    public float mdeviation=0;
    public double correlation=0;
    public double energy=0;
    public double merit = 0;
    public ArrayList<View> views = new ArrayList<View>(10);
    public ArrayList<Atom> atoms = new ArrayList<Atom>(500);
    public ArrayList<Ring> rings=new ArrayList<Ring>();
    public ArrayList<Bond> bonds=new ArrayList<Bond>();

    public float[] imarray;
    public float[][] relative_coordlist;
    public static Run_lammps lammps;
    public float factor=-1;
    public float fov;
    public int impWidth;
    public String xyzpath;
    public boolean qstem =false;

    public Master(int atomnumber, int impWidth, float fov)
    {
        this.atomnumber=atomnumber;
        this.impWidth=impWidth;
        this.fov=fov;
    }

    public void set_xyzPath(String path)
    {
        xyzpath=path;	
    }

    public void shift_coordinate(int atomid, float[] shift, boolean bcorrelation, boolean benergy)
    {
        this.atoms.get(atomid).shift_coordinate(shift, bcorrelation, benergy);
    }

    public void set_coordinates(int atomid, float x, float y, float z, boolean bcorrelation, boolean benergy)
    {
        this.atoms.get(atomid).set_coordinates(x, y, z,bcorrelation,benergy);
    }
    
    /*public void complete_update(boolean update2D, boolean correlation, boolean energy)
    {
        if(update2D)
        {
            for(Atom atom:atoms)
            {
                for(Atom2D atom2d:atom.exp_atoms2D)
                {
                    if (atom2d.atom.viewcounts==views.size())
                    {
                        atom2d.update_rotation();
                        atom2d.shear_2Datom();
                        atom2d.compress();
                    }
                }
            }
            
            for (View view: views)
            {
                view.translate_allatoms();
                view.update_deviation();
            }
        }
        
        if(energy)
        {
            if (factor>0)
            {
                for (Atom at : this.atoms)
                {
                    if(at.viewcounts==views.size())
                    {
                        if(at.lammpsid>0)
                            lammps.set_coordinate(at.lammpsid,at.x/impWidth*fov*10,at.y/impWidth*fov*10,at.z/impWidth*fov*10);
                        else
                            System.out.println("WARNING: Atom "+at.id+ " not included in lammps");
                    }            
                }
                lammps.set_coordinate(-2,0,0,0);
                lammps.requestEnergy();
            }
        }
        
        
        if(correlation)
        {
            for(View view:views)
            {
                view.update_pixels(true);
                view.update_correlation();
            }
        }
        
        if(energy)
        {
            if (factor>0)
                this.energy=lammps.getEnergy()/seen_atoms;
        }
        
        if(correlation||energy)
            update_merit();
    }*/
    
    public void update(boolean update2D, boolean correlation, boolean energy)
    {
        if(update2D)
        {
            for(Atom atom:atoms)
            {
                for(Atom2D atom2d:atom.exp_atoms2D)
                {
                    if (atom2d.atom.viewcounts==views.size())
                    {
                        atom2d.update_rotation();
                        atom2d.shear_2Datom();
                        atom2d.compress();
                        if (atom.id==0)
                            {atom2d.view.update_translation();}
                        //atom2d.view.update_translation2();
                        atom2d.translate_2Datom();
                        //atom2d.update_deviation();
                    }
                }
            }
            update_deviations();
        }
        
        if(energy)
        {
            if (factor>0)
            {
                for (Atom at : this.atoms)
                {
                    if(at.viewcounts==views.size())
                    {
                        if(at.lammpsid>0)
                            lammps.set_coordinate(at.lammpsid,at.x/impWidth*fov*10,at.y/impWidth*fov*10,at.z/impWidth*fov*10);
                        else
                            System.out.println("WARNING: Atom "+at.id+ " not included in lammps");
                    }            
                }
                lammps.requestEnergy();
            }
        }
        
        
        if(correlation)
        {
            for(View view:views)
            {
                view.update_pixels(true);
                view.update_correlation();
            }
        }
        
        if(energy)
        {
            if (factor>0)
                this.energy=lammps.getEnergy()/seen_atoms;
        }
        
        if(correlation||energy)
            update_merit();
    }

   
    
    public void norm_z()
    {
        float zoffset=atoms.get(0).z;
        for(int atomi=0;atomi<atomnumber;atomi++)
        {
            if (atoms.get(atomi).viewcounts==views.size())
                atoms.get(atomi).z-=zoffset;
        }
        for(View view : views)
        {
            view.update_translation();
            //view.update_translation2();
            view.update_pixels(true);
        }
    }
    public void update_deviations()
    {
        this.mdeviation=0;
        for(Atom atom:atoms)
        {
            atom.sumOfDistances=0;
            for(Atom2D atom2d :atom.exp_atoms2D)
            {
                atom2d.ldeviation=(float)(atom2d.expx-atom2d.realx)*(atom2d.expx-atom2d.realx)+(atom2d.expy-atom2d.realy)*(atom2d.expy-atom2d.realy);
                atom.sumOfDistances+=atom2d.ldeviation;
            }
            this.mdeviation+=atom.sumOfDistances;
        }
        for (View view: views)
        {
            view.update_deviation();
        }
    }

    public void update_correlations()
    {
        for(View view:views)
        {
            view.update_pixels(true);
            view.update_correlation();
        }
    }

    public void update_rings()
    {
        for (Ring ring : this.rings)
        {
            ring.update();
        }
    }

    public void change_fov(float fov)
    {
        this.fov=fov;
        for(View view:this.views)
        {
            view.sim.fov=fov;
            view.sim.beam.initBeam();
        }
        update(true,true,true);
    }
    
    float[] norm_z(float[] coordlist)
    {
        System.out.println(Arrays.toString(coordlist));
        float minz=0;
        for (float z:coordlist)
        {
            if(z<minz)
            {
                minz=z;
            }
        }
        for (float z:coordlist)
        {
            z-=minz;
        }
        return coordlist;
        
    }

    
    /*
    public void update_energy()
    {
        write_xyz();
        //lammps.doProcess();
        laamps.shift_coordinate()
        this.energy=lammps.getEnergy()/seen_atoms;
        //System.out.println(energy);
    }
    */
    public void init_lammps(String path)
    {     
        write_xyz();
        if(lammps==null)
        {
            System.out.println("Starting LAAMPS");
            lammps = new Run_lammps(path); 
            System.out.println("LAMMPS started");
        }
        else
        {
            update(false,false,true);
            return;
        }
        this.energy=lammps.energy/seen_atoms;        
    }

    public void update_merit()
    {
        merit=factor*energy+(1-factor)*correlation;
    }
    public void write_xyz() 
    {

        //invert image x,y for proper 3D view in jmol
        try
        {
            File xyzfile = new File(xyzpath+"master.xyz");
            FileWriter fw = new FileWriter(xyzfile);
            PrintWriter pw = new PrintWriter(fw);
            int counter=0;
            for(int i = 0; i < atoms.size(); ++i)
            {
                if(atoms.get(i).viewcounts==views.size())
                {
                    counter+=1;
                }
            }

            pw.println("" + counter);
            pw.println("series");
            char element='C';
            for(int i = 0; i < atoms.size(); ++i)
            {
                if(atoms.get(i).viewcounts==views.size())
                {
                    element='C';/*
                    for(Ring ring :atoms.get(i).rings)
                    {
                        if(ring.atoms.size()==5)
                        {
                            if(element=='C'||element=='N')
                            {
                                element='N';
                            }
                            else
                                element='O';
                        }
                        if(ring.atoms.size()==7||element=='B')
                        {
                            if(element=='C')
                            {
                                element='B';
                            }
                            else
                                element='O';
                        }
                        
                    }*/
                    float[] pos = new float[] {atoms.get(i).x,
                        atoms.get(i).y,atoms.get(i).z};

                    pw.println(atoms.get(i).element + " "+(pos[0]/impWidth*fov*10) + " " + 
                                (pos[1]/impWidth*fov*10) + " " + 
                                (pos[2]/impWidth*fov*10)	);
                }
            }
            pw.close();
        }
        catch(Exception e)
            {
                e.printStackTrace();
            }
    }

    public void write_xyz(String appendix) 
    {
        boolean success;
        success = (new File(xyzpath+"masters")).mkdir();

        //invert image x,y for proper 3D view in jmol
        try
        {
            File xyzfile = new File(xyzpath+"masters/master"+appendix+".xyz");
            FileWriter fw = new FileWriter(xyzfile);
            PrintWriter pw = new PrintWriter(fw);
            int counter=0;
            for(int i = 0; i < atoms.size(); ++i)
            {
                if(atoms.get(i).viewcounts==views.size())
                {
                    counter+=1;
                }
            }

            pw.println("" + counter);
            pw.println("series");
            char element='C';
            for(int i = 0; i < atoms.size(); ++i)
            {
                if(atoms.get(i).viewcounts==views.size())
                {
                    element='C';/*
                    for(Ring ring :atoms.get(i).rings)
                    {
                        if(ring.atoms.size()==5)
                        {
                            if(element=='C'||element=='N')
                            {
                                element='N';
                            }
                            else
                                element='O';
                        }
                        if(ring.atoms.size()==7||element=='B')
                        {
                            if(element=='C')
                            {
                                element='B';
                            }
                            else
                                element='O';
                        }
                        
                    }*/
                    float[] pos = new float[] {atoms.get(i).x,
                        atoms.get(i).y,atoms.get(i).z};

                    pw.println(atoms.get(i).element + " "+(pos[0]/impWidth*fov*10) + " " + 
                                (pos[1]/impWidth*fov*10) + " " + 
                                (pos[2]/impWidth*fov*10));
                }
            }
            pw.close();
        }
        catch(Exception e)
            {
                e.printStackTrace();
            }
    }



    public class Atom
    {
        public float x;
        public float y;
        public float z;
        public int id;
        public int lammpsid;
        public double sumOfStats=0;
        public int viewcounts;
        public float sumOfDistances=0;
        public ArrayList<Atom> neighbors = new ArrayList<Atom>(3);
        public ArrayList<Atom2D> exp_atoms2D = new ArrayList<Atom2D>(10);
        public ArrayList<Ring> rings =  new ArrayList<Ring>(3);
        public float temp_z;
        public float[] force=new float[]{Float.NaN,Float.NaN,Float.NaN}; 
        public float[] prefered_dir=new float[3];
        public int element=6;




        public Atom(float x, float y, float z, int id, int viewcounts)
        {
            this.x=x;
            this.y=y;
            this.z=z;
            this.id=id;    
            this.viewcounts=viewcounts;
        }

        public void shift_coordinate(float[] shift, boolean bcorrelation, boolean benergy)
        {
            this.set_coordinates(this.x+shift[0],this.y+shift[1],this.z+shift[2],bcorrelation, benergy);
        }
        /*
        public void set_coordinates(float x, float y, float z, boolean bcorrelation, boolean benergy)
        {
            if (viewcounts==views.size())
            {
                this.x=x;
                this.y=y;
                this.z=z;
                if (benergy)
                {
                    lammps.set_coordinate(this.lammpsid,x/impWidth*fov*10,y/impWidth*fov*10,z/impWidth*fov*10);//convert pixel in Angstrom
                    lammps.requestEnergy();
                }
                this.update2D(); 
                if(bcorrelation)
                    update_correlations();
                //if(benergy)
                //    update_energy();
                if (benergy)
                    energy=lammps.getEnergy()/seen_atoms;
                update_merit();
            }
        }*/
        
        public void set_coordinates(float x, float y, float z, boolean bcorrelation, boolean benergy)
        {
            if (viewcounts==views.size())
            {
                this.x=x;
                this.y=y;
                this.z=z;
                //if (factor==0)
                //    benergy=false;
                if (benergy)
                {
                    lammps.set_coordinate(this.lammpsid,x/impWidth*fov*10,y/impWidth*fov*10,z/impWidth*fov*10);//convert pixel in Angstrom
                    lammps.requestEnergy();
                }
                update(true,bcorrelation,false);//energy already calculated
                //if(benergy)
                //    update_energy();
                if (benergy)
                    energy=lammps.getEnergy()/seen_atoms;
                update_merit();
            }
        }

        public float sum_distances()
        {
            float sum=0;
            for(Atom2D atom2d:exp_atoms2D)
            {
                sum+=atom2d.ldeviation;
            }
            return sum;
        }
        
        public float sum_newdistances(float[] shift)
        {
            float sum=0;
            float oldx=x;
            float oldy=y;
            float oldz=z;
            
            set_coordinates(x+shift[0],y+shift[1],z+shift[2],false,false);
            sum=sum_distances();
            set_coordinates(oldx,oldy,oldz,false,false);
            return sum;
        }

        public void update2D()
        {
            for(Atom2D atom2d:exp_atoms2D)
            {
                if (atom2d.atom.viewcounts==views.size())
                {
                    atom2d.update_rotation();
                    atom2d.shear_2Datom();
                    atom2d.compress();
                    if (id==0)
                        {atom2d.view.update_translation();}
                    //atom2d.view.update_translation2();
                    atom2d.translate_2Datom();
                    atom2d.update_deviation();
                }
            }
            /*if(switch_)
            {
                for(View view:master.views)
                {
                    view.update_pixels(true);
                    view.update_correlation();
                }
            }*/
            // this is too time consuming for each atom, only for shifting
            //  one atom needed           
        }
        
        public float get_avgBondlength()
        {
            float bl=0;
            int counter=0;
            for ( Atom n: neighbors)
            {
                if(n.viewcounts==views.size())
                {
                    bl+=(float)Math.sqrt((this.x-n.x)*(this.x-n.x)+(this.y-n.y)*(this.y-n.y)+(this.z-n.z)*(this.z-n.z));
                    counter+=1;
                }                        
            }
            if (counter>=3)
                return bl/(neighbors.size());
            else
                return -1;
        }

        public float[] get_centerDirection(float length) //1=full length
        {
            if (neighbors.size()<3)
                return new float[]{Float.NaN,Float.NaN,Float.NaN};
            else
            {
                float Ox=0;
                float Oy=0;
                float Oz=0;
                
                for (Atom neighbor:neighbors)
                {
                    Ox+=neighbor.x;
                    Oy+=neighbor.y;
                    Oz+=neighbor.z;
                }
                Ox/=3;
                Oy/=3;
                Oz/=3;
                return new float[]{(Ox-this.x)*length,(Oy-this.y)*length,(Oz-this.z)*length};
            }
        }
        
        

        @Override
        public String toString()
        {
            return "Atom has coordinates "+this.x+" "+this.y+" "+this.z+" [xyz] ";
        }
    }

    public class Atom2D
    {
        public float realx;
        public float realy;
        public float realz;
        public float expx;
        public float expy;
        public float ldeviation=0;
        public float localstat=0;
        public float intensity_factor=1;
        public Atom atom;
        public View view;



        public Atom2D(float expx, float expy, Atom atom, View view)
        {
            this.expx=expx;
            this.expy=expy;
            this.atom=atom;
            this.view=view;

            if(atom.viewcounts==views.size())
            {
                update_rotation();
                shear_2Datom();
                compress();
                if (atom.id==0)
                    view.update_translation();
                //view.update_translation2();
                translate_2Datom();
                update_deviation();                   
            }
        }

        public void update2D()
        {
            update_rotation();
            shear_2Datom();
            compress();
            if (atom.id==0)
                view.update_translation();
            //view.update_translation2();
            translate_2Datom();
            update_deviation(); 
        }
        
        void update_rotation()
        {
            float[] xyz = view.quaternion.rotate(new float[]{
                atom.x-atoms.get(0).x,
                atom.y-atoms.get(0).y,
                atom.z});
            realx=xyz[0]+atoms.get(0).x;
            realy=xyz[1]+atoms.get(0).y; 
            realz=xyz[2];

        }

        void compress()
        {
            float dx=this.realx-impWidth/2;
            float dy=this.realy-impWidth/2;
            this.realx+=(view.stretch_factor-1)*dx;
            this.realy+=(view.stretch_factor-1)*dy;
        }

        void update_deviation()
        {
            if (Float.isNaN(this.expx) ||Float.isNaN(this.expx) )
                return;
            mdeviation-=this.ldeviation;
            this.ldeviation=(float)(expx-realx)*(expx-realx)+(expy-realy)*(expy-realy);
            mdeviation+=this.ldeviation;
            /*
            master.deviation=-deviation;
            this.deviation=(float)Math.sqrt((expx-realx)*(expx-realx)+(expy-realy)*(expy-realy));
            master.deviation+=deviation;
            */

        }

        void shear_2Datom()
        {
            float x = this.realx*view.skewmatrix.map[0][0]+this.realy*view.skewmatrix.map[0][1];
            float y = this.realx*view.skewmatrix.map[1][0]+this.realy*view.skewmatrix.map[1][1];
            this.realx=x;
            this.realy=y;
        }

        void translate_2Datom()
        {
            this.realx+=view.translation[0]+view.translation2[0];
            this.realy+=view.translation[1]+view.translation2[1];
        }

        public void get_local_stats(ImagePlus imp, int radius)
        {
            int x=(int)(atom.exp_atoms2D.get(view.slicenumber-1).realx+0.5);
            int y=(int)(atom.exp_atoms2D.get(view.slicenumber-1).realy+0.5);
            imp.setRoi(x-radius, y-radius, 2*radius+1, 2*radius+1);
            imp.setSlice(view.slicenumber);
            localstat=(float) imp.getStatistics().mean;
        }

        public float distance2D(Atom2D other)
        {
            float x = (this.realx-other.realx)*(this.realx-other.realx);
            float y = (this.realy-other.realy)*(this.realy-other.realy);
            return (float)Math.sqrt(x+y);
        }

    }


    public class View
    {
        public ArrayList<Atom2D> atoms2D = new ArrayList<Atom2D>(500); 
        public ArrayList<Minus1Atom> oneatoms =new ArrayList<Minus1Atom>(300);
        public Quaternion quaternion;
        public String label;
        public float[] translation = new float[3];
        public float[] translation2 = new float[3];
        public float[] alltranslate=new float[3];
        
        public Shear skewmatrix=new Shear();
        public int slicenumber;
        public int atom0index;
        public float[] imarray=new float[impWidth*impWidth];
        public float[] imarray2=new float[impWidth*impWidth];
        
        public float[] simarray=new float[impWidth*impWidth];
        public float[] diffarray=new float[impWidth*impWidth];//((squared))
        public float[] diffarray2=new float[impWidth*impWidth];
        public float[] expIm=new float[impWidth*impWidth];
        public float[] directionVector=new float[3];
        public float first_angle=0;
        public float inplaneangle;
        public Fakesimulation sim;
        public QSTEM_simulator qsim; 
        public double lcorrelation=0;
        public float stretch_factor=1;
        public float deviation=0;
        //public float deltaintensity=1;
        public float backgroundlevel=1;
        public float[] offset=new float[2]; // for better matching of translation
        
        
        public ArrayList<float[]> relative_coordlist = new ArrayList<float[]> ();


        public View(int slicenumber, float[] expIm)
        {
            this.slicenumber=slicenumber;
            qsim= new QSTEM_simulator(fov*10,impWidth,xyzpath,"qstem.qsc");
            sim=new Fakesimulation(imarray2,fov, impWidth);
            this.expIm=ImageCalculations.normalize_image(expIm,expIm);
        }

        public void update_pixels(boolean correlation)
        {
            
            Arrays.fill(imarray, backgroundlevel);
            float px;
            float py;
            if(!atoms2D.isEmpty())
            {
                for(Atom2D atom2D : atoms2D )  
                {
                    float intensity=(float) Math.pow(atom2D.atom.element,1.6);// MAADF
                    
                    intensity*=atom2D.intensity_factor;//should correct intensity fluctuations
                    if(atom2D.atom.viewcounts==views.size())
                    {
                        px= atom2D.realx;
                        py= atom2D.realy;
                        
                    }
                    else
                    {
                        px= atom2D.expx;
                        py= atom2D.expy;
                    }
                     
                    if((int)px>0.5 && (int)(px+1)<impWidth  && ((int)py>0.5 && (int)(py+1)<impWidth ))
                    {
                        imarray[(int)px+(int)py*impWidth ]+=(1-(px-(int)px))*(1-(py-(int)py))*intensity;
                        imarray[(int)px+1+(int)py*impWidth ]+=(px-(int)px)*(1-(py-(int)py))*intensity;
                        imarray[(int)px+1+((int)py+1)*impWidth ]+=(px-(int)px)*(py-(int)py)*intensity;
                        imarray[(int)px+((int)py+1)*impWidth ]+=(1-(px-(int)px))*(py-(int)py)*intensity;
                    }
                }
                
                
                for(Minus1Atom atom1 : oneatoms)
                {
                    if((int)atom1.x>0.5 && (int)(atom1.x+1)<impWidth  && ((int)atom1.y>0.5 && (int)(atom1.y+1)<impWidth ))
                        imarray[(int)(atom1.x+atom1.y*impWidth) ]+=(atom1.intensity_factor*(float) Math.pow(6,1.6));
                }
                if(correlation)
                    update_model();
            }
        }

        public void update_model()
        {
            imarray2=imarray.clone();
            /*Arrays.fill(imarray2, this.backgroundlevel);
            
            float px;
            float py;
            if(!atoms2D.isEmpty())
            {
                
                for(Atom2D atom2D : atoms2D )
                {
                    float intensity=(float) Math.pow(atom2D.atom.element,1.6);// MAADF
                    intensity*=atom2D.intensity_factor;//should correct intensity fluctuations
                    if(atom2D.atom.viewcounts==views.size())
                    {
                        px= atom2D.realx;
                        py= atom2D.realy;
                    }
                    else
                    {
                        px= atom2D.expx;
                        py= atom2D.expy;
                    }
                        

                    if((int)px>0.5 && (int)(px+1)<impWidth  && ((int)py>0.5 && (int)(py+1)<impWidth ))
                    {
                        imarray2[(int)px+(int)py*impWidth ]+=(1-(px-(int)px))*(1-(py-(int)py))*intensity;
                        imarray2[(int)px+1+(int)py*impWidth ]+=(px-(int)px)*(1-(py-(int)py))*intensity;
                        imarray2[(int)px+1+((int)py+1)*impWidth ]+=(px-(int)px)*(py-(int)py)*intensity;
                        imarray2[(int)px+((int)py+1)*impWidth ]+=(1-(px-(int)px))*(py-(int)py)*intensity;
                    }

                  
                        
                    
                }
                
                   
                for(Minus1Atom atom1 : oneatoms)
                {
                    px=atom1.x;
                    py=atom1.y;
                    float intensity=(float) Math.pow(6,1.6);// MAADF
                    intensity*=atom1.intensity_factor;//should correct intensity fluctuations
                    
                    if((int)px>0.5 && (int)(px+1)<impWidth  && ((int)py>0.5 && (int)(py+1)<impWidth ))
                    {
                        imarray2[(int)px+(int)py*impWidth ]+=(1-(px-(int)px))*(1-(py-(int)py))*intensity;
                        imarray2[(int)px+1+(int)py*impWidth ]+=(px-(int)px)*(1-(py-(int)py))*intensity;
                        imarray2[(int)px+1+((int)py+1)*impWidth ]+=(px-(int)px)*(py-(int)py)*intensity;
                        imarray2[(int)px+((int)py+1)*impWidth ]+=(1-(px-(int)px))*(py-(int)py)*intensity;//6 for carbon
                    }
                }
               
                if( deltaintensity!=1)
                {
                    for( float val : imarray2)
                    {
                        if (val!=0)
                        {
                            val*=deltaintensity;
                        }
                    }
                }
            */ 
               
                //imarray2=imarray.
                simarray=ImageCalculations.normalize_image(this.simulate_image(qstem),expIm);
                diffarray=calc_difference();
                diffarray2=calc_difference2();
            //} 
             

        }
        
        public void update_relative_coordlist()
        {
            float px;
            float py;
            float pz;
            relative_coordlist.clear();
            if(!atoms2D.isEmpty())
            {
                
                for(Atom2D atom2D : atoms2D )
                {
                    if(atom2D.atom.viewcounts==views.size())
                    {
                        px= atom2D.realx;
                        py= atom2D.realy;
                        pz= atom2D.realz;

                        relative_coordlist.add(new float[] {px/impWidth,py/impWidth,pz/impWidth});
                    }
                    else
                    {
                        px= atom2D.expx;
                        py= atom2D.expy;

                        relative_coordlist.add(new float[] {px/impWidth,py/impWidth,0});
                    }
                }
                
                   
                for(Minus1Atom atom1 : oneatoms)
                {
                    px=atom1.x;
                    py=atom1.y;
                    relative_coordlist.add(new float[] {px/impWidth,py/impWidth,0});
                }
            }  
            
        }
        

        public float[] calc_difference()//squared difference
        {
            float[] ret= new float[impWidth *impWidth ];
            for (int i=0;i<impWidth *impWidth ;++i)
            {
                ret[i]=(simarray[i]-expIm[i])*(simarray[i]-expIm[i]);
            }
            return ret;
        }
        public float[] calc_difference2()
        {
            float[] ret= new float[impWidth *impWidth ];
            for (int i=0;i<impWidth *impWidth ;++i)
            {
                ret[i]=(simarray[i]-expIm[i]);
            }
            return ret;
        }


        public void update_correlation()
        {
            correlation-=this.lcorrelation;
            //this.correlation=(float)img_match(simarray,expIm);
            this.lcorrelation=(double)ImageCalculations.img_match2(simarray,expIm);
            correlation+=this.lcorrelation;
        }
        
        //align to atom with id=0
        public void update_translation() 
        {
            if(!atoms2D.isEmpty())
            { 
                translation[0]=(atoms2D.get(atom0index).expx-atoms2D.get(atom0index).realx)+this.offset[0];
                translation[1]=(atoms2D.get(atom0index).expy-atoms2D.get(atom0index).realy)+this.offset[1];
                /*
                for(Atom2D atom2d: atoms2D)
                {
                    atom2d.realx+=translation[0];
                    atom2d.realy+=translation[1];
                }
                */
                //System.out.println(Arrays.toString(translation));
            }
        }
        
        //minimize deviation by translating
        public void update_translation2()
        {
            translation2[0]=0;
            translation2[1]=0;
            if(!atoms2D.isEmpty())
            { 
                for(Atom2D atom2d: atoms2D)
                {
                    if (Float.isNaN(atom2d.expx) || Float.isNaN(atom2d.realx)||
                            Float.isNaN(atom2d.expy)||Float.isNaN(atom2d.realy))
                        continue;
                    translation2[0]+=atom2d.expx-atom2d.realx;
                    translation2[1]+=atom2d.expy-atom2d.realy;
                }
                translation2[0]=translation2[0]/atoms2D.size();
                translation2[1]=translation2[1]/atoms2D.size();
                //System.out.println(Arrays.toString(translation2));
            }
        }
        
       /* public void translate_allatoms()
        {
            alltranslate[0]=0;
            alltranslate[1]=0;
            if(!atoms2D.isEmpty())
            { 
                for(Atom2D atom2d: atoms2D)
                {
                    if(Float.isNaN(atom2d.expx)||Float.isNaN(atom2d.expy))
                        continue;
                    alltranslate[0]+=atom2d.expx-atom2d.realx;
                    alltranslate[1]+=atom2d.expy-atom2d.realy;
                }
                alltranslate[0]=alltranslate[0]/atoms2D.size();
                alltranslate[1]=alltranslate[1]/atoms2D.size();
                //System.out.println(Arrays.toString(translation2));
            }
            
            for(Atom2D atom2d: atoms2D)
            {
                atom2d.realx+=alltranslate[0];
                atom2d.realy+=alltranslate[1];
            }
            
        }*/
        
        public void update_deviation()
        {
            //mdeviation-=this.deviation;
            this.deviation=0;
            for(Atom2D atom2d: atoms2D)
            {
                if (Float.isNaN(atom2d.expx) ||Float.isNaN(atom2d.expx) )
                continue;
                atom2d.ldeviation=(float)(atom2d.expx-atom2d.realx)*(atom2d.expx-atom2d.realx)+(atom2d.expy-atom2d.realy)*(atom2d.expy-atom2d.realy);
                this.deviation+=atom2d.ldeviation;
            }
            //mdeviation+=this.deviation;
        }

        public float[] simulate_image(boolean qstem)
        {
            if(qstem)
            {
                update_relative_coordlist();
                qsim.qsc.change_aberrations(this.sim.beam.aberrations);
                qsim.cfg.write_cfg(relative_coordlist,true);
                return qsim.simulate_image(sim.sigma);
            }
            else
            {
                return sim.simulate_image(this.imarray2);
            }
        }



        public void get_directionVector()
        {
            float[] z=new float[]{0,0,1};
            directionVector=quaternion.rotate(z);
            //System.out.println(quaternion.x0+ " "+quaternion.x1+ " "+quaternion.x2+ " "+quaternion.x3);
            //System.out.println(Arrays.toString(directionVector));
        }      
        public void apply_shear(boolean forward)
        {
            if(!atoms2D.isEmpty())
            {
                if(forward==true)
                {
                    for(Atom2D atom2d : atoms2D)
                    {
                        float x = atom2d.realx*skewmatrix.map[0][0]+atom2d.realy*skewmatrix.map[0][1];
                        float y = atom2d.realx*skewmatrix.map[1][0]+atom2d.realy*skewmatrix.map[1][1];
                        atom2d.realx=x;
                        atom2d.realy=y;
                    }
                }
                else
                {
                    for(Atom2D atom2d : atoms2D)
                    {
                        float x = atom2d.realx*skewmatrix.map[1][1]-atom2d.realy*skewmatrix.map[0][1];
                        float y = -atom2d.realx*skewmatrix.map[1][0]+atom2d.realy*skewmatrix.map[0][0];
                        atom2d.realx=x;
                        atom2d.realy=y;
                    }
                }
            }
        }

    }

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

    public class Ring
    {
        public ArrayList<Atom> atoms;
        public float[] center= new float[3];
        public int length;
        public float radius;

        public Ring(ArrayList<Atom> atoms)
        {
            this.atoms=(ArrayList<Atom>) atoms;
            for(Atom at : atoms)
            {
                at.rings.add(this);
            }
            this.length= atoms.size();
            rings.add(this); 
            update();
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
                z+=atom.z;
                counter+=1;
                //System.out.println(atom.toString());
            }
            center[0]=x/counter;
            center[1]=y/counter;
            center[2]=z/counter;

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
                dz=(atom.z-center[2])*(atom.z-center[2]);
                counter+=1;
                radius+=(float)Math.sqrt(dx+dy+dz);
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

    }

    public class Minus1Atom
    {
        public int x;
        public int y;
        public float intensity_factor=1;

        public Minus1Atom(int x, int y)
        {
            
            this.x=x;
            this.y=y;
        }
    }
    

}
