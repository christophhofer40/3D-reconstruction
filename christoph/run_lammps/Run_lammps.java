package christoph.run_lammps;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.OutputStream;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author christoph
 */
public class Run_lammps
{
    //path to xyz model
    String path;
    String thispath;
    Process process=null;
    BufferedInputStream is;
    InputStreamReader isr;
    BufferedReader br;
    OutputStream os;
    BufferedOutputStream bos;
    PrintWriter pwos;
    boolean hasEnergy=false;
    //OutputStreamWriter osw;
    //String npath;
    public double energy;
    public Run_lammps(String path)
    {
        thispath=this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"christoph/run_lammps/";
        this.path=path;
        energy=do_Process();
    }
    
    public void setPath(String path)
    {
        this.path=path;
    }
    
    
    public double do_Process()
    {
        try 
        {
            xyzToData(new File(path+"master.xyz"));
        } catch (IOException ex) {
            Logger.getLogger(Run_lammps.class.getName()).log(Level.SEVERE, null, ex);
        }
        run_lammps();
        return energy;
    }
    public void requestEnergy()
    {
        if(hasEnergy)
            System.out.println("WARNING: dublicated energy request");
        pwos.println("run 0");
        pwos.println("variable teng equal c_eatoms");
        pwos.println("print \"%% Energy = ${teng}\"");
        pwos.flush(); 
        hasEnergy=true;
    }
    public double getEnergy()
    {        
        try 
        {
            if (hasEnergy)
            {
                while(true)
                {
                    String line2=br.readLine();
                    //System.out.println(line2);
                    if(line2.startsWith("%% Energy"))
                    {
                        String ss[] = line2.split(" ");
                        energy = Double.parseDouble(ss[ss.length-1]);
                        while(br.ready()){
                            br.readLine();
                        }
                        break;
                    }
                    if(line2.startsWith("ERROR"))
                    {
                        //System.out.println(line2);
                        while(br.ready()){
                            br.readLine();
                        }
                        break;
                    }
                }
            }
            else
            {
                System.out.println("WARNING: reading energy second time");
                
            }
        } catch (IOException ex) {
            Logger.getLogger(Run_lammps.class.getName()).log(Level.SEVERE, null, ex);
        }
            
        hasEnergy=false;
        return this.energy;
    }

    void run_lammps()           
    {
        if(process==null)
        {
            try
            {
                process = new ProcessBuilder(thispath+"src/lmp_serial2", "-log","v2.log").start();
                is = new BufferedInputStream(process.getInputStream());
                isr = new InputStreamReader(is);
                br = new BufferedReader(isr);
                os=process.getOutputStream();
                bos=new BufferedOutputStream(os,1000000000);
                pwos=new PrintWriter(bos);
                
                File f = new File(path+"in.grGB");
                if(!f.exists() || !f.canRead() )
                    System.out.println("cannot read " + f.getPath());
                System.out.println("reading "+f.getPath());
        
                Scanner s= new Scanner(f);
                while(s.hasNextLine())
                {
                    String line1 = s.nextLine();
                    if(line1.length()>0 && !line1.startsWith("#"))
                    {
                        if(line1.startsWith("read_data"))
                        {
                            pwos.println("read_data "+path+"modelmaster.data");
                        }
                        else
                        {
                            pwos.println(line1);
                        }
                        pwos.flush();
                        System.out.println("cmd: " + line1);
                        
                    }
                    
                }
                System.out.println("reading energy");
                //pwos.println("run 0");
                //pwos.println("dummy");
                //pwos.flush();
                while(true)
                {
                    String line2=br.readLine();
                    System.out.println(line2);
                    if(line2.startsWith("%% Energy"))
                    {
                       String ss[] = line2.split(" ");
                       energy = Double.parseDouble(ss[ss.length-1]);
                       System.out.println("done1 ... energy = "+energy); 
                       while(br.ready()){
                           line2=br.readLine();                         
                           System.out.println(line2);
                       }
                       break;
                    }    
                    if(line2.startsWith("ERROR"))
                    {
                        //System.out.println(line2);
                        while(br.ready()){
                           br.readLine();
                       }
                        break;
                    }
                }
                System.out.println("done2 ... energy = "+energy);                
            }
            catch(Exception e)
            {
                e.printStackTrace();
            }
            hasEnergy=false;
        }
    }

    public void set_coordinate(int id, float x, float y, float z)     
    {
        try
        {
            if (id==-2)
            {
                pwos.flush();
                
            }
            else
            {
                pwos.println("group atomi id "+(id+1));
                //System.out.println(id+1);
                pwos.println("set atom "+id+" x "+x+ " y "+y+ " z "+z);
                pwos.flush();
            }
            //pwos.println("dummy");

            /*
            while(true)
            {
                String line2=br.readLine();
                //System.out.println(line2);
                if(line2.startsWith("%% Energy"))
                {
                   String ss[] = line2.split(" ");                      
                   energy = Double.parseDouble(ss[ss.length-1]);
                   while(br.ready()){
                       br.readLine();
                   }
                   break;
                }    
                if(line2.startsWith("ERROR"))
                {
                    System.out.println(line2);
                    while(br.ready()){
                       br.readLine();
                    }
                    break;
                }
            }*/           
        } 
        catch (Exception ex) 
        {
            Logger.getLogger(Run_lammps.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    

    
    /*
    void read_log() throws FileNotFoundException
    {
        File f = new File(thispath.replace("plugins/christoph/run_lammps/","")+"log.lammps");
        //FileReader fr= new FileReader(f);
        final Scanner s= new Scanner(f).useDelimiter(" ");
        while(s.hasNextLine())
        {    
            String line = s.nextLine();
            if(line.startsWith("%%"))
            {
                Scanner ls = new Scanner(line.replace(";", ""));
                ls.next();
                ls.next();
                ls.next();
                energy=ls.nextDouble();
            }
        }
    }
    */
    
    void xyzToData(File xyz) throws IOException
    {
      
        FileWriter fw = new FileWriter(path+"modelmaster.data",false);     
        PrintWriter pw= new PrintWriter(fw);
        final Scanner s= new Scanner(xyz).useDelimiter(" ");
        pw.println("LAMMPS data file. CGCMM style. atom_style atomic generated by Run_lammps.java");
        String el;
        float x,y,z;
        String line = s.nextLine();
        Scanner ls= new Scanner(line);
        pw.println(" "+ls.nextInt()+ " atoms");
        pw.println(" 0 bonds");
        pw.println(" 0 angles");
        pw.println(" 0 dihedrals");
        pw.println(" 0 impropers");
        pw.println(" 1 atom types");
        pw.println(" 0 bond types");
        pw.println(" 0 angle types");
        pw.println(" 0 dihedral types");
        pw.println(" 0 improper types");
        pw.println(" -10 200  xlo xhi");
        pw.println(" -0.500000 200  ylo yhi");
        pw.println(" -100 100  zlo zhi");
        pw.println("");
        pw.println("# Pair Coeffs");
        pw.println("#");
        pw.println("# 1  C");
        pw.println("");
        pw.println(" Masses");
        pw.println("");
        pw.println(" 1 12.010700 # C");
        pw.println("");
        pw.println("Atoms # atomic");
        pw.println("");
        s.nextLine();
        int counter=1;
        while(s.hasNextLine())
        {           
            line = s.nextLine();
            ls= new Scanner(line);
            el=ls.next();
            x=ls.nextFloat();
            y=ls.nextFloat();
            z=ls.nextFloat();
            pw.println(counter+" 1 "+x+" "+y+" "+z+" # "+el);
            counter+=1;
        }
        pw.close();
    }
    
}
