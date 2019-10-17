package christoph;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

public class Read_img
{
    public Read_img()
    {
        
    }
    /*
    Float[] read_img(String imgFilename)           
    {
        FileReader fr=null;
        try 
        {
            String imgpath=this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"qstem_images/";
            File imgFile = new File(imgpath+imgFilename);
            fr = new FileReader(imgFile);
        } 
        catch (FileNotFoundException ex) 
        {
            Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
        } 
        

        
        
        
        int[] integer_array=new int[8];
        // String[] integer_names= new String[]{"header_size","param_size", "comment_size", 
        //                                "Nx","Ny","is_complex", "data_size", "version"};
        double[] double_array=new double[3];
        // String[] double_names= new String[]{"t","dx","dy","aux_data"};
        
        String comment;
        for (int i=0;i<8;++i)
        {
            try {
                char[] char_array=new char[4];
                fr.read(char_array, 0,4);
                System.out.println(Arrays.toString(char_array));
                integer_array[i]=Integer.parseInt(String.valueOf(char_array));
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        for (int i=0;i<4;++i)
        {
            try {
                char[] char_array=new char[8];
                fr.read(char_array, 0,8);
                double_array[i]=Double.parseDouble(String.valueOf(char_array));
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        double[] auxdata=new double[integer_array[1]];
        for (int i=0;i<integer_array[1];++i)
        {   
            try {
                char[] char_array=new char[8];
                fr.read(char_array, 0,8);
                auxdata[i]=Double.parseDouble(String.valueOf(char_array));
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        
        for (int i=0;i<integer_array[2];++i)
        {   
            char[] char_array=new char[1];
            try {
                fr.read(char_array, 0,8);
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
            comment=String.valueOf(char_array);
        }
        
        ArrayList<Float> data =new ArrayList<Float>(integer_array[3]*integer_array[4]);
        for (int i =0;i<integer_array[3]*integer_array[4];++i)
        {
            char[] char_array=new char[integer_array[6]];
            try {
                fr.read(char_array, 0,integer_array[6]);
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
            data.add(Float.parseFloat(String.valueOf(char_array)));
        }
        
        try 
        {
            fr.close();
        } catch (IOException ex) {
            Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
        }

        return (Float[])data.toArray();
    }
    
    void read_img2(String imgFilename)
    {
        FileInfo fi= new FileInfo();
        fi.directory=this.getClass().getProtectionDomain().getCodeSource().getLocation().getPath()+"qstem_images/";
        fi.fileName=imgFilename;
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
    }
    */
    
    public static float[] read_img(String imgFilename)
    {
        FileInputStream fis = null;
        float[] fdata = null;
        try 
        {
            //  read the file into a byte array
            File file = new File(imgFilename);
            fis = new FileInputStream(file);
            byte [] arr = new byte[(int)file.length()];
            fis.read(arr);
            //  create a byte buffer and wrap the array
            ByteBuffer bb = ByteBuffer.wrap(arr);
            //  if the file uses little endian as apposed to network
            //  (big endian, Java's native) format,
            //  then set the byte order of the ByteBuffer
            bb.order(ByteOrder.LITTLE_ENDIAN);
            //  read your integers using ByteBuffer's getInt().
            //  four bytes converted into an integer!
            int[] integer_array=new int[8];
            // String[] integer_names= new String[]{"header_size","param_size", "comment_size", 
            //                                "Nx","Ny","is_complex", "data_size", "version"};
            double[] double_array=new double[3];
            // String[] double_names= new String[]{"t","dx","dy","aux_data"};

            String comment;
            for (int i=0;i<8;++i)
            {
                integer_array[i]=bb.getInt();
            }
            
            for (int i=0;i<3;++i)
            {
                double_array[i]=bb.getDouble();
            }
            double[] auxdata=new double[integer_array[1]];
            for (int i=0;i<integer_array[1];++i)
            {
                auxdata[i]=bb.getDouble();
            }
            
            for (int i=0;i<integer_array[2];++i)
            {   
                byte icomment=bb.get();
            }
            
            fdata= new float[integer_array[3]*integer_array[4]];
            
            for (int i =0;i<integer_array[3]*integer_array[4];++i)
            {
                if (integer_array[6]==8)  
                    fdata[i]=(float)bb.getDouble();
                else
                    fdata[i]=bb.getFloat();
            }
        
            

        } catch (FileNotFoundException ex) {
            Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                fis.close();
            } catch (IOException ex) {
                Logger.getLogger(Read_img.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        if (fdata==null)
        {
            System.out.println("something wrong with reading .img");
        }
       

        
        return fdata;
        
        
    }
    
}


