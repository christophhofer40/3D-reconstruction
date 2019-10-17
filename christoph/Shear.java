package christoph;

public class Shear
    {
        public float[][] map = {{1.0f,0.0f},{0.0f,1.0f}};

        void apply(float[] pos)
        {
            final float x =(float) (pos[0]*map[0][0] + pos[1]*map[0][1]);
            final float y =(float) (pos[0]*map[1][0] + pos[1]*map[1][1]);
            pos[0] = x;
            pos[1] = y;
        }
        
                

        void apply_inv(float[] pos) //assumes det == 1
        {
            final float x = (float) ( pos[0]*map[1][1] - pos[1]*map[0][1]);
            final float y = (float) (-pos[0]*map[1][0] + pos[1]*map[0][0]);
            pos[0] = x;
            pos[1] = y;
        }

        public void normalize()
        {
            final double idet = Math.sqrt(1.0/(map[0][0]*map[1][1]-map[1][0]*map[0][1]));
            map[0][0] *= idet;
            map[0][1] *= idet;
            map[1][0] *= idet;
            map[1][1] *= idet;
        }

        void reset()
        {
            map[0][0] = 1.0f;
            map[0][1] = 0.0f;
            map[1][0] = 0.0f;
            map[1][1] = 1.0f;
        }

        public String getStr()
        {
            final double det = (map[0][0]*map[1][1]-map[1][0]*map[0][1]);
            return "{{"+map[0][0]+","+map[0][1]+"},{"+map[1][0]+","+map[1][1]+"}} det: " + det;
        }
    }
