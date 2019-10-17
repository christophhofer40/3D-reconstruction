package christoph;

public class Quaternion 
{
    public float x0, x1, x2, x3; 


    // create a new object with the given components
    public Quaternion(float x0, float x1, float x2, float x3) {
        this.x0 = x0;
        this.x1 = x1;
        this.x2 = x2;
        this.x3 = x3;
    }

    //create new 3D rotation
    public Quaternion(float phi, float[] vec)
    {
        angleToQuaternion(phi,vec);
    }

    public float getAngle()
    {

        return (float)Math.acos(this.x0)*2;
    }

    public float[] getAxisVector()
    {
        float[] v = {(this.x1/norm()),(this.x2/norm()),(this.x3/norm())};
        return v;
    }

    public void angleToQuaternion(float phi, float[] vec)
    {
        if(vec[0]==0 && vec[1]==0 && vec[2]==0)
        {
            this.x0 = 1;
            this.x1 = 0;
            this.x2 = 0;
            this.x3 = 0;   
        }
        else 
        {
            final float cp = (float)Math.cos(0.5*phi);
            final float sp = (float)Math.sin(0.5*phi);
            final float inorm = (1.0f/(float)(Math.sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])));
            this.x0 = cp;
            this.x1 = inorm * sp * vec[0];
            this.x2 = inorm * sp * vec[1];
            this.x3 = inorm * sp * vec[2];
        }
    }


    // return a string representation of the invoking object
    @Override
    public String toString() {
        return x0 + " \t" + x1 + " \t" + x2 + " \t" + x3;
    }



    // return the quaternion norm
    public float norm() {
        return (float)Math.sqrt(x0*x0 + x1*x1 + x2*x2 + x3*x3);
    }

    public void normalize()
    {
        x0=x0/norm();
        x1=x1/norm();
        x2=x2/norm();
        x3=x3/norm();
    }

    // return the quaternion conjugate
    public Quaternion conjugate() {
        return new Quaternion(x0, -x1, -x2, -x3);
    }

    // return a new Quaternion whose value is (this + b)
    public Quaternion plus(Quaternion b) {
        //Quaternion a = this;
        return new Quaternion(x0+b.x0, x1+b.x1, x2+b.x2, x3+b.x3);
    }

    // return a new Quaternion whose value is (this - b)
    public Quaternion minus(Quaternion b) {
        //Quaternion a = this;
        return new Quaternion(x0-b.x0, x1-b.x1, x2-b.x2, x3-b.x3);
    }



    // return a new Quaternion whose value is (this * b)
    public Quaternion times(Quaternion b) {
        Quaternion a = this;
        final float y0 = a.x0*b.x0 - a.x1*b.x1 - a.x2*b.x2 - a.x3*b.x3;
        final float y1 = a.x0*b.x1 + a.x1*b.x0 + a.x2*b.x3 - a.x3*b.x2;
        final float y2 = a.x0*b.x2 - a.x1*b.x3 + a.x2*b.x0 + a.x3*b.x1;
        final float y3 = a.x0*b.x3 + a.x1*b.x2 - a.x2*b.x1 + a.x3*b.x0;
        return new Quaternion(y0, y1, y2, y3);
    }

    // return a new Quaternion whose value is (this * b.conjugate())
    private Quaternion times_conjugate(Quaternion b) {
        Quaternion a = this;
        final float y0 =  a.x0*b.x0 + a.x1*b.x1 + a.x2*b.x2 + a.x3*b.x3;
        final float y1 = -a.x0*b.x1 + a.x1*b.x0 - a.x2*b.x3 + a.x3*b.x2;
        final float y2 = -a.x0*b.x2 + a.x1*b.x3 + a.x2*b.x0 - a.x3*b.x1;
        final float y3 = -a.x0*b.x3 - a.x1*b.x2 + a.x2*b.x1 + a.x3*b.x0;
        return new Quaternion(y0, y1, y2, y3);
    }


    // return a new Quaternion whose value is the inverse of this
    public Quaternion inverse() {
        final float id = 1.0f/(x0*x0 + x1*x1 + x2*x2 + x3*x3);
        return new Quaternion(x0*id, -x1*id, -x2*id, -x3*id);
    }


    // return a / b
    // we use the definition a * b^-1 (as opposed to b^-1 a)
    public Quaternion divides(Quaternion b) {
        Quaternion a = this;
        return a.times(b.inverse());
    }


    //rotates a given 3D point around this normalized Quaternion
    public float[] rotate( float[] pos)
    {
            Quaternion loc = new Quaternion(0.0f, pos[0], pos[1], pos[2]);
            loc = this.times( loc.times_conjugate(this));
            final float[] pos2 = {(float)loc.x1,(float)loc.x2,(float)loc.x3};
            return pos2;
    }

    //rotates a given 3D vector Quaternion around this normalized Quaternion
    public Quaternion rotate( Quaternion pos)
    {
            return this.times( pos.times_conjugate(this) );
    }


}
