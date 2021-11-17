using System.Runtime.InteropServices;

namespace NextAfterNS;

public static class NextAfter
{
    [StructLayout(LayoutKind.Explicit)]
    private struct FloatIntUnion
    {
        [FieldOffset(0)]
        public int i;
        [FieldOffset(0)]
        public float f;
    }

    //  Returns the next float after x in the direction of y.
    public static float nextafter(float x, float y)
    {
        if (float.IsNaN(x) || float.IsNaN(y))
        {
            return x + y;
        }

        if (x == y)
        {
            return y;  // nextafter(0, -0) = -0
        }

        FloatIntUnion u;
        u.i = 0; u.f = x;  // shut up the compiler

        switch (x)
        {
            case 0:
                u.i = 1;
                return y > 0 ? u.f : -u.f;
        }

        if (x > 0 == y > x)
        {
            u.i++;
        }
        else
        {
            u.i--;
        }

        return u.f;
    }        
        
    [StructLayout(LayoutKind.Explicit)]
    private struct DoubleIntUnion
    {
        [FieldOffset(0)]
        public int i;
        [FieldOffset(0)]
        public double f;
    }

    //  Returns the next double after x in the direction of y.
    public static double nextafter(double x, double y)
    {
        if (double.IsNaN(x) || double.IsNaN(y))
        {
            return x + y;
        }

        if (x == y)
        {
            return y;  // nextafter(0, -0) = -0
        }

        DoubleIntUnion u;
        u.i = 0; u.f = x;  // shut up the compiler

        switch (x)
        {
            case 0:
                u.i = 1;
                return y > 0 ? u.f : -u.f;
        }

        if (x > 0 == y > x)
        {
            u.i++;
        }
        else
        {
            u.i--;
        }

        return u.f;
    }        
}