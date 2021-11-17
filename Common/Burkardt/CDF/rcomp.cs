using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double rcomp ( double a, double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RCOMP evaluates exp(-X) * X**A / Gamma(A).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Input:
        //
        //    double *A, *X, arguments of the quantity to be computed.
        //
        //  Output:
        //
        //    double RCOMP, the value of exp(-X) * X^A / Gamma(A).
        //
        //  Local:
        //
        //    double RT2PIN = 1/SQRT(2*PI)
        //
    {
        double rt2pin = .398942280401433e0;
        double rcomp,t,t1,u;
        rcomp = 0.0e0;
        switch (a)
        {
            case >= 20.0e0:
                goto S20;
        }

        t = a*Math.Log(x)-x;
        switch (a)
        {
            case >= 1.0e0:
                goto S10;
        }

        rcomp = a*Math.Exp(t)*(1.0e0+gam1(a));
        return rcomp;
        S10:
        rcomp = Math.Exp(t)/ gamma_x(a);
        return rcomp;
        S20:
        u = x/ a;
        switch (u)
        {
            case 0.0e0:
                return rcomp;
        }

        t = Math.Pow(1.0e0/ a,2.0);
        t1 = (((0.75e0*t-1.0e0)*t+3.5e0)*t-105.0e0)/(a*1260.0e0);
        t1 -= a*rlog(u);
        rcomp = rt2pin*Math.Sqrt(a)*Math.Exp(t1);
        return rcomp;
    }
}