using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double fpser ( double a, double b, double x, double eps )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FPSER evaluates IX(A,B)(X) for very small B.
            //
            //  Discussion:
            //
            //    This routine is appropriate for use when
            //
            //      B < min ( EPS, EPS * A )
            //
            //    and
            //
            //      X <= 0.5.
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
            //  Parameters:
            //
            //    Input, double *A, *B, parameters of the function.
            //
            //    Input, double *X, the point at which the function is to
            //    be evaluated.
            //
            //    Input, double *EPS, a tolerance.
            //
            //    Output, double FPSER, the value of IX(A,B)(X).
            //
        {
            double an;
            double c;
            double fpser;
            int K1 = 1;
            double s;
            double t;
            double tol;

            fpser = 1.0e0;
            if(a <= 1e-3*eps) goto S10;
            fpser = 0.0e0;
            t = a*Math.Log(x);
            if(t < exparg(K1)) return fpser;
            fpser = Math.Exp(t);
            S10:
            //
            //                NOTE THAT 1/B(A,B) = B
            //
            fpser = b/ a*fpser;
            tol = eps/ a;
            an = a+1.0e0;
            t = x;
            s = t/an;
            S20:
            an = an + 1.0e0;
            t = x*t;
            c = t/an;
            s = s + c;
            if(Math.Abs(c) > tol) goto S20;
            fpser = fpser * (1.0e0+a*s);

            return fpser;
        }
    }
}