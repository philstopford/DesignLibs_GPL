using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double dbetrm ( double a, double b )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DBETRM computes the Sterling remainder for the complete beta function.
            //
            //  Discussion:
            //
            //    Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
            //    where Lgamma is the log of the (complete) gamma function
            //
            //    Let ZZ be approximation obtained if each log gamma is approximated
            //    by Sterling's formula, i.e.,
            //    Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
            //
            //    The Sterling remainder is Log(Beta(A,B)) - ZZ.
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
            //    Input, double *A, *B, the parameters of the Beta function.
            //
            //    Output, double DBETRM, the Sterling remainder.
            //
        {
            //
            //  Sum from smallest to largest
            //
            double T1 = a+b;
            double dbetrm = -dstrem(T1);
            double T2 = Math.Max(a,b);
            dbetrm = dbetrm + dstrem(T2);
            double T3 = Math.Min(a,b);
            dbetrm = dbetrm + dstrem(T3);

            return dbetrm;
        }
    }
}