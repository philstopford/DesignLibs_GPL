using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double dlanor ( double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DLANOR evaluates the logarithm of the asymptotic Normal CDF.
        //
        //  Discussion:
        //
        //    This routine computes the logarithm of the cumulative normal distribution
        //    from abs ( x ) to infinity for  5 <= abs ( X ).
        //
        //    The relative error at X = 5 is about 0.5E-5.
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
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions
        //    1966, Formula 26.2.12.
        //
        //  Parameters:
        //
        //    Input, double *X, the value at which the Normal CDF is to be
        //    evaluated.  It is assumed that 5 <= abs ( X ).
        //
        //    Output, double DLANOR, the logarithm of the asymptotic
        //    Normal CDF.
        //
    {
        double[] coef = {
            -1.0e0,3.0e0,-15.0e0,105.0e0,-945.0e0,10395.0e0,-135135.0e0,2027025.0e0,
            -34459425.0e0,654729075.0e0,-13749310575e0,316234143225.0e0
        };
        const double dlsqpi = 0.91893853320467274177e0;
        const int K1 = 12;

        double xx = Math.Abs ( x );
        switch (xx)
        {
            case < 5.0e0:
                throw new Exception( "DLANOR: Argument too small.");
        }
        double approx = -dlsqpi-0.5e0*xx*xx-Math.Log(xx);
        double xx2 = xx*xx;
        double T2 = 1.0e0/xx2;
        double correc = eval_pol ( coef, K1, T2 ) / xx2;
        correc = alnrel ( correc );
        double dlanor = approx+correc;

        return dlanor;
    }
}