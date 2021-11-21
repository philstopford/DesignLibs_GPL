using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double dt1 ( double p, double q, double df )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DT1 computes an approximate inverse of the cumulative T distribution.
        //
        //  Discussion:
        //
        //    Returns the inverse of the T distribution function, i.e.,
        //    the integral from 0 to INVT of the T density is P. This is an
        //    initial approximation.
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
        //    Input, double *P, *Q, the value whose inverse from the
        //    T distribution CDF is desired, and the value (1-P).
        //
        //    Input, double *DF, the number of degrees of freedom of the
        //    T distribution.
        //
        //    Output, double DT1, the approximate value of X for which
        //    the T density CDF with DF degrees of freedom has value P.
        //
    {
        double[][] coef = new double[4][];

        coef[0] = new [] {1.0e0, 1.0e0, 0.0e0 , 0.0e0, 0.0e0};
        coef[1] = new [] {3.0e0, 16.0e0, 5.0e0, 0.0e0, 0.0e0};
        coef[2] = new[] {-15.0e0, 17.0e0, 19.0e0, 3.0e0, 0.0e0};
        coef[3] = new[] {-945.0e0,-1920.0e0,1482.0e0,776.0e0,79.0e0};
            
        double[] denom = {
            4.0e0,96.0e0,384.0e0,92160.0e0
        };
        int i;
        int[] ideg = {
            2,3,4,5
        };

        double x = Math.Abs ( dinvnr ( p, q ) );
        double xx = x * x;
        double sum = x;
        double denpow = 1.0e0;
        for ( i = 0; i < 4; i++ )
        {
            double term = eval_pol ( coef[i], ideg[i], xx ) * x;
            denpow *= df;
            sum += term / ( denpow * denom[i] );
        }

        double xp = p switch
        {
            >= 0.5e0 => sum,
            _ => -sum
        };
        double dt1 = xp;

        return dt1;
    }
}