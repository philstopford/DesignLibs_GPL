using System;
using Burkardt.AppliedStatistics;

namespace Burkardt;

public static class BetaNC
{
    public static double beta_noncentral_cdf ( double a, double b, double lambda, double x,
            double error_max )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //  BETA_NONCENTRAL_CDF evaluates the noncentral Beta CDF.
        //
        //  Discussion:
        //
        //    The reference mistakenly phrases the opposite of the correct
        //    stopping criterion, that is, it says:
        //
        //      "stop when PSUM < 1 - ERROR"
        //
        //    but this must be:
        //
        //      "stop when 1 - ERROR < PSUM."
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Harry Posten,
        //    An Effective Algorithm for the Noncentral Beta Distribution Function,
        //    The American Statistician,
        //    Volume 47, Number 2, May 1993, pages 129-131.
        //
        //  Parameters:
        //
        //    Input, double A, B, the shape parameters.
        //
        //    Input, double LAMBDA, the noncentrality parameter.
        //
        //    Input, double X, the argument of the function.
        //
        //    Input, double ERROR_MAX, the error control.
        //
        //    Output, double BETA_NONCENTRAL_CDF, the value of the noncentral Beta CDF.
        //
    {
        int ifault = 0;

        int i = 0;
        double pi = Math.Exp ( - lambda / 2.0 );

        double beta_log = Algorithms.alogam ( a, ref ifault )
                          + Algorithms.alogam ( b, ref ifault )
                          - Algorithms.alogam ( a + b, ref ifault );

        double bi = Algorithms.betain ( x, a, b, beta_log, ref ifault );

        double si = Math.Exp (
            a * Math.Log ( x )
            + b * Math.Log ( 1.0 - x )
            - beta_log
            - Math.Log ( a ) );

        double p_sum = pi;
        double pb_sum = Math.PI * bi;

        while ( p_sum < 1.0 - error_max )
        {
            double pj = pi;
            double bj = bi;
            double sj = si;

            i += 1;
            pi = 0.5 * lambda * pj / i;
            bi = bj - sj;
            si = x * ( a + b + i - 1 ) * sj / ( a + i );

            p_sum += pi;
            pb_sum += Math.PI * bi;
        }

        double value = pb_sum;

        return value;
    }
}