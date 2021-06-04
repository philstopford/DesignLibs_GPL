using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double dt1(double p, double q, double df)

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
            double[][] coef = new double[5][];
            coef[0] = new[] {1.0e0, 1.0e0, 0.0e0, 0.0e0};
            coef[1] = new[] {0.0e0, 3.0e0, 16.0e0, 5.0e0};
            coef[2] = new[] {0.0e0, 0.0e0, -15.0e0, 17.0e0};
            coef[3] = new [] {19.0e0,3.0e0,0.0e0,-945.0e0};
            coef[4] = new[] {-1920.0e0, 1482.0e0, 776.0e0, 79.0e0};

            double[] denom =  {4.0e0,96.0e0,384.0e0,92160.0e0};
            double denpow;
            double dt1;
            int i;
            int[] ideg =  {
                2,3,4,5
            }
            ;
            double sum;
            double term;
            double x;
            double xp;
            double xx;

            x = Math.Abs(dinvnr(p, q));
            xx = x * x;
            sum = x;
            denpow = 1.0e0;
            for (i = 0; i < 4; i++)
            {
                term = eval_pol(coef[i], ideg[i], xx) * x;
                denpow = denpow * df;
                sum = sum + (term / (denpow * denom[i]));
            }

            if (0.5e0 <= p)
            {
                xp = sum;
            }
            else
            {
                xp = -sum;
            }

            dt1 = xp;

            return dt1;
        }
    }
}