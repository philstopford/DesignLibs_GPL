using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static void beta_noncentral_cdf_values(ref int n_data, ref double a, ref double b,
                                                ref double lambda, ref double x, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
        //
        //  Discussion:
        //
        //    The values presented here are taken from the reference, where they
        //    were given to a limited number of decimal places.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    R Chattamvelli, R Shanmugam,
        //    Algorithm AS 310:
        //    Computing the Non-central Beta Distribution Function,
        //    Applied Statistics,
        //    Volume 46, Number 1, 1997, pages 146-156.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0
        //    before the first call.  On each call, the routine increments N_DATA by 1,
        //    and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double *A, *B, the shape parameters.
        //
        //    Output, double *LAMBDA, the noncentrality parameter.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
        {
            int N_MAX = 25;

            double[] a_vec =  {
                5.0,
                5.0,
                5.0,
                10.0,
                10.0,
                10.0,
                20.0,
                20.0,
                20.0,
                10.0,
                10.0,
                15.0,
                20.0,
                20.0,
                20.0,
                30.0,
                30.0,
                10.0,
                10.0,
                10.0,
                15.0,
                10.0,
                12.0,
                30.0,
                35.0
            };
            
            double[] b_vec =  {
                5.0,
                5.0,
                5.0,
                10.0,
                10.0,
                10.0,
                20.0,
                20.0,
                20.0,
                20.0,
                10.0,
                5.0,
                10.0,
                30.0,
                50.0,
                20.0,
                40.0,
                5.0,
                10.0,
                30.0,
                20.0,
                5.0,
                17.0,
                30.0,
                30.0
            };
            
            double[] fx_vec =  {
                0.4563021,
                0.1041337,
                0.6022353,
                0.9187770,
                0.6008106,
                0.0902850,
                0.9998655,
                0.9925997,
                0.9641112,
                0.9376626573,
                0.7306817858,
                0.1604256918,
                0.1867485313,
                0.6559386874,
                0.9796881486,
                0.1162386423,
                0.9930430054,
                0.0506899273,
                0.1030959706,
                0.9978417832,
                0.2555552369,
                0.0668307064,
                0.0113601067,
                0.7813366615,
                0.8867126477
            };
            
            double[] lambda_vec =  {
                54.0,
                140.0,
                170.0,
                54.0,
                140.0,
                250.0,
                54.0,
                140.0,
                250.0,
                150.0,
                120.0,
                80.0,
                110.0,
                65.0,
                130.0,
                80.0,
                130.0,
                20.0,
                54.0,
                80.0,
                120.0,
                55.0,
                64.0,
                140.0,
                20.0
            };
            
            double[] x_vec =  {
                0.8640,
                0.9000,
                0.9560,
                0.8686,
                0.9000,
                0.9000,
                0.8787,
                0.9000,
                0.9220,
                0.868,
                0.900,
                0.880,
                0.850,
                0.660,
                0.720,
                0.720,
                0.800,
                0.644,
                0.700,
                0.780,
                0.760,
                0.795,
                0.560,
                0.800,
                0.670
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                b = 0.0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }


        public static double betanc(double x, double a, double b, double lambda, ref int ifault)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BETANC computes the tail of the noncentral Beta distribution.
        //
        //  Discussion:
        //
        //    This routine returns the cumulative probability of X for the non-central
        //    Beta distribution with parameters A, B and non-centrality LAMBDA.
        //
        //    Note that if LAMBDA = 0, the standard Beta distribution is defined.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Russell Lenth.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Russell Lenth,
        //    Algorithm AS 226:
        //    Computing Noncentral Beta Probabilities,
        //    Applied Statistics,
        //    Volume 36, Number 2, 1987, pages 241-244.
        //
        //    H Frick,
        //    Algorithm AS R84:
        //    A Remark on Algorithm AS 226:
        //    Computing Noncentral Beta Probabilities,
        //    Applied Statistics,
        //    Volume 39, Number 2, 1990, pages 311-312.
        //
        //  Parameters:
        //
        //    Input, double X, the value defining the cumulative
        //    probability lower tail.  Normally, 0 <= X <= 1, but any value
        //    is allowed.
        //
        //    Input, double A, B, the parameters of the distribution.
        //    0 < A, 0 < B.
        //
        //    Input, double LAMBDA, the noncentrality parameter
        //    of the distribution.  0 <= LAMBDA.  The program can produce reasonably
        //    accurate results for values of LAMBDA up to about 100.
        //
        //    Output, int *IFAULT, error flag.
        //    0, no error occurred.
        //    nonzero, an error occurred.
        //
        //    Output, double BETANC, the cumulative probability
        //    of X.
        //
        {
            double ax;
            double beta;
            double c;
            double errbd;
            double errmax = 1.0E-07;
            double gx;
            int itrmax = 150;
            double q;
            double sumq;
            double temp;
            double value;
            double xj;

            ifault = 0;

            if (lambda < 0.0 ||
                a <= 0.0 ||
                b <= 0.0)
            {
                ifault = 2;
                value = -1.0;
                return value;
            }

            if (x <= 0.0)
            {
                value = 0.0;
                return value;
            }

            if (1.0 <= x)
            {
                value = 1.0;
                return value;
            }

            c = 0.5 * lambda;
            //
            //  Initialize the series.
            //
            beta = Helpers.LogGamma(a)
                   + Helpers.LogGamma(b)
                   - Helpers.LogGamma(a + b);

            temp = betain(x, a, b, beta, ref ifault);

            gx = Math.Exp(a * Math.Log(x) + b * Math.Log(1.0 - x) - beta - Math.Log(a));

            q = Math.Exp(-c);

            xj = 0.0;
            ax = q * temp;
            sumq = 1.0 - q;
            value = ax;
            //
            //  Recur over subsequent terms until convergence is achieved.
            //
            ifault = 1;

            for (;;)
            {
                xj = xj + 1.0;
                temp = temp - gx;
                gx = x * (a + b + xj - 1.0) * gx / (a + xj);
                q = q * c / xj;
                sumq = sumq - q;
                ax = temp * q;
                value = value + ax;
                //
                //  Check for convergence and act accordingly.
                //
                errbd = Math.Abs((temp - gx) * sumq);

                if (errbd <= errmax)
                {
                    ifault = 0;
                    break;
                }

                if (itrmax < (int) xj)
                {
                    break;
                }
            }

            return value;
        }
    }
}