using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static double log_normal_cdf(double x, double mu, double sigma)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOG_NORMAL_CDF evaluates the Log Normal CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the PDF.
            //    0.0 < X.
            //
            //    Input, double MU, SIGMA, the parameters of the PDF.
            //    0.0 < SIGMA.
            //
            //    Output, double LOG_NORMAL_CDF, the value of the CDF.
            //
        {
            double cdf;
            double logx;

            if (x <= 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                logx = Math.Log(x);

                cdf = normal_cdf(logx, mu, sigma);
            }

            return cdf;
        }

        public static double log_normal_cdf_inv(double cdf, double mu, double sigma)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOG_NORMAL_CDF_INV inverts the Log Normal CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double CDF, the value of the CDF.
            //    0.0 <= CDF <= 1.0.
            //
            //    Input, double MU, SIGMA, the parameters of the PDF.
            //    0.0 < SIGMA.
            //
            //    Input, double LOG_NORMAL_CDF_INV, the corresponding argument.
            //
        {
            double logx;
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            logx = normal_cdf_inv(cdf, mu, sigma);

            x = Math.Exp(logx);

            return x;
        }

        public static void log_normal_cdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = Log NormalDistribution [ mu, sigma ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the shape parameter of the distribution.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 12;

            double[] fx_vec =
            {
                0.2275013194817921E-01,
                0.2697049307349095E+00,
                0.5781741008028732E+00,
                0.7801170895122241E+00,
                0.4390310097476894E+00,
                0.4592655190218048E+00,
                0.4694258497695908E+00,
                0.4755320473858733E+00,
                0.3261051056816658E+00,
                0.1708799040927608E+00,
                0.7343256357952060E-01,
                0.2554673736161761E-01
            };

            double[] mu_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
            };

            double[] sigma_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01
            };

            double[] x_vec =
            {
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}