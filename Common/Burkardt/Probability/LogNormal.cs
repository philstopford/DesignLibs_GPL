using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class LogNormal
    {
        public static double log_normal_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CDF evaluates the Lognormal CDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                double logx = Math.Log(x);

                cdf = Normal.normal_cdf(logx, a, b);
            }

            return cdf;
        }

        public static double log_normal_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CDF_INV inverts the Lognormal CDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Input, double LOG_NORMAL_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double logx = Normal.normal_cdf_inv(cdf, a, b);

            double x = Math.Exp(logx);

            return x;
        }

        public static void log_normal_cdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double x, ref double fx )
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
        //      dist = LogNormalDistribution [ mu, sigma ]
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

        public static bool log_normal_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CHECK checks the parameters of the Lognormal PDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, bool LOG_NORMAL_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LOG_NORMAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            return true;
        }

        public static double log_normal_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_MEAN returns the mean of the Lognormal PDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = Math.Exp(a + 0.5 * b * b);

            return mean;
        }

        public static double log_normal_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_PDF evaluates the Lognormal PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X)
        //      = exp ( - 0.5 * ( ( log ( X ) - A ) / B )^2 )
        //        / ( B * X * sqrt ( 2 * PI ) )
        //
        //    The Lognormal PDF is also known as the Cobb-Douglas PDF,
        //    and as the Antilog_normal PDF.
        //
        //    The Lognormal PDF describes a variable X whose logarithm
        //    is normally distributed.
        //
        //    The special case A = 0, B = 1 is known as Gilbrat's PDF.
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
        //    0.0 < X
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= 0.0)
            {
                pdf = 0.0;
            }
            else
            {
                double y = (Math.Log(x) - a) / b;
                pdf = Math.Exp(-0.5 * y * y) / (b * x * Math.Sqrt(2.0 * r8_pi));
            }

            return pdf;
        }

        public static double log_normal_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_SAMPLE samples the Lognormal PDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double LOG_NORMAL_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = log_normal_cdf_inv(cdf, a, b);

            return x;
        }

        public static double log_normal_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_VARIANCE returns the variance of the Lognormal PDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            double variance = Math.Exp(2.0 * a + b * b) * (Math.Exp(b * b) - 1.0);

            return variance;
        }
    }
}