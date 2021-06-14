using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Logistic
    {
        public static double logistic_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_CDF evaluates the Logistic CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double LOGISTIC_CDF, the value of the CDF.
        //
        {
            double cdf = 1.0 / (1.0 + Math.Exp((a - x) / b));

            return cdf;
        }

        public static double logistic_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_CDF_INV inverts the Logistic CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, double LOGISTIC_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("LOGISTIC_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double x = a - b * Math.Log((1.0 - cdf) / cdf);

            return x;
        }

        public static void logistic_cdf_values(ref int n_data, ref double mu, ref double beta, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = LogisticDistribution [ mu, beta ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 August 2004
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
        //    Output, double &BETA, the shape parameter of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 12;

            double[] beta_vec =
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

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.8807970779778824E+00,
                0.9820137900379084E+00,
                0.9975273768433652E+00,
                0.6224593312018546E+00,
                0.5825702064623147E+00,
                0.5621765008857981E+00,
                0.5498339973124779E+00,
                0.6224593312018546E+00,
                0.5000000000000000E+00,
                0.3775406687981454E+00,
                0.2689414213699951E+00
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
                beta = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                mu = mu_vec[n_data - 1];
                beta = beta_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool logistic_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_CHECK checks the parameters of the Logistic CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, bool LOGISTIC_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("LOGISTIC_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            return true;
        }

        public static double logistic_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_MEAN returns the mean of the Logistic PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, double LOGISTIC_MEAN, the mean of the PDF.
        //
        {
            double mean = a;

            return mean;
        }

        public static double logistic_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_PDF evaluates the Logistic PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = exp ( ( A - X ) / B ) /
        //      ( B * ( 1 + exp ( ( A - X ) / B ) )^2 )
        //
        //    The Logistic PDF is also known as the Sech-Squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double LOGISTIC_PDF, the value of the PDF.
        //
        {
            double temp = Math.Exp((a - x) / b);

            double pdf = temp / (b * (1.0 + temp) * (1.0 + temp));

            return pdf;
        }

        public static double logistic_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_SAMPLE samples the Logistic PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, double LOGISTIC_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = logistic_cdf_inv(cdf, a, b);

            return x;
        }

        public static double logistic_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOGISTIC_VARIANCE returns the variance of the Logistic PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, double LOGISTIC_VARIANCE, the variance of the PDF.
        //
        {
            

            double variance = Math.PI * Math.PI * b * b / 3.0;

            return variance;
        }
    }
}