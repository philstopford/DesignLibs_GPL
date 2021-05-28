using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Cauchy
    {
        static double cauchy_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_CDF evaluates the Cauchy CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, double CDF, the value of the CDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double y = (x - a) / b;

            double cdf = 0.5 + Math.Atan(y) / r8_pi;

            return cdf;
        }

        static double cauchy_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_CDF_INV inverts the Cauchy CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, double CAUCHY_CDF_INV, the corresponding argument.
        //
        {
            const double r8_pi = 3.14159265358979323;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("CAUCHY_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            }

            double x = a + b * Math.Tan(r8_pi * (cdf - 0.5));

            return x;
        }

        static void cauchy_cdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
                                        ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = CauchyDistribution [ mu, sigma ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2004
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
        //    Output, double &SIGMA, the variance of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 12;

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.8524163823495667E+00,
                0.9220208696226307E+00,
                0.9474315432887466E+00,
                0.6475836176504333E+00,
                0.6024163823495667E+00,
                0.5779791303773693E+00,
                0.5628329581890012E+00,
                0.6475836176504333E+00,
                0.5000000000000000E+00,
                0.3524163823495667E+00,
                0.2500000000000000E+00
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

            return;
        }

        static bool cauchy_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_CHECK checks the parameters of the Cauchy CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, bool CAUCHY_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("CAUCHY_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            return true;
        }

        static double cauchy_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_MEAN returns the mean of the Cauchy PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
            double mean = a;

            return mean;
        }

        static double cauchy_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_PDF evaluates the Cauchy PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = 1 / ( PI * B * ( 1 + ( ( X - A ) / B )^2 ) )
        //
        //    The Cauchy PDF is also known as the Breit-Wigner PDF.  It
        //    has some unusual properties.  In particular, the integrals for the
        //    expected value and higher order moments are "singular", in the
        //    sense that the limiting values do not exist.  A result can be
        //    obtained if the upper and lower limits of integration are set
        //    equal to +T and -T, and the limit as T=>INFINITY is taken, but
        //    this is a very weak and unreliable sort of limit.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, double PDF, the value of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double y = (x - a) / b;

            double pdf = 1.0 / (r8_pi * b * (1.0 + y * y));

            return pdf;
        }

        static double cauchy_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_SAMPLE samples the Cauchy PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, double CAUCHY_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = cauchy_cdf_inv(cdf, a, b);

            return x;
        }

        static double cauchy_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CAUCHY_VARIANCE returns the variance of the Cauchy PDF.
        //
        //  Discussion:
        //
        //    The variance of the Cauchy PDF is not well defined.  This routine
        //    is made available for completeness only, and simply returns
        //    a "very large" number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Output, double VARIANCE, the mean of the PDF.
        //
        {
            const double r8_huge = 1.0E+30;

            double variance = r8_huge;

            return variance;
        }
    }
}