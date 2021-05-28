using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Geometric
    {
        public static double geometric_cdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_CDF evaluates the Geometric CDF.
        //
        //  Discussion:
        //
        //    CDF(X,P) is the probability that there will be at least one
        //    successful trial in the first X Bernoulli trials, given that
        //    the probability of success in a single trial is P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the maximum number of trials.
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double GEOMETRIC_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= 0)
            {
                cdf = 0.0;
            }
            else if (a == 0.0)
            {
                cdf = 0.0;
            }
            else if (a == 1.0)
            {
                cdf = 1.0;
            }
            else
            {
                cdf = 1.0 - Math.Pow((1.0 - a), x);
            }

            return cdf;
        }

        public static int geometric_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_CDF_INV inverts the Geometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, int GEOMETRIC_CDF_INV, the corresponding value of X.
        //
        {
            int x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GEOMETRIC_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            if (a == 1.0)
            {
                x = 1;
            }
            else if (a == 0.0)
            {
                x = typeMethods.i4_huge();
            }
            else
            {
                x = 1 + (int) (Math.Log(1.0 - cdf) / Math.Log(1.0 - a));
            }

            return x;
        }

        public static void geometric_cdf_values(ref int n_data, ref int x, ref double p, ref double cdf )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
        //
        //  Discussion:
        //
        //    The geometric or Pascal probability density function gives the
        //    probability that the first success will happen on the X-th Bernoulli
        //    trial, given that the probability of a success on a single trial is P.
        //
        //    The value of CDF ( X, P ) is the probability that the first success
        //    will happen on or before the X-th trial.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = GeometricDistribution [ p ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //    Daniel Zwillinger and Stephen Kokoska,
        //    CRC Standard Probability and Statistics Tables and Formulae,
        //    Chapman and Hall / CRC Press, 2000.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &X, the number of trials.
        //
        //    Output, double &P, the probability of success
        //    on one trial.
        //
        //    Output, double &CDF, the cumulative density function.
        //
        {
            int N_MAX = 14;

            double[] cdf_vec =
            {
                0.1900000000000000E+00,
                0.2710000000000000E+00,
                0.3439000000000000E+00,
                0.6861894039100000E+00,
                0.3600000000000000E+00,
                0.4880000000000000E+00,
                0.5904000000000000E+00,
                0.9141006540800000E+00,
                0.7599000000000000E+00,
                0.8704000000000000E+00,
                0.9375000000000000E+00,
                0.9843750000000000E+00,
                0.9995117187500000E+00,
                0.9999000000000000E+00
            };

            double[] p_vec =
            {
                0.1E+00,
                0.1E+00,
                0.1E+00,
                0.1E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.3E+00,
                0.4E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.9E+00
            };

            int[] x_vec =
            {
                1, 2, 3, 10, 1,
                2, 3, 10, 3, 3,
                3, 5, 10, 3
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0;
                p = 0.0;
                cdf = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                p = p_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
            }
        }

        public static bool geometric_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_CHECK checks the parameter of the Geometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, bool GEOMETRIC_CHECK, is true if the parameters are legal.
        //
        {
            if (a < 0.0 || 1.0 < a)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  A < 0 or 1 < A.");
                return false;
            }

            return true;
        }

        public static double geometric_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_MEAN returns the mean of the Geometric PDF.
        //
        //  Discussion:
        //
        //    MEAN is the expected value of the number of trials required
        //    to obtain a single success.
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
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = 1.0 / a;

            return mean;
        }

        public static double geometric_pdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_PDF evaluates the Geometric PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = A * ( 1 - A )^(X-1)
        //
        //    PDF(A;X) is the probability that exactly X Bernoulli trials, each
        //    with probability of success A, will be required to achieve
        //    a single success.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the number of trials.
        //    0 < X
        //
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;
            //
            //  Special cases.
            //
            if (x < 1)
            {
                pdf = 0.0;
            }
            else if (a == 0.0)
            {
                pdf = 0.0;
            }
            else if (a == 1.0)
            {
                if (x == 1)
                {
                    pdf = 1.0;
                }
                else
                {
                    pdf = 0.0;
                }
            }
            else
            {
                pdf = a * Math.Pow((1.0 - a), (x - 1));

            }

            return pdf;
        }

        public static int geometric_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_SAMPLE samples the Geometric PDF.
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
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int GEOMETRIC_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            int x = geometric_cdf_inv(cdf, a);

            return x;
        }

        public static double geometric_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_VARIANCE returns the variance of the Geometric PDF.
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
        //    Input, double A, the probability of success on one trial.
        //    0.0 <= A <= 1.0.
        //
        //    Output, double GEOMETRIC_VARIANCE, the variance of the PDF.
        //
        {
            double variance = (1.0 - a) / (a * a);

            return variance;
        }
    }
}