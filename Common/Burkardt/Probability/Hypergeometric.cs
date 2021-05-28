using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Hypergeometric
    {
        public static double hypergeometric_cdf(int x, int n, int m, int l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the CDF.
        //
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Output, double HYPERGEOMETRIC_CDF, the value of the CDF.
        //
        {
            int x2;

            double c1_log = typeMethods.i4_choose_log(l - m, n);
            double c2_log = typeMethods.i4_choose_log(l, n);

            double pdf = Math.Exp(c1_log - c2_log);
            double cdf = pdf;

            for (x2 = 0; x2 <= x - 1; x2++)
            {
                pdf = pdf * (double) ((m - x2) * (n - x2))
                      / (double) ((x2 + 1) * (l - m - n + x2 + 1));

                cdf = cdf + pdf;
            }

            return cdf;
        }

        public static void hypergeometric_cdf_values(ref int n_data, ref int sam, ref int suc, ref int pop,
                ref int n, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
        //
        //  Discussion:
        //
        //    CDF(X)(A,B) is the probability of at most X successes in A trials,
        //    given that the probability of success on a single trial is B.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = HypergeometricDistribution [ sam, suc, pop ]
        //      CDF [ dist, n ]
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
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition, CRC Press, 1996, pages 651-652.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &SAM, int &SUC, int &POP, the sample size,
        //    success size, and population parameters of the function.
        //
        //    Output, int &N, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.6001858177500578E-01,
                0.2615284665839845E+00,
                0.6695237889132748E+00,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.5332595856827856E+00,
                0.1819495964117640E+00,
                0.4448047017527730E-01,
                0.9999991751316731E+00,
                0.9926860896560750E+00,
                0.8410799901444538E+00,
                0.3459800113391901E+00,
                0.0000000000000000E+00,
                0.2088888139634505E-02,
                0.3876752992448843E+00,
                0.9135215248834896E+00
            };

            int[] n_vec =
            {
                7, 8, 9, 10,
                6, 6, 6, 6,
                6, 6, 6, 6,
                0, 0, 0, 0
            };

            int[] pop_vec =
            {
                100, 100, 100, 100,
                100, 100, 100, 100,
                100, 100, 100, 100,
                90, 200, 1000, 10000
            };

            int[] sam_vec =
            {
                10, 10, 10, 10,
                6, 7, 8, 9,
                10, 10, 10, 10,
                10, 10, 10, 10
            };

            int[] suc_vec =
            {
                90, 90, 90, 90,
                90, 90, 90, 90,
                10, 30, 50, 70,
                90, 90, 90, 90
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                sam = 0;
                suc = 0;
                pop = 0;
                n = 0;
                fx = 0.0;
            }
            else
            {
                sam = sam_vec[n_data - 1];
                suc = suc_vec[n_data - 1];
                pop = pop_vec[n_data - 1];
                n = n_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool hypergeometric_check(int n, int m, int l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_CHECK checks the parameters of the Hypergeometric CDF.
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
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Output, bool HYPERGEOMETRIC_CHECK, is true if the parameters are legal.
        //
        {
            if (n < 0 || l < n)
            {
                Console.WriteLine(" ");
                Console.WriteLine("HYPERGEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  Input N is out of range.");
                return false;
            }

            if (m < 0 || l < m)
            {
                Console.WriteLine(" ");
                Console.WriteLine("HYPERGEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  Input M is out of range.");
                return false;
            }

            if (l < 0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("HYPERGEOMETRIC_CHECK - Warning!");
                Console.WriteLine("  Input L is out of range.");
                return false;
            }

            return true;
        }

        public static double hypergeometric_mean(int n, int m, int l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_MEAN returns the mean of the Hypergeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Output, double HYPERGEOMETRIC_MEAN, the mean of the PDF.
        //
        {
            double mean = (double) (n * m) / (double) (l);

            return mean;
        }

        public static double hypergeometric_pdf(int x, int n, int m, int l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.
        //
        //  Discussion:
        //
        //    PDF(N,M,L;X) = C(M,X) * C(L-M,N-X) / C(L,N).
        //
        //    PDF(N,M,L;X) is the probability of drawing X white balls in a
        //    single random sample of size N from a population containing
        //    M white balls and a total of L balls.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the desired number of white balls.
        //    0 <= X <= N, usually, although any value of X can be given.
        //
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Output, double HYPERGEOMETRIC_PDF, the probability of exactly K white balls.
        //
        {
            double pdf;
            //
            //  Special cases.
            //
            if (x < 0)
            {
                pdf = 1.0;
            }
            else if (n < x)
            {
                pdf = 0.0;
            }
            else if (m < x)
            {
                pdf = 0.0;
            }
            else if (l < x)
            {
                pdf = 0.0;
            }
            else if (n == 0)
            {
                if (x == 0)
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
                double c1 = typeMethods.i4_choose_log(m, x);
                double c2 = typeMethods.i4_choose_log(l - m, n - x);
                double c3 = typeMethods.i4_choose_log(l, n);

                double pdf_log = c1 + c2 - c3;

                pdf = Math.Exp(pdf_log);

            }

            return pdf;
        }

        public static int hypergeometric_sample(int n, int m, int l, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_SAMPLE samples the Hypergeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 165.
        //
        //  Parameters:
        //
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int HYPERGEOMETRIC_SAMPLE, a sample of the PDF.
        //
        {
            double c1_log = typeMethods.i4_choose_log(l - m, n);
            double c2_log = typeMethods.i4_choose_log(l, n);

            double a = Math.Exp(c1_log - c2_log);
            double b = a;

            double u = UniformRNG.r8_uniform_01(ref seed);

            int x = 0;

            while (a < u)
            {
                b = b * (double) ((m - x) * (n - x))
                    / (double) ((x + 1) * (l - m - n + x + 1));

                a = a + b;

                x = x + 1;

            }

            return x;
        }

        public static double hypergeometric_variance(int n, int m, int l)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERGEOMETRIC_VARIANCE returns the variance of the Hypergeometric PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of balls selected.
        //    0 <= N <= L.
        //
        //    Input, int M, the number of white balls in the population.
        //    0 <= M <= L.
        //
        //    Input, int L, the number of balls to select from.
        //    0 <= L.
        //
        //    Output, double HYPERGEOMETRIC_VARIANCE, the variance of the PDF.
        //
        {
            double variance = (double) (n * m * (l - m) * (l - n))
                              / (double) (l * l * (l - 1));

            return variance;
        }
    }
}