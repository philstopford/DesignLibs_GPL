using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Runs
    {
        public static double runs_mean(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RUNS_MEAN returns the mean of the Runs PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //
        //    Output, double RUNS_MEAN, the mean of the PDF.
        //
        {
            double mean = (double) (m + 2 * m * n + n)
                          / (double) (m + n);

            return mean;
        }

        public static double runs_pdf(int m, int n, int r)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RUNS_PDF evaluates the Runs PDF.
        //
        //  Discussion:
        //
        //    Suppose we have M symbols of one type and N of another, and we consider
        //    the various possible permutations of these symbols.
        //
        //    Let "R" be the number of runs in a given permutation.  By a "run", we
        //    mean a maximal sequence of identical symbols.  Thus, for instance,
        //    the permutation
        //
        //      ABBBAAAAAAAA
        //
        //    has three runs.
        //
        //    The probability that a permutation of M+N symbols, with M of one kind
        //    and N of another, will have exactly R runs is:
        //
        //      PDF(M,N)(R) = 2 * C(M-1,R/2-1) * C(N-1,R/2-1)
        //                    / C(M+N,N) for R even;
        //
        //                  = ( C(M-1,(R-1)/2) * C(N-1,(R-3)/2 )
        //                    + C(M-1,(R-3)/2) * C(N-1,(R-1)/2 )
        //                    ) / C(M+N,N) for R odd.
        //
        //    Note that the maximum number of runs for a given M and N is:
        //
        //      M + N,                if M = N
        //      2 * min ( M, N ) + 1  otherwise
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Kalimutha Krishnamoorthy,
        //    Handbook of Statistical Distributions with Applications,
        //    Chapman and Hall, 2006,
        //    ISBN: 1-58488-635-8,
        //    LC: QA273.6.K75.
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //
        //    Input, int R, the number of runs.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;

            if (m < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("RUN_PDF - Fatal error!");
                Console.WriteLine("  M must be at least 0.");
                Console.WriteLine("  The input value of M = " + m + "");
                return(1);
            }

            if (n < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("RUN_PDF - Fatal error!");
                Console.WriteLine("  N must be at least 0.");
                Console.WriteLine("  The input value of N = " + n + "");
                return(1);
            }

            if (n + m <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("RUN_PDF - Fatal error!");
                Console.WriteLine("  M+N must be at least 1.");
                Console.WriteLine("  The input value of M+N = " + m + n + "");
                return(1);
            }

            //
            //  If all the symbols are of one type, there is always 1 run.
            //
            if (m == 0 || n == 0)
            {
                if (r == 1)
                {
                    pdf = 1.0;
                }
                else
                {
                    pdf = 0.0;
                }

                return pdf;
            }

            //
            //  Take care of extreme values of R.
            //
            if (r < 2 || m + n < r)
            {
                pdf = 0.0;
                return pdf;
            }

            //
            //  The normal cases.
            //
            if ((r % 2) == 0)
            {
                pdf = (double) (2 * typeMethods.i4_choose(m - 1, (r / 2) - 1)
                                  * typeMethods.i4_choose(n - 1, (r / 2) - 1))
                      / (double) (typeMethods.i4_choose(m + n, n));
            }
            else
            {
                pdf = (double) (typeMethods.i4_choose(m - 1, (r - 1) / 2)
                                * typeMethods.i4_choose(n - 1, (r - 3) / 2)
                                + typeMethods.i4_choose(m - 1, (r - 3) / 2)
                                * typeMethods.i4_choose(n - 1, (r - 1) / 2))
                      / (double) (typeMethods.i4_choose(m + n, n));
            }

            return pdf;
        }

        public static int runs_sample(int m, int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RUNS_SAMPLE samples the Runs PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int RUNS_SAMPLE, the number of runs.
        //
        {
            int[] a = runs_simulate(m, n, ref seed);

            int r = typeMethods.i4vec_run_count(m + n, a);
            
            return r;
        }

        public static int[] runs_simulate(int m, int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RUNS_SIMULATE simulates a case governed by the Runs PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int RUNS_SIMULATE[M+N], a sequence of M 0's and N 1's chosen
        //    uniformly at random.
        //
        {
            int j;
            int k;

            int[] a = new int[m + n];

            for (int i = 0; i < m; i++)
            {
                a[i] = 0;
            }

            for (int i = m; i < m + n; i++)
            {
                a[i] = 1;
            }

            for (int i = 1; i <= m + n - 1; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, m + n, ref seed);

                k = a[i - 1];
                a[i - 1] = a[j - 1];
                a[j - 1] = k;
            }

            return a;
        }

        public static double runs_variance(int m, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RUNS_VARIANCE returns the variance of the Runs PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the parameters of the PDF.
        //
        //    Output, double RUNS_VARIANCE, the variance of the PDF.
        //
        {
            double variance = (double) (2 * m * n * (2 * m * n - m - n))
                              / (double) ((m + n) * (m + n) * (m + n - 1));

            return variance;
        }
    }
}