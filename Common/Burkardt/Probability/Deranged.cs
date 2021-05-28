using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Deranged
    {
        public static double deranged_cdf(int x, int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_CDF evaluates the Deranged CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the maximum number of items in their correct places.
        //    0 <= X <= A.
        //
        //    Input, int A, the number of items.
        //    1 <= A.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x < 0 || a < x)
            {
                cdf = 0.0;
            }
            else
            {
                int sum2 = 0;
                int x2;
                for (x2 = 0; x2 <= x; x2++)
                {
                    int cnk = typeMethods.i4_choose(a, x2);
                    int dnmk = deranged_enum(a - x2);
                    sum2 = sum2 + cnk * dnmk;
                }

                cdf = (double) (sum2) / typeMethods.r8_factorial(a);
            }

            return cdf;
        }

        public static int deranged_cdf_inv(double cdf, int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_CDF_INV inverts the Deranged CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Input, int A, the number of items.
        //    1 <= A.
        //
        //    Output, int DERANGED_CDF_INV, the corresponding argument.
        //
        {
            int x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("DERANGED_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return(1);
            }

            double cdf2 = 0.0;

            for (int x2 = 0; x2 <= a; x2++)
            {
                double pdf = deranged_pdf(x2, a);

                cdf2 = cdf2 + pdf;

                if (cdf <= cdf2)
                {
                    x = x2;
                    return x;
                }

            }

            x = a;

            return x;
        }

        public static bool deranged_check(int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_CHECK checks the parameter of the Deranged PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the total number of items.
        //    1 <= A.
        //
        //    Output, bool DERANGED_CHECK, is true if the parameters are legal.
        //
        {
            if (a < 1)
            {
                Console.WriteLine(" ");
                Console.WriteLine("DERANGED_CHECK - Warning!");
                Console.WriteLine("  A < 1.");
                return false;
            }

            return true;
        }

        public static int deranged_enum(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_ENUM returns the number of derangements of N objects.
        //
        //  Discussion:
        //
        //    A derangement of N objects is a permutation with no fixed
        //    points.  If we symbolize the permutation operation by "P",
        //    then for a derangment, P(I) is never equal to I.
        //
        //  Recursion:
        //
        //      D(0) = 1
        //      D(1) = 0
        //      D(2) = 1
        //      D(N) = (N-1) * ( D(N-1) + D(N-2) )
        //
        //    or
        //
        //      D(0) = 1
        //      D(1) = 0
        //      D(N) = N * D(N-1) + (-1)^N
        //
        //    D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
        //
        //    Based on the inclusion/exclusion law.
        //
        //    D(N) is the number of ways of placing N non-attacking rooks on
        //    an N by N chessboard with one diagonal deleted.
        //
        //    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
        //
        //    The number of permutations with exactly K items in the right
        //    place is COMB(N,K) * D(N-K).
        //
        //  First values:
        //
        //     N         D(N)
        //     0           1
        //     1           0
        //     2           1
        //     3           2
        //     4           9
        //     5          44
        //     6         265
        //     7        1854
        //     8       14833
        //     9      133496
        //    10     1334961
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects to be permuted.
        //
        //    Output, int DERANGED_ENUM, the number of derangements of N objects.
        //
        {
            int dn;

            if (n < 0)
            {
                dn = 0;
            }
            else if (n == 0)
            {
                dn = 1;
            }
            else if (n == 1)
            {
                dn = 0;
            }
            else if (n == 2)
            {
                dn = 1;
            }
            else
            {
                int dnm1 = 0;
                dn = 1;

                for (int i = 3; i <= n; i++)
                {
                    int dnm2 = dnm1;
                    dnm1 = dn;
                    dn = (i - 1) * (dnm1 + dnm2);
                }

            }

            return dn;
        }

        public static double deranged_mean(int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_MEAN returns the mean of the Deranged CDF.
        //
        //  Discussion:
        //
        //    The mean is computed by straightforward summation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of items.
        //    1 <= A.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = 0.0;
            for (int x = 0; x <= a; x++)
            {
                double pdf = deranged_pdf(x, a);
                mean = mean + pdf * x;
            }

            return mean;
        }

        public static double deranged_pdf(int x, int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_PDF evaluates the Deranged PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) is the probability that exactly X items will occur in
        //    their proper place after a random permutation of A items.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the number of items in their correct places.
        //    0 <= X <= A.
        //
        //    Input, int A, the total number of items.
        //    1 <= A.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < 0 || a < x)
            {
                pdf = 0.0;
            }
            else
            {
                int cnk = typeMethods.i4_choose(a, x);
                int dnmk = deranged_enum(a - x);
                pdf = (double) (cnk * dnmk) / typeMethods.r8_factorial(a);
            }

            return pdf;
        }

        public static int deranged_sample(int a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_SAMPLE samples the Deranged PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of items.
        //    1 <= A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int DERANGED_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            int x = deranged_cdf_inv(cdf, a);

            return x;
        }

        public static double deranged_variance(int a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DERANGED_VARIANCE returns the variance of the Deranged CDF.
        //
        //  Discussion:
        //
        //    The variance is computed by straightforward summation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the number of items.
        //    1 <= A.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            double mean = deranged_mean(a);

            double variance = 0.0;
            for (int x = 0; x <= a; x++)
            {
                double pdf = deranged_pdf(x, a);
                variance = variance + pdf * Math.Pow((x - mean), 2);
            }

            return variance;
        }
    }
}