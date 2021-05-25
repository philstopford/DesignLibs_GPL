using System;
using Burkardt.AppliedStatistics;

namespace ASA007Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA007_TEST.
        //
        //  Discussion:
        //
        //    ASA007_TEST tests the ASA007 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA007_TEST:");
            Console.WriteLine("  Test the ASA007 library.");

            test01 ( );
            test02 ( );

            Console.WriteLine("");
            Console.WriteLine("ASA007_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of SYMINV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int N_MAX = 15;

            double[] a = new double[(N_MAX * (N_MAX + 1)) / 2];
            double[] afull = new double[N_MAX * N_MAX];
            double[] c = new double[(N_MAX * (N_MAX + 1)) / 2];
            double[] cfull = new double[N_MAX * N_MAX];
            int ifault = 0;
            int nullty = 0;
            int k;
            double[] w = new double[N_MAX];

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  SYMINV computes the inverse of a positive");
            Console.WriteLine("  definite symmetric matrix.");
            Console.WriteLine("  A compressed storage format is used");
            Console.WriteLine("");
            Console.WriteLine("  Here we look at the matrix A which is");
            Console.WriteLine("  N+1 on the diagonal and");
            Console.WriteLine("  N   on the off diagonals.");

            for (int n = 1; n <= N_MAX; n++)
            {
                //
                //  Set A to the lower triangle of the matrix which is N+1 on the diagonal
                //  and N on the off diagonals.
                //
                k = 0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 1; j < i; j++)
                    {
                        a[k] = (double) (n);
                        k = k + 1;
                    }

                    a[k] = (double) (n + 1);
                    k = k + 1;
                }

                Algorithms.syminv(a, n, ref c, w, ref nullty, ref ifault);

                Console.WriteLine("");
                Console.WriteLine("  Matrix order N = " + n + "");
                Console.WriteLine("  Matrix nullity NULLTY = " + nullty + "");

                k = 0;
                for (int j = 1; j <= n; j++)
                {
                    for (int i = 1; i < j; i++)
                    {
                        cfull[i - 1 + (j - 1) * n] = c[k];
                        cfull[j - 1 + (i - 1) * n] = c[k];
                        k = k + 1;
                    }

                    cfull[j - 1 + (j - 1) * n] = c[k];
                    k = k + 1;
                }

                k = 0;
                for (int j = 1; j <= n; j++)
                {
                    for (int i = 1; i < j; i++)
                    {
                        afull[i - 1 + (j - 1) * n] = a[k];
                        afull[j - 1 + (i - 1) * n] = a[k];
                        k = k + 1;
                    }

                    afull[j - 1 + (j - 1) * n] = a[k];
                    k = k + 1;
                }

                //
                //  Compute C * A - I.
                //
                double diff = 0.0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 1; j <= i; j++)
                    {
                        double cta = 0.0;
                        for (k = 1; k <= n; k++)
                        {
                            cta = cta + cfull[i - 1 + (k - 1) * n] * afull[k - 1 + (j - 1) * n];
                        }

                        if (i == j)
                        {
                            diff = diff + Math.Pow(1.0 - cta, 2);
                        }
                        else
                        {
                            diff = diff + cta * cta;
                        }
                    }
                }

                diff = Math.Sqrt(diff);

                Console.WriteLine("  RMS ( C * A - I ) = " + diff + "");
            }
        }


        static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 demonstrates the use of SYMINV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int N_MAX = 15;

            double[] a = new double[(N_MAX * (N_MAX + 1)) / 2];
            double[] afull = new double[N_MAX * N_MAX];
            double[] c = new double[(N_MAX * (N_MAX + 1)) / 2];
            double[] cfull = new double[N_MAX * N_MAX];
            int ifault = 0;
            int nullty = 0;
            double[] w = new double[N_MAX];

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  SYMINV computes the inverse of a positive");
            Console.WriteLine("  definite symmetric matrix.");
            Console.WriteLine("  A compressed storage format is used");
            Console.WriteLine("");
            Console.WriteLine("  Here we look at the Hilbert matrix");
            Console.WriteLine("  A(I,J) = 1/(I+J-1)");
            Console.WriteLine("");
            Console.WriteLine("  For this matrix, we expect errors to grow quickly.");

            for (int n = 1; n <= N_MAX; n++)
            {
                //
                //  Set A to the Hilbert matrix.
                //
                int k = 0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 1; j <= i; j++)
                    {
                        a[k] = 1.0 / (double) (i + j - 1);
                        k = k + 1;
                    }
                }

                Algorithms.syminv(a, n, ref c, w, ref nullty, ref ifault);

                Console.WriteLine("");
                Console.WriteLine("  Matrix order N = " + n + "");
                Console.WriteLine("  Maxtrix nullity NULLTY = " + nullty + "");

                k = 0;
                for (int j = 1; j <= n; j++)
                {
                    for (int i = 1; i < j; i++)
                    {
                        cfull[i - 1 + (j - 1) * n] = c[k];
                        cfull[j - 1 + (i - 1) * n] = c[k];
                        k = k + 1;
                    }

                    cfull[j - 1 + (j - 1) * n] = c[k];
                    k = k + 1;
                }

                k = 0;
                for (int j = 1; j <= n; j++)
                {
                    for (int i = 1; i < j; i++)
                    {
                        afull[i - 1 + (j - 1) * n] = a[k];
                        afull[j - 1 + (i - 1) * n] = a[k];
                        k = k + 1;
                    }

                    afull[j - 1 + (j - 1) * n] = a[k];
                    k = k + 1;
                }

                //
                //  Compute C * A - I.
                //
                double diff = 0.0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 1; j <= i; j++)
                    {
                        double cta = 0.0;
                        for (k = 1; k <= n; k++)
                        {
                            cta = cta + cfull[i - 1 + (k - 1) * n] * afull[k - 1 + (j - 1) * n];
                        }

                        if (i == j)
                        {
                            diff = diff + Math.Pow(1.0 - cta, 2);
                        }
                        else
                        {
                            diff = diff + cta * cta;
                        }
                    }
                }

                diff = Math.Sqrt(diff);

                Console.WriteLine("  RMS ( C * A - I ) = " + diff + "");
            }
        }

    }
}