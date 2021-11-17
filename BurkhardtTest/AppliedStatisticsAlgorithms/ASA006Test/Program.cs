using System;
using Burkardt.AppliedStatistics;

namespace ASA006Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA006_TEST.
        //
        //  Discussion:
        //
        //    ASA006_TEST tests the ASA006 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine();
        Console.WriteLine("ASA006_TEST:");
        Console.WriteLine("  Test the ASA006 library.");

        test01 ( );
        test02 ( );
        test03 ( );

        Console.WriteLine();
        Console.WriteLine("ASA006_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine();
    }


    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of CHOLESKY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 15;

        double[] a = new double [N_MAX * (N_MAX + 1) / 2];
        int ifault = 0;
        int nullty = 0;
        double[] u = new double [N_MAX * (N_MAX + 1) / 2];
        double[] ufull = new double [N_MAX * N_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  CHOLESKY computes the Cholesky factorization");
        Console.WriteLine("  of a positive definite symmetric matrix.");
        Console.WriteLine("  A compressed storage format is used");
        Console.WriteLine("");
        Console.WriteLine("  Here we look at the matrix A which is");
        Console.WriteLine("  N+1 on the diagonal and");
        Console.WriteLine("  N   on the off diagonals.");

        for (int n = 1; n <= N_MAX; n++)
        {
            int nn = n * (n + 1) / 2;
            //
            //  Set A to the lower triangle of the matrix which is N+1 on the diagonal
            //  and N on the off diagonals.
            //
            int k = 0;
            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j < i; j++)
                {
                    a[k] = n;
                    k += 1;
                }

                a[k] = n + 1;
                k += 1;
            }

            Algorithms.cholesky(a, n, nn, ref u, ref nullty, ref ifault);

            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Matrix nullity NULLTY = " + nullty + "");

            k = 0;
            for (int j = 1; j <= n; j++)
            {
                for (int i = 1; i <= j; i++)
                {
                    ufull[i - 1 + (j - 1) * n] = u[k];
                    k += 1;
                }

                for (int i = j + 1; i <= n; i++)
                {
                    ufull[i - 1 + (j - 1) * n] = 0.0;
                }
            }

            //
            //  Compute U' * U and compare to A.
            //
            k = 0;
            double diff = 0.0;
            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    double utu = 0.0;
                    for (int l = 1; l <= n; l++)
                    {
                        utu += ufull[l - 1 + (i - 1) * n] * ufull[l - 1 + (j - 1) * n];
                    }

                    diff += (a[k] - utu) * (a[k] - utu);
                    k += 1;
                }
            }

            diff = Math.Sqrt(diff);

            Console.WriteLine("  RMS ( A - U'*U ) = " + diff + "");
        }
    }


    private static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 demonstrates the use of CHOLESKY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 15;

        double[] a = new double [N_MAX * (N_MAX + 1) / 2];
        int ifault = 0;
        int nullty = 0;
        double[] u = new double[N_MAX * (N_MAX + 1) / 2];
        double[] ufull = new double [N_MAX * N_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  CHOLESKY computes the Cholesky factorization");
        Console.WriteLine("  of a positive definite symmetric matrix.");
        Console.WriteLine("  A compressed storage format is used");
        Console.WriteLine("");
        Console.WriteLine("  Here we look at the Hilbert matrix");
        Console.WriteLine("  A(I,J) = 1/(I+J-1)");
        Console.WriteLine("");
        Console.WriteLine("  For this matrix, we expect errors to grow quickly.");

        for (int n = 1; n <= N_MAX; n++)
        {
            int nn = n * (n + 1) / 2;
            //
            //  Set A to the Hilbert matrix.
            //
            int k = 0;
            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    a[k] = 1.0 / (i + j - 1);
                    k += 1;
                }
            }

            Algorithms.cholesky(a, n, nn, ref u, ref nullty, ref ifault);

            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Maxtrix nullity NULLTY = " + nullty + "");

            k = 0;
            for (int j = 1; j <= n; j++)
            {
                for (int i = 1; i <= j; i++)
                {
                    ufull[i - 1 + (j - 1) * n] = u[k];
                    k += 1;
                }

                for (int i = j + 1; i <= n; i++)
                {
                    ufull[i - 1 + (j - 1) * n] = 0.0;
                }
            }

            //
            //  Compute U' * U and compare to A.
            //
            k = 0;
            double diff = 0.0;
            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    double utu = 0.0;
                    for (int l = 1; l <= n; l++)
                    {
                        utu += ufull[l - 1 + (i - 1) * n] * ufull[l - 1 + (j - 1) * n];
                    }

                    diff += (a[k] - utu) * (a[k] - utu);
                    k += 1;
                }
            }

            diff = Math.Sqrt(diff);

            Console.WriteLine("  RMS ( A - U'*U ) = " + diff + "");
        }
    }


    private static void test03()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 demonstrates the use of SUBCHL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 15;
        int NN_MAX = N_MAX * (N_MAX + 1) / 2;

        double[] a = new double[NN_MAX];
        int[] b = new int[N_MAX];
        int ifault = 0;
        int nullty = 0;
        double[] u = new double[NN_MAX];
        double[] ufull = new double[N_MAX * N_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  SUBCHL computes the Cholesky factor");
        Console.WriteLine("  of a submatrix");
        Console.WriteLine("  of a positive definite symmetric matrix.");
        Console.WriteLine("  A compressed storage format is used.");
        Console.WriteLine("");
        Console.WriteLine("  Here we look at the Hilbert matrix");
        Console.WriteLine("  A(I,J) = 1/(I+J-1).");
        Console.WriteLine("");
        Console.WriteLine("  For this particular matrix, we expect the");
        Console.WriteLine("  errors to grow rapidly.");
        //
        //  Set A to the N_MAX order Hilbert matrix.
        //
        int k = 0;
        for (int i = 1; i <= N_MAX; i++)
        {
            for (int j = 1; j <= i; j++)
            {

                a[k] = 1.0 / (i + j - 1);
                k += 1;
            }
        }

        for (int n = 1; n <= N_MAX; n++)
        {
            for (int i = 1; i <= n; i++)
            {
                b[i - 1] = i;
            }

            double det = 0;
            Algorithms.subchl(a, b, n, ref u, ref nullty, ref ifault, NN_MAX, ref det);

            Console.WriteLine("");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Maxtrix nullity NULLTY = " + nullty + "");
            Console.WriteLine("  Matrix determinant DET = " + det + "");

            k = 0;
            for (int j = 1; j <= n; j++)
            {
                for (int i = 1; i <= j; i++)
                {
                    k += 1;
                    ufull[i - 1 + (j - 1) * n] = u[k - 1];
                }

                for (int i = j + 1; i <= n; i++)
                {
                    ufull[i - 1 + (j - 1) * n] = 0.0;
                }
            }

            //
            //  Compute U' * U and compare to A.
            //
            k = 0;
            double diff = 0.0;
            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= i; j++)
                {
                    k += 1;
                    double utu = 0.0;
                    for (int l = 1; l <= n; l++)
                    {
                        utu += ufull[l - 1 + (i - 1) * n] * ufull[l - 1 + (j - 1) * n];
                    }

                    diff += Math.Pow(a[k - 1] - utu, 2);
                }
            }

            diff = Math.Sqrt(diff);
            Console.WriteLine("  RMS ( A - U'*U ) = " + diff + "");
        }
    }

}