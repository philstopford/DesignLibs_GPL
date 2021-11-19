using System;
using Burkardt.SolveNS;
using Burkardt.Uniform;

namespace RestartedGeneralizedMinimumResidualTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MGMRES_TEST.
        //
        //  Discussion:
        //
        //    MGMRES_TEST tests the MGMRES library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("MGMRES_TEST:");
        Console.WriteLine("  Test the MGMRES library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("MGMRES_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests MGMRES_ST on the simple -1,2-1 matrix.
        //
        //  Discussion:
        //
        //    This is a very weak test, since the matrix has such a simple
        //    structure, is diagonally dominant (though not strictly), 
        //    and is symmetric.
        //
        //    To make the matrix bigger, simply increase the value of N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;
        int NZ_NUM = 3 * N - 2;

        double[] a = new double[NZ_NUM];
        int i;
        int[] ia = new int[NZ_NUM];
        int itr_max = 0;
        int[] ja = new int[NZ_NUM];
        int k;
        int mr = 0;
        int n = N;
        int nz_num = NZ_NUM;
        double[] rhs = new double[N];
        int test;
        double tol_abs = 0;
        double tol_rel = 0;
        double x_error = 0;
        double[] x_estimate = new double[N];
        double[] x_exact = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test MGMRES_ST on the simple -1,2-1 matrix.");
        //
        //  Set the matrix.
        //  Note that we use zero based index values in IA and JA.
        //
        k = 0;

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case > 0:
                    ia[k] = i;
                    ja[k] = i - 1;
                    a[k] = -1.0;
                    k += 1;
                    break;
            }

            ia[k] = i;
            ja[k] = i;
            a[k] = 2.0;
            k += 1;

            if (i < n - 1)
            {
                ia[k] = i;
                ja[k] = i + 1;
                a[k] = -1.0;
                k += 1;
            }

        }

        //
        //  Set the right hand side:
        //
        for (i = 0; i < n - 1; i++)
        {
            rhs[i] = 0.0;
        }

        rhs[N - 1] = n + 1;
        //
        //  Set the exact solution.
        //
        for (i = 0; i < n; i++)
        {
            x_exact[i] = i + 1;
        }

        for (test = 1; test <= 3; test++)
        {
            //
            //  Set the initial solution estimate.
            //
            for (i = 0; i < n; i++)
            {
                x_estimate[i] = 0.0;
            }

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            switch (test)
            {
                case 1:
                    itr_max = 1;
                    mr = 20;
                    break;
                case 2:
                    itr_max = 2;
                    mr = 10;
                    break;
                case 3:
                    itr_max = 5;
                    mr = 4;
                    break;
            }

            tol_abs = 1.0E-08;
            tol_rel = 1.0E-08;

            Console.WriteLine("");
            Console.WriteLine("  Test " + test + "");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Inner iteration limit = " + mr + "");
            Console.WriteLine("  Outer iteration limit = " + itr_max + "");
            Console.WriteLine("  Initial X_ERROR = " + x_error + "");

            RestartedGeneralizedMinimumResidual.mgmres_st(n, nz_num, ia, ja, a, ref x_estimate, rhs, itr_max, mr,
                tol_abs, tol_rel);

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            Console.WriteLine("  Final X_ERROR = " + x_error + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests MGMRES_ST on a 9 by 9 matrix.
        //
        //  Discussion:
        //
        //    A = 
        //      2  0  0 -1  0  0  0  0  0
        //      0  2 -1  0  0  0  0  0  0
        //      0 -1  2  0  0  0  0  0  0
        //     -1  0  0  2 -1  0  0  0  0
        //      0  0  0 -1  2 -1  0  0  0
        //      0  0  0  0 -1  2 -1  0  0
        //      0  0  0  0  0 -1  2 -1  0
        //      0  0  0  0  0  0 -1  2 -1
        //      0  0  0  0  0  0  0 -1  2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 9;
        int NZ_NUM = 23;

        double[] a =
        {
            2.0, -1.0,
            2.0, -1.0,
            -1.0, 2.0,
            -1.0, 2.0, -1.0,
            -1.0, 2.0, -1.0,
            -1.0, 2.0, -1.0,
            -1.0, 2.0, -1.0,
            -1.0, 2.0, -1.0,
            -1.0, 2.0
        };
        int i;
        int[] ia =
        {
            0, 0,
            1, 1,
            2, 2,
            3, 3, 3,
            4, 4, 4,
            5, 5, 5,
            6, 6, 6,
            7, 7, 7,
            8, 8
        };
        int itr_max;
        int[] ja =
        {
            0, 3,
            1, 2,
            1, 2,
            0, 3, 4,
            3, 4, 5,
            4, 5, 6,
            5, 6, 7,
            6, 7, 8,
            7, 8
        };
        int mr;
        int n = N;
        int nz_num = NZ_NUM;
        double[] rhs =
        {
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0
        };
        int seed = 123456789;
        int test;
        double tol_abs;
        double tol_rel;
        double x_error;
        double[] x_estimate = new double[1];
        double[] x_exact =
        {
            3.5,
            1.0,
            1.0,
            6.0,
            7.5,
            8.0,
            7.5,
            6.0,
            3.5
        };

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test MGMRES_ST on matrix that is not quite the -1,2,-1 matrix,");
        Console.WriteLine("  of order N = " + n + "");

        for (test = 1; test <= 2; test++)
        {
            switch (test)
            {
                case 1:
                {
                    Console.WriteLine("");
                    Console.WriteLine("  First try, use zero initial vector:");
                    Console.WriteLine("");

                    x_estimate = new double[n];
                    for (i = 0; i < n; i++)
                    {
                        x_estimate[i] = 0.0;
                    }

                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("  Second try, use random initial vector:");
                    Console.WriteLine("");

                    UniformRNG.r8vec_uniform_01(n, ref seed, ref x_estimate);
                    break;
            }

            //
            //  Set the initial solution estimate.
            //
            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            Console.WriteLine("  Before solving, X_ERROR = " + x_error + "");

            itr_max = 20;
            mr = n - 1;
            tol_abs = 1.0E-08;
            tol_rel = 1.0E-08;

            RestartedGeneralizedMinimumResidual.mgmres_st(n, nz_num, ia, ja, a, ref x_estimate, rhs, itr_max, mr,
                tol_abs, tol_rel);

            x_error = 0.0;
            for (i = 0; i < N; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            Console.WriteLine("  After solving, X_ERROR = " + x_error + "");

            Console.WriteLine("");
            Console.WriteLine("  Final solution estimate:");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + x_estimate[i].ToString().PadLeft(12) + "");
            }
        }
    }

    private static void test03()

        //******************************************************************************
        //
        //  Purpose:
        //
        //    TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.
        //
        //  Discussion:
        //
        //    This is a very weak test, since the matrix has such a simple
        //    structure, is diagonally dominant (though not strictly), 
        //    and is symmetric.
        //
        //    To make the matrix bigger, simply increase the value of N.
        //
        //    Note that PGMRES_ILU_CR expects the matrix to be stored using the
        //    sparse compressed row format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;
        int NZ_NUM = 3 * N - 2;

        double[] a = new double[NZ_NUM];
        int i;
        int[] ia = new int[N + 1];
        int itr_max = 0;
        int[] ja = new int[NZ_NUM];
        int k;
        int mr = 0;
        int n = N;
        int nz_num = NZ_NUM;
        double[] rhs = new double[N];
        int test;
        double tol_abs = 0;
        double tol_rel = 0;
        double x_error = 0;
        double[] x_estimate = new double[N];
        double[] x_exact = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.");
        //
        //  Set the matrix.
        //  Note that we use zero based index valuesin IA and JA.
        //
        k = 0;
        ia[0] = 0;

        Console.WriteLine("");
        Console.WriteLine("  ia[" + 0 + "] = " + ia[0] + "");
        for (i = 0; i < n; i++)
        {
            ia[i + 1] = ia[i];
            switch (i)
            {
                case > 0:
                    ia[i + 1] += 1;
                    ja[k] = i - 1;
                    a[k] = -1.0;
                    k += 1;
                    break;
            }

            ia[i + 1] += 1;
            ja[k] = i;
            a[k] = 2.0;
            k += 1;

            if (i < N - 1)
            {
                ia[i + 1] += 1;
                ja[k] = i + 1;
                a[k] = -1.0;
                k += 1;
            }

            Console.WriteLine("  ia[" + i + 1 + "] = " + ia[i + 1] + "");
        }

        //
        //  Set the right hand side:
        //
        for (i = 0; i < n - 1; i++)
        {
            rhs[i] = 0.0;
        }

        rhs[n - 1] = n + 1;
        //
        //  Set the exact solution.
        //
        for (i = 0; i < n; i++)
        {
            x_exact[i] = i + 1;
        }

        for (test = 1; test <= 3; test++)
        {
            //
            //  Set the initial solution estimate.
            //
            for (i = 0; i < n; i++)
            {
                x_estimate[i] = 0.0;
            }

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            switch (test)
            {
                case 1:
                    itr_max = 1;
                    mr = 20;
                    break;
                case 2:
                    itr_max = 2;
                    mr = 10;
                    break;
                case 3:
                    itr_max = 5;
                    mr = 4;
                    break;
            }

            tol_abs = 1.0E-08;
            tol_rel = 1.0E-08;

            Console.WriteLine("");
            Console.WriteLine("  Test " + test + "");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Inner iteration limit = " + mr + "");
            Console.WriteLine("  Outer iteration limit = " + itr_max + "");
            Console.WriteLine("  Initial X_ERROR = " + x_error + "");

            RestartedGeneralizedMinimumResidual.pmgmres_ilu_cr(n, nz_num, ia, ja, a, ref x_estimate, rhs, itr_max,
                mr,
                tol_abs, tol_rel);

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            Console.WriteLine("  Final X_ERROR = " + x_error + "");
        }
    }

    private static void test04()

        //******************************************************************************
        //
        //  Purpose:
        //
        //    TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int NZ_NUM = 9;

        double[] a =
        {
            1.0, 2.0, 1.0,
            2.0,
            3.0, 3.0,
            4.0,
            1.0, 5.0
        };
        int i;
        int[] ia = {0, 3, 4, 6, 7, 9};
        int itr_max = 0;
        int[] ja =
        {
            0, 3, 4,
            1,
            0, 2,
            3,
            1, 4
        };
        int mr = 0;
        int n = N;
        int nz_num = NZ_NUM;
        double[] rhs = {14.0, 4.0, 12.0, 16.0, 27.0};
        int test;
        double tol_abs = 0;
        double tol_rel = 0;
        double x_error = 0;
        double[] x_estimate = new double[N];
        double[] x_exact = {1.0, 2.0, 3.0, 4.0, 5.0};

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.");

        Console.WriteLine("");
        for (i = 0; i < n + 1; i++)
        {
            Console.WriteLine("  ia[" + i + "] = " + ia[i] + "");
        }

        for (test = 1; test <= 3; test++)
        {
            //
            //  Set the initial solution estimate.
            //
            for (i = 0; i < n; i++)
            {
                x_estimate[i] = 0.0;
            }

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            switch (test)
            {
                case 1:
                    itr_max = 1;
                    mr = 20;
                    break;
                case 2:
                    itr_max = 2;
                    mr = 10;
                    break;
                case 3:
                    itr_max = 5;
                    mr = 4;
                    break;
            }

            tol_abs = 1.0E-08;
            tol_rel = 1.0E-08;

            Console.WriteLine("");
            Console.WriteLine("  Test " + test + "");
            Console.WriteLine("  Matrix order N = " + n + "");
            Console.WriteLine("  Inner iteration limit = " + mr + "");
            Console.WriteLine("  Outer iteration limit = " + itr_max + "");
            Console.WriteLine("  Initial X_ERROR = " + x_error + "");

            RestartedGeneralizedMinimumResidual.pmgmres_ilu_cr(n, nz_num, ia, ja, a, ref x_estimate, rhs, itr_max,
                mr,
                tol_abs, tol_rel);

            x_error = 0.0;
            for (i = 0; i < n; i++)
            {
                x_error += Math.Pow(x_exact[i] - x_estimate[i], 2);
            }

            x_error = Math.Sqrt(x_error);

            Console.WriteLine("  Final X_ERROR = " + x_error + "");
        }
    }
}