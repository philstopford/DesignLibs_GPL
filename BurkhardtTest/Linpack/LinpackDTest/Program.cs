using System;
using Burkardt.Linpack;
using Burkardt.MatrixNS;
using Burkardt.Uniform;

namespace LinpackDTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LINPACK_D_TEST.
        //
        //  Discussion:
        //
        //    LINPACK_D_TEST tests the LINPACK_D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LINPACK_D_TEST");
            
        Console.WriteLine("  Test the LINPACK_D library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();
        test13();
        test14();
        test15();
        test16();
        test17();
        test18();
        test19();

        test20();
        test21();
        test22();
        dqrdc_test();
        dqrsl_test();
        test24();
        test25();
        test26();
        test27();
        dsvdc_test();
        test29();

        test30();
        test31();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LINPACK_D_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DCHDC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] b = new double[LDA * N];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;
        int k;
        double[] work = new double[N];
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  For a general matrix,");
        Console.WriteLine("  DCHDC computes the Cholesky decomposition.");
        Console.WriteLine("");
        Console.WriteLine("  The number of equations is N = " + N + "");
        //
        //  Set the values of the matrix A.
        //
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= N; i++)
            {
                if (i == j - 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else if (i == j)
                {
                    a[i - 1 + (j - 1) * LDA] = 2.0;
                }
                else if (i == j + 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else
                {
                    a[i - 1 + (j - 1) * LDA] = 0.0;
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Decompose the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Decompose the matrix.");

        job = 0;

        for (i = 0; i < N; i++)
        {
            ipvt[i] = 0;
        }

        info = DCHDC.dchdc(ref a, LDA, N, work, ref ipvt, job);

        if (info != N)
        {
            Console.WriteLine("");
            Console.WriteLine("  DCHDC returned INFO = " + info + "");
            Console.WriteLine("  This means the matrix is not positive definite.");
            return;
        }

        //
        //  Zero out the lower diagonal.
        //
        for (i = 2; i <= N; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        //
        //  Print the factorization.
        //
        Console.WriteLine("");
        Console.WriteLine("  The Cholesky factor U:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Compute the Cholesky product.
        //
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                b[i - 1 + (j - 1) * LDA] = 0.0;
                for (k = 1; k <= N; k++)
                {
                    b[i - 1 + (j - 1) * LDA] += a[k - 1 + (i - 1) * LDA] * a[k - 1 + (j - 1) * LDA];
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The product U' * U:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + b[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests DCHEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = N;
        int NZ = 1;

        double[] a = new double[LDA * N];
        double[] b = new double[LDA * N];
        double[] c = new double[N];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;
        int k;
        int l;
        double[] s = new double[N];
        double[] work = new double[N];
        double[] z = new double[N];
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For a general matrix,");
        Console.WriteLine("  DCHEX can shift columns in a Cholesky factorization.");
        Console.WriteLine("");
        Console.WriteLine("  The number of equations is N = " + N + "");
        //
        //  Set the values of the matrix A.
        //
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= N; i++)
            {
                if (i == j - 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else if (i == j)
                {
                    a[i - 1 + (j - 1) * LDA] = 2.0;
                }
                else if (i == j + 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else
                {
                    a[i - 1 + (j - 1) * LDA] = 0.0;
                }
            }
        }

        for (i = 1; i <= N; i++)
        {
            z[i - 1] = i;
        }

        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The vector Z:");
        Console.WriteLine("");

        cout = "";
        for (i = 1; i <= N; i++)
        {
            cout += "  " + z[i - 1].ToString().PadLeft(12);
        }

        //
        //  Decompose the matrix.
        //
        Console.WriteLine(cout);
        Console.WriteLine("  Decompose the matrix.");

        job = 0;

        for (i = 0; i < N; i++)
        {
            ipvt[i] = 0;
        }

        info = DCHDC.dchdc(ref a, LDA, N, work, ref ipvt, job);

        if (info != N)
        {
            Console.WriteLine("");
            Console.WriteLine("  DCHDC returned INFO = " + info + "");
            Console.WriteLine("  This means the matrix is not positive definite.");
            return;
        }

        //
        //  Zero out the lower diagonal.
        //
        for (i = 2; i <= N; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        //
        //  Print the factorization.
        //
        Console.WriteLine("");
        Console.WriteLine("  The Cholesky factor U:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Right circular shift columns L through K.
        //
        k = 1;
        l = 3;

        Console.WriteLine("");
        Console.WriteLine("  Right circular shift columns K  = " + k +
                          " through L = " + l + "");

        job = 1;
        DCHEX.dchex(ref a, LDA, N, k, l, ref z, N, NZ, ref c, ref s, job);
        //
        //  Left circular shift columns K+1 through L.
        //
        Console.WriteLine("");
        Console.WriteLine("  Left circular shift columns K+1 = " + k + 1 +
                          " through L = " + l + "");

        job = 2;
        DCHEX.dchex(ref a, LDA, N, k + 1, l, ref z, N, NZ, ref c, ref s, job);
        //
        //  Print the factorization.
        //
        Console.WriteLine("");
        Console.WriteLine("  The shifted Cholesky factor U:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The shifted vector Z:");
        Console.WriteLine("");

        cout = "";
        for (i = 1; i <= N; i++)
        {
            cout += "  " + z[i - 1].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        //
        // Compute the Cholesky product.
        //
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                b[i - 1 + (j - 1) * LDA] = 0.0;
                for (k = 1; k <= N; k++)
                {
                    b[i - 1 + (j - 1) * LDA] += a[k - 1 + (i - 1) * LDA] * a[k - 1 + (j - 1) * LDA];
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The shifted product U' * U:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + b[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests DCHUD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int P = 20;
        int LDR = P;
        int NZ = 1;

        double[] b = new double[P];
        double[] c = new double[P];
        int i;
        int j;
        int job;
        double[] r = new double[LDR * P];
        double[] rho = new double[NZ];
        double[] row = new double[P];
        double[] s = new double[P];
        int seed;
        double[] x = new double[P];
        double[] y = new double[NZ];
        double[] z = new double[P * NZ];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For a general matrix,");
        Console.WriteLine("  DCHUD updates a Cholesky decomposition.");
        Console.WriteLine("");
        Console.WriteLine("  In this example, we use DCHUD to solve a");
        Console.WriteLine("  least squares problem R * b = z.");
        Console.WriteLine("");
        Console.WriteLine("  The number of equations is P = " + P + "");
        //
        //  Initialize.
        //
        for (j = 1; j <= P; j++)
        {
            for (i = 1; i <= P; i++)
            {
                r[i - 1 + (j - 1) * LDR] = 0.0;
            }
        }

        for (j = 1; j <= NZ; j++)
        {
            for (i = 1; i <= P; i++)
            {
                z[i - 1 + (j - 1) * P] = 0.0;
            }
        }

        for (i = 1; i <= P; i++)
        {
            x[i - 1] = i;
        }

        //
        //  Use DCHUD to form R, Z and RHO by adding X and Y a row at a time.
        //  X is a row of the least squares matrix and Y the right hand side.
        //
        seed = 123456789;

        for (i = 1; i <= P; i++)
        {
            row = UniformRNG.r8mat_uniform_01(1, P, ref seed);
            y[0] = 0.0;
            for (j = 1; j <= P; j++)
            {
                y[0] += row[j - 1] * x[j - 1];
            }

            rho[0] = 0.0;
            DCHUD.dchud(ref r, LDR, P, row, ref z, P, NZ, y, ref rho, ref c, ref s);
        }

        //
        //  Generate the least squares solution, b = inverse ( R ) * Z.
        //
        for (j = 1; j <= NZ; j++)
        {
            for (i = 1; i <= P; i++)
            {
                b[i - 1] = z[i - 1 + (j - 1) * P];
            }

            job = 01;

            DTRSL.dtrsl(r, LDR, P, ref b, job);

            Console.WriteLine("");
            Console.WriteLine("  Solution vector # " + j + "");
            Console.WriteLine("  (Should be (1,2,3...,n))");
            Console.WriteLine("");

            for (i = 1; i <= P; i++)
            {
                if (i <= 5 || P - 5 < i)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(6)
                                           + "  " + b[i - 1].ToString().PadLeft(14) + "");
                }

                switch (i)
                {
                    case 5:
                        Console.WriteLine("  ......  ..............");
                        break;
                }
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests DGBCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        int ML = 1;
        int MU = 1;
        int LDA = 2 * ML + MU + 1;

        double[] a = new double[LDA * N];
        int[] ipivot = new int[N];
        int j;
        int m;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGBCO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        m = ML + MU + 1;
        Console.WriteLine("  The bandwidth of the matrix is " + m + "");

        for (j = 1; j <= N; j++)
        {
            a[m - 2 + (j - 1) * LDA] = -1.0;
            a[m - 1 + (j - 1) * LDA] = 2.0;
            a[m + (j - 1) * LDA] = -1.0;
        }

        //
        //  Estimate the condition.
        //
        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition.");

        rcond = DGBCO.dgbco(ref a, LDA, N, ML, MU, ref ipivot, z);

        Console.WriteLine("");
        Console.WriteLine("  Estimated reciprocal condition = " + rcond + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests DGBFA and DGBSL.
        //
        //  Discussion:
        //
        //    The problem solved here is a larger version of this one:
        //
        //    Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
        //                (-1  2 -1  0  0)                          (0)
        //                ( 0 -1  2 -1  0)                          (0)
        //                ( 0  0 -1  2 -1)                          (0)
        //                ( 0  0  0 -1  2)                          (1)
        //
        //
        //    Solution is   (1)
        //                  (1)
        //                  (1)
        //                  (1)
        //                  (1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local Parameters:
        //
        //    N is the number of equations.
        //
        //    ML is the number of subdiagonals,
        //    MU the number of superdiagonals.
        //
        //    LDA is the leading dimension of the array used to store the
        //    matrix, which must be at least 2*ML+MU+1.
        //
    {
        int N = 10;
        int ML = 1;
        int MU = 1;
        int LDA = 2 * ML + MU + 1;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int info;
        int[] ipivot = new int[N];
        int j;
        int job;
        int m;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGBFA computes the LU factors,");
        Console.WriteLine("  DGBSL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the right hand side B.
        //
        b[0] = 1.0;
        for (i = 2; i <= N - 1; i++)
        {
            b[i - 1] = 0.0;
        }

        b[N - 1] = 1.0;
        //
        //  Set the matrix A.
        //
        m = ML + MU + 1;
        Console.WriteLine("  The bandwidth of the matrix is " + m + "");

        for (j = 1; j <= N; j++)
        {
            a[m - 2 + (j - 1) * LDA] = -1.0;
            a[m - 1 + (j - 1) * LDA] = 2.0;
            a[m + (j - 1) * LDA] = -1.0;
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DGBFA.dgbfa(ref a, LDA, N, ML, MU, ref ipivot);

        if (info != 0)
        {
            Console.WriteLine("  Error!  DGBFA returns INFO = " + info + "");
            return;
        }

        //
        //  Call DGBSL to solve the linear system.  The solution
        //  is returned in B, that is, it overwrites the right hand side.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        job = 0;
        DGBSL.dgbsl(a, LDA, N, ML, MU, ipivot, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,1,1,1,1,...,1,1))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests DGBFA and DGBDI.
        //
        //  Discussion:
        //
        //    Matrix A is ( 2 -1  0  0  0)
        //                (-1  2 -1  0  0)
        //                ( 0 -1  2 -1  0)
        //                ( 0  0 -1  2 -1)
        //                ( 0  0  0 -1  2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local Parameters:
        //
        //    N is the number of equations.
        //
        //    ML is the number of subdiagonals,
        //    MU the number of superdiagonals.
        //
        //    LDA is the leading dimension of the array used to store the
        //    matrix, which must be at least 2*ML+MU+1.
        //
    {
        const int N_MAX = 128;
        int ML = 1;
        int MU = 1;
        int LDA = 2 * ML + MU + 1;

        double[] a = new double[LDA * N_MAX];
        double[] det = new double[2];
        int i;
        int info;
        int[] ipivot = new int[N_MAX];
        int j;
        int m;
        int n;
        int n_log;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGBFA factors the matrix,");
        Console.WriteLine("  DGBDI computes the determinant as");
        Console.WriteLine("    det = MANTISSA * 10^EXPONENT");
        Console.WriteLine("");
        Console.WriteLine("  Find the determinant of the -1,2,-1 matrix");
        Console.WriteLine("  for N = 2, 4, 8, 16, 32, 64, 128.");
        Console.WriteLine("");
        Console.WriteLine("  (For this matrix, det ( A ) = N + 1.)");
        //
        //  Set the matrix A.
        //
        m = ML + MU + 1;
        Console.WriteLine("  The bandwidth of the matrix is " + m + "");
        Console.WriteLine("");
        Console.WriteLine("       N    Mantissa       Exponent");
        Console.WriteLine("");

        n = 1;

        for (n_log = 1; n_log <= 7; n_log++)
        {
            n = 2 * n;

            for (j = 1; j <= n; j++)
            {
                for (i = 1; i <= LDA; i++)
                {
                    a[i - 1 + (j - 1) * LDA] = 0.0;
                }
            }

            for (j = 1; j <= n; j++)
            {
                i = j;
                a[i - j + ML + MU + (j - 1) * LDA] = 2.0;
            }

            for (j = 2; j <= n; j++)
            {
                i = j - 1;
                a[i - j + ML + MU + (j - 1) * LDA] = -1.0;
            }

            for (j = 1; j <= n - 1; j++)
            {
                i = j + 1;
                a[i - j + ML + MU + (j - 1) * LDA] = -1.0;
            }

            info = DGBFA.dgbfa(ref a, LDA, n, ML, MU, ref ipivot);

            if (info != 0)
            {
                Console.WriteLine("  Error!  DGBFA returns INFO = " + info + "");
                return;
            }

            DGBDI.dgbdi(a, LDA, n, ML, MU, ipivot, ref det);

            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + det[0].ToString().PadLeft(14)
                                   + "  " + det[1].ToString().PadLeft(14) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests DGBFA and DGBSL.
        //
        //  Discussion:
        //
        //    DGBFA and DGBSL are for general banded matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;
        int ML = 25;
        int MU = 25;
        int LDA = 2 * ML + MU + 1;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int ihi;
        int ilo;
        int info;
        int[] ipivot = new int[N];
        int j;
        int job;
        int m;
        double temp;

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGBFA computes the LU factors,");
        Console.WriteLine("  DGBSL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to matrix A and right hand side B.
        //
        //  We want to try a problem with a significant bandwidth.
        //
        m = ML + MU + 1;
        Console.WriteLine("  The bandwidth of the matrix is " + m + "");

        for (j = 1; j <= N; j++)
        {
            ilo = Math.Max(1, j - MU);
            ihi = Math.Min(N, j + ML);

            temp = 0.0;
            for (i = ilo; i <= ihi; i++)
            {
                a[i - j + m - 1 + (j - 1) * LDA] = -1.0;
                temp -= 1.0;
            }

            temp += 1.0;
            a[m - 1 + (j - 1) * LDA] = 4.0 - temp;
            b[j - 1] = 4.0;
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DGBFA.dgbfa(ref a, LDA, N, ML, MU, ref ipivot);

        if (info != 0)
        {
            Console.WriteLine("  Error!  DGBFA returns INFO = " + info + "");
            return;
        }

        //
        //  Call DGBSL to solve the linear system.  The solution
        //  is returned in B, that is, it overwrites the right hand side.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        job = 0;
        DGBSL.dgbsl(a, LDA, N, ML, MU, ipivot, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,1,1,1,1,...,1,1))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 calls DGECO and DGESL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local Parameters:
        //
        //    LDA defines the maximum matrix size we will use.
        //
    {
        int LDA = 10;

        double[] a = new double[LDA * LDA];
        double[] b = new double[LDA];
        int i;
        int[] ipvt = new int[LDA];
        int job;
        int n;
        double rcond;
        double[] z = new double[LDA];

        n = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  For a general matrix,");
        Console.WriteLine("  DGECO computes the LU factors and computes");
        Console.WriteLine("  its reciprocal condition number;");
        Console.WriteLine("  DGESL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + n + "");
        //
        //  Set the values of the matrix A.
        //
        a[0 + 0 * LDA] = 1.0;
        a[0 + 1 * LDA] = 2.0;
        a[0 + 2 * LDA] = 3.0;

        a[1 + 0 * LDA] = 4.0;
        a[1 + 1 * LDA] = 5.0;
        a[1 + 2 * LDA] = 6.0;

        a[2 + 0 * LDA] = 7.0;
        a[2 + 1 * LDA] = 8.0;
        a[2 + 2 * LDA] = 0.0;
        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        rcond = DGECO.dgeco(ref a, LDA, n, ref ipvt, ref z);

        Console.WriteLine("  The reciprocal matrix condition number = " + rcond + "");

        switch (rcond + 1.0)
        {
            case 1.0:
                Console.WriteLine("  Error!  The matrix is nearly singular!");
                return;
        }

        //
        //  Set a right hand side.
        //
        b[0] = 6.0;
        b[1] = 15.0;
        b[2] = 15.0;
        //
        //  Solve the linear system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        job = 0;
        DGESL.dgesl(a, LDA, n, ipvt, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solution returned by DGESL");
        Console.WriteLine("  (Should be (1,1,1))");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + b[i - 1].ToString().PadLeft(14) + "");
        }

        //
        //  A second right hand side can be solved without refactoring a.
        //
        Console.WriteLine("");
        Console.WriteLine("  Call DGESL for a new right hand");
        Console.WriteLine("  side for the same, factored matrix.");
        //
        //  Set the right hand side.
        //
        b[0] = 1.0;
        b[1] = 4.0;
        b[2] = 7.0;
        //
        //  Solve the system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve a linear system.");

        job = 0;
        DGESL.dgesl(a, LDA, n, ipvt, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solution returned by DGESL");
        Console.WriteLine("  (should be (1,0,0))");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + b[i - 1].ToString().PadLeft(14) + "");
        }

        //
        //  The transposed problem A'*X = B can be solved by DGESL
        //  as well, without any refactoring.
        //
        Console.WriteLine("");
        Console.WriteLine("  Call DGESL for transposed problem.");
        //
        //  Set the right hand side.
        //
        b[0] = 6.0;
        b[1] = 6.0;
        b[2] = -3.0;
        //
        //  Solve the transposed system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Call DGESL to solve a transposed linear system.");

        job = 1;
        DGESL.dgesl(a, LDA, n, ipvt, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solution returned by DGESL");
        Console.WriteLine("  (should be (-1,0,1))");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            Console.WriteLine("  " + b[i - 1].ToString().PadLeft(14) + "");
        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests DGEFA and DGEDI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;
        int LDA = 3;
        //
        //  Matrix listed by columns.
        //
        double[] a =
        {
            1.0, 4.0, 7.0,
            2.0, 5.0, 8.0,
            3.0, 6.0, 0.0
        };
        double[] det = new double[2];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;
        double[] work = new double[N];
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGEFA computes the LU factors;");
        Console.WriteLine("  DGEDI computes the inverse and determinant.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DGEFA.dgefa(ref a, LDA, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("  Error!  The matrix is nearly singular!");
            return;
        }

        //
        //  Get the inverse and determinant.
        //
        Console.WriteLine("");
        Console.WriteLine("  Get the inverse and determinant.");

        job = 11;
        DGEDI.dgedi(ref a, LDA, N, ipvt, ref det, work, job);

        Console.WriteLine("");
        Console.WriteLine("  The determinant = " + det[0] + " * 10^" + det[1] + "");

        Console.WriteLine("");
        Console.WriteLine("  The inverse matrix:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests DGEFA and DGESL.
        //
        //  Discussion:
        //
        //    Solve A*x = b where A is a given matrix, and B a right hand side.
        //
        //    We will also assume that A is stored in the simplest
        //    possible way.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;
        int LDA = N;
        //
        //  The entries of the matrix A are listed by columns.
        //
        double[] a =
        {
            1.0, 4.0, 7.0,
            2.0, 5.0, 8.0,
            3.0, 6.0, 0.0
        };
        double[] b = {6.0, 15.0, 15.0};
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  DGEFA computes the LU factors;");
        Console.WriteLine("  DGESL solves a factored linear system;");
        Console.WriteLine("  tor a general banded matrix.");
        Console.WriteLine("");
        Console.WriteLine("  The number of equations is N = " + N + "");

        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The right hand side B:");
        Console.WriteLine("");
        cout = "";
        for (i = 1; i <= N; i++)
        {
            cout += "  " + b[i - 1].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DGEFA.dgefa(ref a, LDA, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("  DGEFA returned an error flag INFO = " + info + "");
            return;
        }

        //
        //  Solve the system.
        //
        job = 0;
        DGESL.dgesl(a, LDA, N, ipvt, ref b, job);

        Console.WriteLine("");
        Console.WriteLine("  DGESL returns the solution:");
        Console.WriteLine("  (Should be (1,1,1))");
        Console.WriteLine("");

        cout = "";
        for (i = 1; i <= N; i++)
        {
            cout += "  " + b[i - 1].ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests DGEFA and DGESL.
        //
        //  Discussion:
        //
        //    In this example, we solve a relatively large linear system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  For a general banded matrix,");
        Console.WriteLine("  DGEFA computes the LU factors;");
        Console.WriteLine("  DGESL solves a factored linear system;");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to the matrix A and the right hand side B.
        //
        //  The problem is just an enlarged version of the
        //  problem for N = 5, which is:
        //
        //  Matrix A is ( n -1 -1 -1 -1)    Right hand side B is  (1)
        //              (-1  n -1 -1 -1)                          (1)
        //              (-1 -1  n -1 -1)                          (1)
        //              (-1 -1 -1  n -1)                          (1)
        //              (-1 -1 -1 -1  n)                          (1)
        //
        //  Solution is   (1)
        //                (1)
        //                (1)
        //                (1)
        //                (1)
        //
        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 1.0;
        }

        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= N; i++)
            {
                if (i == j)
                {
                    a[i - 1 + (j - 1) * LDA] = N;
                }
                else
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DGEFA.dgefa(ref a, LDA, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("  DGEFA returned an error flag INFO = " + info + "");
            return;
        }

        //
        //  Solve the system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the factored system.");

        job = 0;
        DGESL.dgesl(a, LDA, N, ipvt, ref b, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,1,1,1,1,...,1,1))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }

    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests DGTSL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;

        double[] b = new double[N];
        double[] c = new double[N];
        double[] d = new double[N];
        double[] e = new double[N];
        int i;
        int info;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  For a general tridiagonal matrix,");
        Console.WriteLine("  DGTSL factors and solves a linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        Console.WriteLine("");
        //
        //  Set up the linear system, by storing the values of the
        //  subdiagonal, diagonal, and superdiagonal in C, D, and E,
        //  and the right hand side in B.
        //
        c[0] = 0.0;
        for (i = 2; i <= N; i++)
        {
            c[i - 1] = -1.0;
        }

        for (i = 1; i <= N; i++)
        {
            d[i - 1] = 2.0;
        }

        for (i = 1; i <= N - 1; i++)
        {
            e[i - 1] = -1.0;
        }

        e[N - 1] = 0.0;

        for (i = 1; i <= N - 1; i++)
        {
            b[i - 1] = 0.0;
        }

        b[N - 1] = N + 1;
        //
        // Factor the matrix and solve the system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix and solve the system.");

        info = DGTSL.dgtsl(N, ref c, ref d, ref e, ref b);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  DGTSL returns nonzero INFO = " + info + "");
            return;
        }

        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests DPBCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        int LDA = 2;

        double[] a = new double[LDA * N];
        int j;
        int m;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  For a positive definite symmetric banded matrix,");
        Console.WriteLine("  DPBCO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the number of nonzero diagonals.
        //
        m = 1;
        //
        //  Set the value of the subdiagonal and diagonal.
        //
        for (j = 1; j <= N; j++)
        {
            a[0 + (j - 1) * LDA] = -1.0;
            a[1 + (j - 1) * LDA] = 2.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition.");

        rcond = DPBCO.dpbco(ref a, LDA, N, m, ref z);

        Console.WriteLine("");
        Console.WriteLine("  Reciprocal condition  = " + rcond + "");

    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests DPBDI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 128;
        int LDA = 2;

        double[] a = new double[LDA * N_MAX];
        double[] det = new double[2];
        int info;
        int j;
        int m;
        int n;
        int n_log;

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  For a positive definite symmetric banded matrix,");
        Console.WriteLine("  DPBDI computes the determinant as");
        Console.WriteLine("    det = MANTISSA * 10^EXPONENT");
        Console.WriteLine("");
        Console.WriteLine("  Find the determinant of the -1,2,-1 matrix");
        Console.WriteLine("  for N = 2, 4, 8, 16, 32, 64, 128.");
        Console.WriteLine("");
        Console.WriteLine("  (For this matrix, det ( A ) = N + 1.)");
        //
        //  Set the number of  nonzero diagonals.
        //
        m = 1;

        Console.WriteLine("");
        Console.WriteLine("       N    Mantissa       Exponent");
        Console.WriteLine("");

        n = 1;

        for (n_log = 1; n_log <= 7; n_log++)
        {
            n = 2 * n;

            a[0 + 0 * LDA] = 0.0;
            for (j = 2; j <= n; j++)
            {
                a[0 + (j - 1) * LDA] = -1.0;
            }

            for (j = 1; j <= n; j++)
            {
                a[1 + (j - 1) * LDA] = 2.0;
            }

            info = DPBFA.dpbfa(ref a, LDA, n, m);

            if (info != 0)
            {
                Console.WriteLine("  Error!  DPBFA returns INFO = " + info + "");
                return;
            }

            DPBDI.dpbdi(a, LDA, n, m, ref det);

            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + det[0].ToString().PadLeft(14)
                                   + "  " + det[1].ToString().PadLeft(14) + "");
        }
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests DPBFA and DPBSL.
        //
        //  Discussion:
        //
        //    DPBFA and DPBSL are for a positive definite symmetric band matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;
        int LDA = 2;
        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int info;
        int j;
        int m;

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  For a positive definite symmetric banded matrix,");
        Console.WriteLine("  DPBFA computes the LU factors.");
        Console.WriteLine("  DPBSL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to matrix A and right hand side B.
        //
        //  The problem is just an enlarged version of the
        //  problem for N = 5, which is:
        //
        //  Matrix A is ( 2 -1  0  0  0)    right hand side B is  (1)
        //              (-1  2 -1  0  0)                          (0)
        //              ( 0 -1  2 -1  0)                          (0)
        //              ( 0  0 -1  2 -1)                          (0)
        //              ( 0  0  0 -1  2)                          (1)
        //
        //
        //  solution is   (1)
        //                (1)
        //                (1)
        //                (1)
        //                (1)
        //
        //  Set the right hand side.
        //
        b[0] = 1.0;
        for (i = 2; i <= N - 1; i++)
        {
            b[i - 1] = 0.0;
        }

        b[N - 1] = 1.0;
        //
        //  Set the number of nonzero diagonals.
        //
        m = 1;
        //
        //  Set the value of the subdiagonal and diagonal.
        //
        for (j = 1; j <= N; j++)
        {
            a[0 + (j - 1) * LDA] = -1.0;
            a[1 + (j - 1) * LDA] = 2.0;
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DPBFA.dpbfa(ref a, LDA, N, m);

        if (info != 0)
        {
            Console.WriteLine("  Error!  DPBFA returns INFO = " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        DPBSL.dpbsl(a, LDA, N, m, ref b);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,1,1,1,1,...,1,1))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests DPOCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = N;

        double[] a = new double[LDA * N];
        int i;
        int j;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  For a positive definite symmetric banded matrix,");
        Console.WriteLine("  DPOCO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a[i + j * LDA] = 0.0;
            }
        }

        for (i = 1; i <= N; i++)
        {
            a[i - 1 + (i - 1) * LDA] = 2.0;
            a[i - 1 + (i - 2) * LDA] = i switch
            {
                > 1 => -1.0,
                _ => a[i - 1 + (i - 2) * LDA]
            };

            if (i < N)
            {
                a[i - 1 + i * LDA] = -1.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition.");

        rcond = DPOCO.dpoco(ref a, LDA, N, ref z);

        Console.WriteLine("");
        Console.WriteLine("  Reciprocal condition  = " + rcond + "");

    }

    private static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests DPOFA and DPODI.
        //
        //  Discussion:
        //
        //    DPOFA factors a positive definite symmetric matrix,
        //    and DPODI can compute the determinant or the inverse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] det = new double[2];
        int i;
        int info;
        int j;
        int job;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  For a positive definite symmetric matrix,");
        Console.WriteLine("  DPOFA computes the LU factors.");
        Console.WriteLine("  DPODI computes the inverse or determinant.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a[i + j * LDA] = 0.0;
            }
        }

        for (i = 1; i <= N; i++)
        {
            a[i - 1 + (i - 1) * LDA] = 2.0;
            a[i - 1 + (i - 2) * LDA] = i switch
            {
                > 1 => -1.0,
                _ => a[i - 1 + (i - 2) * LDA]
            };

            if (i < N)
            {
                a[i - 1 + i * LDA] = -1.0;
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DPOFA.dpofa(ref a, LDA, N);

        if (info != 0)
        {
            Console.WriteLine("  Error, DPOFA returns INFO = " + info + "");
            return;
        }

        //
        //  Invert the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Get the determinant and inverse.");

        job = 11;
        DPODI.dpodi(ref a, LDA, N, ref det, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  Determinant = " + det[0] + " * 10^" + det[1] + "");
        //
        //  DPODI produces only the 'upper half triangle' of the inverse,
        //  which is actually symmetric.  Thus, the lower half could be
        //  produced by copying from the upper half.  However, the first row
        //  of A, as returned, is exactly the first row of the inverse.
        //
        Console.WriteLine("");
        Console.WriteLine("  First row of inverse:");
        Console.WriteLine("");
        cout = "";
        for (j = 1; j <= N; j++)
        {
            cout += "  " + a[0 + (j - 1) * LDA].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

    }

    private static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests DPOFA and DPOSL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int info;
        int j;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  For a positive definite symmetric matrix,");
        Console.WriteLine("  DPOFA computes the LU factors.");
        Console.WriteLine("  DPOSL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        for (j = 0; j < N; j++)
        {
            for (i = 0; i < N; i++)
            {
                a[i + j * LDA] = 0.0;
            }
        }

        for (i = 1; i <= N; i++)
        {
            a[i - 1 + (i - 1) * LDA] = 2.0;
            a[i - 1 + (i - 2) * LDA] = i switch
            {
                > 1 => -1.0,
                _ => a[i - 1 + (i - 2) * LDA]
            };

            if (i < N)
            {
                a[i - 1 + i * LDA] = -1.0;
            }
        }

        //
        //  Set the right hand side.
        //
        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
            for (j = 1; j <= N; j++)
            {
                b[i - 1] += a[i - 1 + (j - 1) * LDA] * x[j - 1];
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DPOFA.dpofa(ref a, LDA, N);

        if (info != 0)
        {
            Console.WriteLine("  Error, DPOFA returns INFO = " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        DPOSL.dposl(a, LDA, N, ref b);
        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }

    }

    private static void test19()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests DPPCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] a = new double[N * (N + 1) / 2];
        int i;
        int j;
        int k;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST19");
        Console.WriteLine("  For a positive definite symmetric packed matrix,");
        Console.WriteLine("  DPPCO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                if (i == j - 1)
                {
                    a[k - 1] = -1.0;
                }
                else if (i == j)
                {
                    a[k - 1] = 2.0;
                }
                else
                {
                    a[k - 1] = 0.0;
                }
            }
        }

        //
        //  Estimate the condition.
        //
        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition number.");

        rcond = DPPCO.dppco(ref a, N, ref z);

        Console.WriteLine("");
        Console.WriteLine("  Reciprocal condition number = " + rcond + "");
    }

    private static void test20()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests DPPFA and DPPDI.
        //
        //  Discussion:
        //
        //    DPPFA factors a packed positive definite symmetric matrix,
        //    and DPPDI can compute the determinant or the inverse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        double[] a = new double[N * (N + 1) / 2];
        double[] b = new double[N * N];
        double[] det = new double[2];
        int i;
        int info;
        int j;
        int job;
        int k;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST20");
        Console.WriteLine("  For a positive definite symmetric packed matrix,");
        Console.WriteLine("  DPPFA computes the LU factors.");
        Console.WriteLine("  DPPDI computes the inverse or determinant.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                if (i == j - 1)
                {
                    a[k - 1] = -1.0;
                }
                else if (i == j)
                {
                    a[k - 1] = 2.0;
                }
                else
                {
                    a[k - 1] = 0.0;
                }
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DPPFA.dppfa(ref a, N);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  Error, DPPFA returns INFO = " + info + "");
            return;
        }

        //
        //  Invert the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Get the determinant and inverse.");

        job = 11;
        DPPDI.dppdi(ref a, N, ref det, job);
        //
        //  Print the results.
        //
        Console.WriteLine("");
        Console.WriteLine("  Determinant = " + det[0] + " * 10^" + det[1] + "");
        //
        //  DPPDI produces only the 'upper half triangle' of the inverse,
        //  which is actually symmetric.  Thus, the lower half could be
        //  produced by copying from the upper half.  However, the first row
        //  of A, as returned, is exactly the first row of the inverse.
        //
        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                b[i - 1 + (j - 1) * N] = a[k - 1];
                b[j - 1 + (i - 1) * N] = a[k - 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The inverse matrix:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + b[i - 1 + (j - 1) * N].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test21()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21 tests DPPFA and DPPSL.
        //
        //  Discussion:
        //
        //    DPOFA factors a positive definite symmetric matrix,
        //    and DPOSL can solve a factored linear system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        double[] a = new double[N * (N + 1) / 2];
        double[] b = new double[N];
        int i;
        int info;
        int j;
        int k;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST21");
        Console.WriteLine("  For a positive definite symmetric packed matrix,");
        Console.WriteLine("  DPPFA computes the LU factors.");
        Console.WriteLine("  DPPSL solves a factored linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
        }

        //
        //  Set the matrix A.
        //
        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                if (i == j - 1)
                {
                    a[k - 1] = -1.0;
                    b[i - 1] += a[k - 1] * x[j - 1];
                    b[j - 1] += a[k - 1] * x[i - 1];
                }
                else if (i == j)
                {
                    a[k - 1] = 2.0;
                    b[i - 1] += a[k - 1] * x[i - 1];
                }
                else
                {
                    a[k - 1] = 0.0;
                }
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DPPFA.dppfa(ref a, N);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  Error, DPPFA returns INFO = " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        DPPSL.dppsl(a, N, ref b);
        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void test22()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests DPTSL.
        //
        //  Discussion:
        //
        //    DPTSL factors and solves a positive definite symmetric tridiagonal system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        double[] b = new double[N];
        double[] d = new double[N];
        double[] e = new double[N];
        int i;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST22");
        Console.WriteLine("  For a positive definite symmetric tridiagonal matrix,");
        Console.WriteLine("  DPTSL factors and solves a linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Set the matrix A.
        //
        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
        }

        for (i = 1; i <= N; i++)
        {
            d[i - 1] = 2.0;
        }

        for (i = 1; i <= N - 1; i++)
        {
            e[i - 1] = -1.0;
        }

        e[N - 1] = 0.0;

        for (i = 1; i <= N; i++)
        {
            switch (i)
            {
                case > 1:
                    b[i - 1] += e[i - 2] * x[i - 2];
                    break;
            }

            b[i - 1] += d[i - 1] * x[i - 1];
            if (i < N)
            {
                b[i - 1] += e[i - 1] * x[i];
            }
        }

        //
        //  Factor and solve the system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix and solve the system.");

        DPTSL.dptsl(N, ref d, e, ref b);
        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }
    }

    private static void dqrdc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DQRDC_TEST tests DQRDC.
        //
        //  Discussion:
        //
        //    DQRDC computes the QR factorization.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;
        int P = 3;
        int LDA = N;

        double[] a =
        {
            1.0, 1.0, 0.0,
            1.0, 0.0, 1.0,
            0.0, 1.0, 1.0
        };
        double[] b = new double[LDA * P];
        int i;
        int info;
        int[] ipvt = new int[P];
        int j;
        int job;
        int k;
        double[] q = new double[N * N];
        double[] qraux = new double[P];
        double[] qty = new double[N];
        double[] qy = new double[N];
        double[] r = new double[N * P];
        double[] rsd = new double[N];
        double[] work = new double[P];
        double[] xb = new double[N];
        double[] y = new double[N];
        string cout;

        Console.WriteLine("");
        Console.WriteLine("DQRDC_TEST");
        Console.WriteLine("  DQRDC computes the QR decomposition of a general");
        Console.WriteLine("  matrix, but does not return Q and R explicitly.");
        Console.WriteLine("");
        Console.WriteLine("  Show how Q and R can be recovered using DQRSL.");
        //
        //  Print the matrix A.
        //
        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= P; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Decompose the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Decompose the matrix.");

        job = 0;
        for (j = 1; j <= P; j++)
        {
            ipvt[j - 1] = 0;
        }

        DQRDC.dqrdc(ref a, LDA, N, P, ref qraux, ref ipvt, work, job);
        //
        //  Print out what DQRDC has stored in A...
        //
        Console.WriteLine("");
        Console.WriteLine("  The packed matrix A which describes Q and R:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= P; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  ...and in QRAUX.
        //
        Console.WriteLine("");
        Console.WriteLine("  The QRAUX vector, containing some additional");
        Console.WriteLine("  information defining Q:");
        Console.WriteLine("");

        cout = "";
        for (i = 1; i <= N; i++)
        {
            cout += "  " + qraux[i - 1].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        //
        //  Print out the resulting R factor.
        //
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= P; j++)
            {
                if (j < i)
                {
                    r[i - 1 + (j - 1) * N] = 0.0;
                }
                else
                {
                    r[i - 1 + (j - 1) * N] = a[i - 1 + (j - 1) * LDA];
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The R factor:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= P; j++)
            {
                cout += "  " + r[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Call DQRSL to extract the information about the Q matrix.
        //  We do this, essentially, by asking DQRSL to tell us the
        //  value of Q*Y, where Y is a column of the identity matrix.
        //
        job = 10000;

        for (i = 1; i <= N; i++)
        {
            //
            //  Set the vector Y.
            //
            for (j = 1; j <= N; j++)
            {
                y[j - 1] = 0.0;
            }

            y[i - 1] = 1.0;
            //
            //  Ask DQRSL to tell us what Q*Y is.
            //
            info = DQRSL.dqrsl(a, LDA, N, P, qraux, y, ref qy, ref qty, ref b, ref rsd, ref xb, job);

            if (info != 0)
            {
                Console.WriteLine("  Error!  DQRSL returns INFO = " + info + "");
                return;
            }

            //
            //  Copy QY into the appropriate column of Q.
            //
            for (j = 1; j <= N; j++)
            {
                q[j - 1 + (i - 1) * N] = qy[j - 1];
            }
        }

        //
        //  Now print out the Q matrix we have extracted.
        //
        Console.WriteLine("");
        Console.WriteLine("  The Q factor:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + q[i - 1 + (j - 1) * N].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Compute Q*R to verify that it equals A.
        //
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= P; j++)
            {
                b[i - 1 + (j - 1) * LDA] = 0.0;
                for (k = 1; k <= N; k++)
                {
                    b[i - 1 + (j - 1) * LDA] += q[i - 1 + (k - 1) * N] * r[k - 1 + (j - 1) * N];
                }
            }
        }

        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The product Q * R:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= P; j++)
            {
                cout += "  " + b[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    private static void dqrsl_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DQRSL_TEST tests DQRSL.
        //
        //  Discussion:
        //
        //    DQRSL can solve a linear system that was factored by DQRDC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int P = 3;
        int LDA = N;

        double[] a =
        {
            1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 2.0, 3.0, 4.0, 5.0,
            1.0, 4.0, 9.0, 16.0, 25.0
        };
        double[] b = new double[P];
        double[] b2 =
        {
            -3.02, 4.4914286, -0.72857143
        };
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int job;
        double[] qraux = new double[P];
        double[] qty = new double[N];
        double[] qy = new double[N];
        double[] r = new double[N];
        double[] work = new double[P];
        double[] xb = new double[N];
        double[] y =
        {
            1.0,
            2.3,
            4.6,
            3.1,
            1.2
        };
        string cout;

        Console.WriteLine("");
        Console.WriteLine("DQRSL_TEST");
        Console.WriteLine("  DQRSL solves a rectangular linear system A*x=b in the");
        Console.WriteLine("  least squares sense after A has been factored by DQRDC.");

        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            cout = "";
            for (j = 0; j < P; j++)
            {
                cout += "  " + a[i + j * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Decompose the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine(" Decompose the matrix.");

        job = 0;

        for (j = 0; j < P; j++)
        {
            ipvt[j] = 0;
        }

        DQRDC.dqrdc(ref a, LDA, N, P, ref qraux, ref ipvt, work, job);
        //
        //  Call DQRSL to compute the least squares solution A*x=b.
        //
        job = 110;

        info = DQRSL.dqrsl(a, LDA, N, P, qraux, y, ref qy, ref qty, ref b, ref r, ref xb, job);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DQRSL_TEST - Warning!");
            Console.WriteLine("  DQRSL returns INFO = " + info + "");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("      X          X(expected):");
        Console.WriteLine("");

        for (i = 0; i < P; i++)
        {
            Console.WriteLine("  " + b[i].ToString().PadLeft(14)
                                   + "  " + b2[i].ToString().PadLeft(14) + "");
        }

    }

    private static void test24()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST24 tests DSICO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;
        int LDA = N;

        double[] a = new double[LDA * N];
        int i;
        int[] ipvt = new int[N];
        int j;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST24");
        Console.WriteLine("  For a symmetric indefinite matrix,");
        Console.WriteLine("  DSICO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to the matrix A and the right hand side B.
        //
        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                if (i == j)
                {
                    a[i - 1 + (j - 1) * LDA] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else
                {
                    a[i - 1 + (j - 1) * LDA] = 0.0;
                }
            }
        }

        //
        //  Estimate the condition.
        //
        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition.");

        rcond = DSICO.dsico(ref a, LDA, N, ref ipvt, ref z);

        Console.WriteLine("");
        Console.WriteLine("  Estimated reciprocal condition = " + rcond + "");

    }

    private static void test25()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST25 tests DSIFA and DSISL.
        //
        //  Discussion:
        //
        //    DSIFA and DSISL are for symmetric indefinite matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST25");
        Console.WriteLine("  DSIFA factor a symmetric indefinite matrix;");
        Console.WriteLine("  DSISL solves a factored linear system,");
        Console.WriteLine("  for a symmetric indefinite matrix.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to the matrix A and the right hand side B.
        //
        for (i = 1; i < N; i++)
        {
            b[i - 1] = 0.0;
        }

        b[N - 1] = N + 1;

        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= N; j++)
            {
                if (i == j)
                {
                    a[i - 1 + (j - 1) * LDA] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[i - 1 + (j - 1) * LDA] = -1.0;
                }
                else
                {
                    a[i - 1 + (j - 1) * LDA] = 0.0;
                }
            }
        }

        //
        //  Factor the matrix A.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DSIFA.dsifa(ref a, LDA, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("  Error!  DSIFA returns INFO = " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        DSISL.dsisl(a, LDA, N, ipvt, ref b);
        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }

    }

    private static void test26()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST26 tests DSPCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;

        double[] a = new double[N * (N + 1) / 2];
        int i;
        int[] ipvt = new int[N];
        int j;
        int k;
        double rcond;
        double[] z = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST26");
        Console.WriteLine("  For a symmetric indefinite packed matrix,");
        Console.WriteLine("  DSPCO estimates the reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to the matrix A.
        //
        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                if (i == j)
                {
                    a[k - 1] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[k - 1] = -1.0;
                }
                else
                {
                    a[k - 1] = 0.0;
                }
            }
        }

        //
        //  Estimate the condition.
        //
        Console.WriteLine("");
        Console.WriteLine("  Estimate the condition.");

        rcond = DSPCO.dspco(ref a, N, ref ipvt, ref z);

        Console.WriteLine("");
        Console.WriteLine("  Estimated reciprocal condition = " + rcond + "");
    }

    private static void test27()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST27 tests DSPFA and DSPSL.
        //
        //  Discussion:
        //
        //    DSPFA and DSPSL are for packed symmetric indefinite matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 100;

        double[] a = new double[N * (N + 1) / 2];
        double[] b = new double[N];
        int i;
        int info;
        int[] ipvt = new int[N];
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST27");
        Console.WriteLine("  DSPFA computes the LU factors,");
        Console.WriteLine("  DSPSL solves a factored linear system,");
        Console.WriteLine("  for a symmetric indefinite packed matrix,");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Assign values to the matrix A and the right hand side B.
        //
        for (i = 1; i <= N - 1; i++)
        {
            b[i - 1] = 0.0;
        }

        b[N - 1] = N + 1;

        k = 0;
        for (j = 1; j <= N; j++)
        {
            for (i = 1; i <= j; i++)
            {
                k += 1;
                if (i == j)
                {
                    a[k - 1] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[k - 1] = -1.0;
                }
                else
                {
                    a[k - 1] = 0.0;
                }
            }
        }

        //
        //  Factor the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Factor the matrix.");

        info = DSPFA.dspfa(ref a, N, ref ipvt);

        if (info != 0)
        {
            Console.WriteLine("  Error!  DSPFA returns INFO = " + info + "");
            return;
        }

        //
        //  Solve the linear system.
        //
        Console.WriteLine("");
        Console.WriteLine("  Solve the linear system.");

        DSPSL.dspsl(a, N, ipvt, ref b);
        //
        //  Print the result.
        //
        Console.WriteLine("");
        Console.WriteLine("  The first and last 5 entries of solution:");
        Console.WriteLine("  (Should be (1,2,3,4,5,...,n-1,n))");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            if (i <= 5 || N - 5 < i)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                       + "  " + b[i - 1].ToString().PadLeft(14) + "");
            }

            switch (i)
            {
                case 5:
                    Console.WriteLine("  ......  ..............");
                    break;
            }
        }

    }

    private static void dsvdc_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSVDC_TEST tests DSVDC.
        //
        //  Discussion:
        //
        //    DSVDC computes the singular value decomposition:
        //
        //      A = U * S * V'
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 6;
        int N = 4;

        double[] a = new double[M * N];
        double[] b = new double[M * N];
        //
        //  E must be dimensioned at least maximum(M+1,N).
        //
        double[] e = new double[M + N];
        int i;
        int info;
        int j;
        int job;
        int k;
        int lda;
        int ldu;
        int ldv;
        //
        //  S must be dimensioned at least maximum(M+1,N).
        //
        double[] s = new double[M + N];
        int seed;
        double[] sigma = new double[M * N];
        double[] u = new double[M * M];
        double[] v = new double[N * N];
        double[] work = new double[M];

        string cout;

        Console.WriteLine("");
        Console.WriteLine("DSVDC_TEST");
        Console.WriteLine("  For an MxN matrix A in general storage,");
        Console.WriteLine("  DSVDC computes the singular value decomposition:");
        Console.WriteLine("    A = U * S * V'");
        Console.WriteLine("");
        Console.WriteLine("  Matrix rows M =    " + M + "");
        Console.WriteLine("  Matrix columns N = " + N + "");
        //
        //  Set A.
        //
        seed = 123456789;

        a = UniformRNG.r8mat_uniform_01(M, N, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  The matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= M; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * M].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        //
        //  Decompose the matrix.
        //
        Console.WriteLine("");
        Console.WriteLine("  Decompose the matrix.");

        job = 11;
        lda = M;
        ldu = M;
        ldv = N;

        info = DSVDC.dsvdc(ref a, lda, M, N, ref s, ref e, ref u, ldu, ref v, ldv, work, job);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DSVDC_TEST - Warning:");
            Console.WriteLine("  DSVDC returned nonzero INFO = " + info + "");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Singular values:");
        Console.WriteLine("");

        for (i = 1; i <= Math.Min(M, N); i++)
        {
            Console.WriteLine("  "
                              + (i + 1).ToString().PadLeft(4) + "  "
                              + s[i - 1].ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Left Singular Vector Matrix U:");
        Console.WriteLine("");

        for (i = 1; i <= M; i++)
        {
            cout = "";
            for (j = 1; j <= M; j++)
            {
                cout += "  " + u[i - 1 + (j - 1) * M].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Right Singular Vector Matrix V:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + v[i - 1 + (j - 1) * N].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

        for (i = 1; i <= M; i++)
        {
            for (j = 1; j <= N; j++)
            {
                if (i == j)
                {
                    sigma[i - 1 + (j - 1) * M] = s[i - 1];
                }
                else
                {
                    sigma[i - 1 + (j - 1) * M] = 0.0;
                }
            }
        }

        //
        //  Verify that A = U * S * V'.
        //
        for (i = 1; i <= M; i++)
        {
            for (k = 1; k <= N; k++)
            {
                b[i - 1 + (k - 1) * M] = 0.0;
                for (j = 1; j <= N; j++)
                {
                    b[i - 1 + (k - 1) * M] += sigma[i - 1 + (j - 1) * M] * v[k - 1 + (j - 1) * N];
                }
            }
        }

        for (i = 1; i <= M; i++)
        {
            for (k = 1; k <= N; k++)
            {
                a[i - 1 + (k - 1) * M] = 0.0;
                for (j = 1; j <= M; j++)
                {
                    a[i - 1 + (k - 1) * M] += u[i - 1 + (j - 1) * M] * b[j - 1 + (k - 1) * M];
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  The product U * S * V' (should equal A):");
        Console.WriteLine("");

        for (i = 1; i <= M; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * M].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

    }

    private static void test29()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST29 tests DTRCO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = N;

        double[] a = new double[LDA * N];
        int i;
        int j;
        int job;
        double rcond;
        int seed = 123456789;
        double[] z = new double[N];
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST29");
        Console.WriteLine("  For a triangular matrix,");
        Console.WriteLine("  DTRCO computes the LU factors and");
        Console.WriteLine("  computes its reciprocal condition number.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Lower triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(LDA, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = i + 1; j <= N; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Lower triangular matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        job = 0;
        rcond = DTRCO.dtrco(a, LDA, N, ref z, job);

        Console.WriteLine("");
        Console.WriteLine("  The reciprocal condition number = " + rcond + "");
        //
        //  Upper triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(LDA, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Upper triangular matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        job = 1;

        rcond = DTRCO.dtrco(a, LDA, N, ref z, job);

        Console.WriteLine("");
        Console.WriteLine("  The reciprocal condition number = " + rcond + "");

    }

    private static void test30()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST30 tests DTRDI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = N;

        double[] a = new double[LDA * N];
        double[] det = new double[2];
        int i;
        int j;
        int job;
        int seed = 123456789;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST30");
        Console.WriteLine("  For a triangular matrix,");
        Console.WriteLine("  DTRDI computes the determinant or inverse.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Lower triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(N, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = i + 1; j <= N; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Lower triangular matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        job = 110;

        DTRDI.dtrdi(ref a, LDA, N, ref det, job);

        Console.WriteLine("");
        Console.WriteLine("  The determinant = " + det[0] + " * 10^(" + det[1] + ").");

        Console.WriteLine("");
        Console.WriteLine("  The inverse matrix:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        //
        //  Upper triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(N, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Upper triangular matrix A:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        job = 111;

        DTRDI.dtrdi(ref a, LDA, N, ref det, job);

        Console.WriteLine("");
        Console.WriteLine("  The determinant = " + det[0] + " * 10^(" + det[1] + ").");

        Console.WriteLine("");
        Console.WriteLine("  The inverse matrix:");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            cout = "";
            for (j = 1; j <= N; j++)
            {
                cout += "  " + a[i - 1 + (j - 1) * LDA].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    private static void test31()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST31 tests DTRSL.
        //
        //  Discussion:
        //
        //    DTRSL solves triangular linear systems.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int LDA = 5;

        double[] a = new double[LDA * N];
        double[] b = new double[N];
        int i;
        int j;
        int job;
        int seed = 123456789;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST31");
        Console.WriteLine("  For a triangular matrix,");
        Console.WriteLine("  DTRSL solves a linear system.");
        Console.WriteLine("  The matrix size is N = " + N + "");
        //
        //  Lower triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(N, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = i + 1; j <= N; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
            for (j = 1; j <= N; j++)
            {
                b[i - 1] += a[i - 1 + (j - 1) * LDA] * x[j - 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  For a lower triangular matrix A,");
        Console.WriteLine("  solve A * x = b");

        job = 00;

        DTRSL.dtrsl(a, LDA, N, ref b, job);

        Console.WriteLine("");
        Console.WriteLine("  The solution (should be 1,2,3,4,5):");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + b[i - 1].ToString().PadLeft(14) + "");
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
            for (j = 1; j <= N; j++)
            {
                b[i - 1] += a[j - 1 + (i - 1) * LDA] * x[j - 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  For a lower triangular matrix A,");
        Console.WriteLine("  solve A' * x = b");

        job = 10;

        DTRSL.dtrsl(a, LDA, N, ref b, job);

        Console.WriteLine("");
        Console.WriteLine("  The solution (should be 1,2,3,4,5):");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + b[i - 1].ToString().PadLeft(14) + "");
        }

        //
        //  Upper triangular matrix A.
        //
        a = UniformRNG.r8mat_uniform_01(N, N, ref seed);

        for (i = 1; i <= N; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * LDA] = 0.0;
            }
        }

        for (i = 1; i <= N; i++)
        {
            x[i - 1] = i;
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
            for (j = 1; j <= N; j++)
            {
                b[i - 1] += a[i - 1 + (j - 1) * LDA] * x[j - 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  For an upper triangular matrix A,");
        Console.WriteLine("  solve A * x = b");

        job = 01;

        DTRSL.dtrsl(a, LDA, N, ref b, job);

        Console.WriteLine("");
        Console.WriteLine("  The solution (should be 1,2,3,4,5):");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + b[i - 1].ToString().PadLeft(14) + "");
        }

        for (i = 1; i <= N; i++)
        {
            b[i - 1] = 0.0;
            for (j = 1; j <= N; j++)
            {
                b[i - 1] += a[j - 1 + (i - 1) * LDA] * x[j - 1];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  For an upper triangular matrix A,");
        Console.WriteLine("  solve A' * x = b");

        job = 11;

        DTRSL.dtrsl(a, LDA, N, ref b, job);

        Console.WriteLine("");
        Console.WriteLine("  The solution (should be 1,2,3,4,5):");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(6) + "  "
                              + b[i - 1].ToString().PadLeft(14) + "");
        }

    }
}