using System;
using Burkardt.Sequence;
using Burkardt.Types;

namespace JacobiEigenValueTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for JACOBI_EIGENVALUE_TEST.
        //
        //  Discussion:
        //
        //    JACOBI_EIGENVALUE_TEST tests the JACOBI_EIGENVALUE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("JACOBI_EIGENVALUE_TEST");
        Console.WriteLine("  Test the JACOBI_EIGENVALUE library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("JACOBI_EIGENVALUE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses a 4x4 test matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 4;

        double[] a =  {
                4.0, -30.0, 60.0, -35.0,
                -30.0, 300.0, -675.0, 420.0,
                60.0, -675.0, 1620.0, -1050.0,
                -35.0, 420.0, -1050.0, 700.0
            }
            ;
        double[] d = new double[N];
        int it_num = 0;
        int rot_num = 0;
        double[] v = new double[N * N];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  For a symmetric matrix A,");
        Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
        Console.WriteLine("  and eigenvectors V so that A * V = D * V.");

        typeMethods.r8mat_print(N, N, a, "  Input matrix A:");

        const int it_max = 100;

        Jacobi.jacobi_eigenvalue(N, a, it_max, ref v, ref d, ref it_num, ref rot_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + it_num + "");
        Console.WriteLine("  Number of rotations  = " + rot_num + "");

        typeMethods.r8vec_print(N, d, "  Eigenvalues D:");

        typeMethods.r8mat_print(N, N, v, "  Eigenvector matrix V:");
        //
        //  Compute eigentest.
        //
        double error_frobenius = typeMethods.r8mat_is_eigen_right(N, N, a, v, d);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                          + error_frobenius + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses a 4x4 test matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 4;

        double[] a =  {
                4.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 3.0, 0.0,
                0.0, 0.0, 0.0, 2.0
            }
            ;
        double[] d = new double[N];
        int it_num = 0;
        int rot_num = 0;
        double[] v = new double[N * N];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For a symmetric matrix A,");
        Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
        Console.WriteLine("  and eigenvectors V so that A * V = D * V.");
        Console.WriteLine("");
        Console.WriteLine("As a sanity check, input a diagonal matrix.");

        typeMethods.r8mat_print(N, N, a, "  Input matrix A:");

        const int it_max = 100;

        Jacobi.jacobi_eigenvalue(N, a, it_max, ref v, ref d, ref it_num, ref rot_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + it_num + "");
        Console.WriteLine("  Number of rotations  = " + rot_num + "");

        typeMethods.r8vec_print(N, d, "  Eigenvalues D:");

        typeMethods.r8mat_print(N, N, v, "  Eigenvector matrix V:");
        //
        //  Compute eigentest.
        //
        double error_frobenius = typeMethods.r8mat_is_eigen_right(N, N, a, v, d);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                          + error_frobenius + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses a 5x5 test matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        double[] a = new double[N * N];
        double[] d = new double[N];
        int it_num = 0;
        int j;
        int rot_num = 0;
        double[] v = new double[N * N];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For a symmetric matrix A,");
        Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
        Console.WriteLine("  and eigenvectors V so that A * V = D * V.");
        Console.WriteLine("");
        Console.WriteLine("  Use the discretized second derivative matrix.");

        for (j = 0; j < N; j++)
        {
            int i;
            for (i = 0; i < N; i++)
            {
                if (i == j)
                {
                    a[i + j * N] = -2.0;
                }
                else if (i == j + 1 || i == j - 1)
                {
                    a[i + j * N] = 1.0;
                }
                else
                {
                    a[i + j * N] = 0.0;
                }
            }
        }

        typeMethods.r8mat_print(N, N, a, "  Input matrix A:");

        const int it_max = 100;

        Jacobi.jacobi_eigenvalue(N, a, it_max, ref v, ref d, ref it_num, ref rot_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations = " + it_num + "");
        Console.WriteLine("  Number of rotations  = " + rot_num + "");

        typeMethods.r8vec_print(N, d, "  Eigenvalues D:");

        typeMethods.r8mat_print(N, N, v, "  Eigenvector matrix V:");
        //
        //  Compute eigentest.
        //
        double error_frobenius = typeMethods.r8mat_is_eigen_right(N, N, a, v, d);
        Console.WriteLine("");
        Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                          + error_frobenius + "");
    }
}