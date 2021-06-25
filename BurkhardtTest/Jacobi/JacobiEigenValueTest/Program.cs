using System;
using Burkardt;
using Burkardt.Types;

namespace JacobiEigenValueTest
{
    class Program
    {
        static void Main(string[] args)
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

        static void test01()

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
            int N = 4;

            double[] a =  {
                4.0, -30.0, 60.0, -35.0,
                -30.0, 300.0, -675.0, 420.0,
                60.0, -675.0, 1620.0, -1050.0,
                -35.0, 420.0, -1050.0, 700.0
            }
            ;
            double[] d = new double[N];
            double error_frobenius;
            int it_max;
            int it_num = 0;
            int n = N;
            int rot_num = 0;
            double[] v = new double[N * N];

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  For a symmetric matrix A,");
            Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
            Console.WriteLine("  and eigenvectors V so that A * V = D * V.");

            typeMethods.r8mat_print(n, n, a, "  Input matrix A:");

            it_max = 100;

            Jacobi.jacobi_eigenvalue(n, a, it_max, ref v, ref d, ref it_num, ref rot_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of iterations = " + it_num + "");
            Console.WriteLine("  Number of rotations  = " + rot_num + "");

            typeMethods.r8vec_print(n, d, "  Eigenvalues D:");

            typeMethods.r8mat_print(n, n, v, "  Eigenvector matrix V:");
            //
            //  Compute eigentest.
            //
            error_frobenius = typeMethods.r8mat_is_eigen_right(n, n, a, v, d);
            Console.WriteLine("");
            Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                 + error_frobenius + "");
        }

        static void test02()

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
            int N = 4;

            double[] a =  {
                4.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 3.0, 0.0,
                0.0, 0.0, 0.0, 2.0
            }
            ;
            double[] d = new double[N];
            double error_frobenius;
            int it_max;
            int it_num = 0;
            int n = N;
            int rot_num = 0;
            double[] v = new double[N * N];

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  For a symmetric matrix A,");
            Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
            Console.WriteLine("  and eigenvectors V so that A * V = D * V.");
            Console.WriteLine("");
            Console.WriteLine("As a sanity check, input a diagonal matrix.");

            typeMethods.r8mat_print(n, n, a, "  Input matrix A:");

            it_max = 100;

            Jacobi.jacobi_eigenvalue(n, a, it_max, ref v, ref d, ref it_num, ref rot_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of iterations = " + it_num + "");
            Console.WriteLine("  Number of rotations  = " + rot_num + "");

            typeMethods.r8vec_print(n, d, "  Eigenvalues D:");

            typeMethods.r8mat_print(n, n, v, "  Eigenvector matrix V:");
            //
            //  Compute eigentest.
            //
            error_frobenius = typeMethods.r8mat_is_eigen_right(n, n, a, v, d);
            Console.WriteLine("");
            Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                 + error_frobenius + "");
        }

        static void test03()

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
            int N = 5;

            double[] a = new double[N * N];
            double[] d = new double[N];
            double error_frobenius;
            int i;
            int it_max;
            int it_num = 0;
            int j;
            int n = N;
            int rot_num = 0;
            double[] v = new double[N * N];

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  For a symmetric matrix A,");
            Console.WriteLine("  JACOBI_EIGENVALUE computes the eigenvalues D");
            Console.WriteLine("  and eigenvectors V so that A * V = D * V.");
            Console.WriteLine("");
            Console.WriteLine("  Use the discretized second derivative matrix.");

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        a[i + j * n] = -2.0;
                    }
                    else if (i == j + 1 || i == j - 1)
                    {
                        a[i + j * n] = 1.0;
                    }
                    else
                    {
                        a[i + j * n] = 0.0;
                    }
                }
            }

            typeMethods.r8mat_print(n, n, a, "  Input matrix A:");

            it_max = 100;

            Jacobi.jacobi_eigenvalue(n, a, it_max, ref v, ref d, ref it_num, ref rot_num);

            Console.WriteLine("");
            Console.WriteLine("  Number of iterations = " + it_num + "");
            Console.WriteLine("  Number of rotations  = " + rot_num + "");

            typeMethods.r8vec_print(n, d, "  Eigenvalues D:");

            typeMethods.r8mat_print(n, n, v, "  Eigenvector matrix V:");
            //
            //  Compute eigentest.
            //
            error_frobenius = typeMethods.r8mat_is_eigen_right(n, n, a, v, d);
            Console.WriteLine("");
            Console.WriteLine("  Frobenius norm error in eigensystem A*V-D*V = "
                 + error_frobenius + "");
        }
    }
}