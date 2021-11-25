using System;
using Burkardt.CholeskyNS;
using Burkardt.Error;
using Burkardt.Laplacian;
using Burkardt.Types;

namespace LaplacianTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAPLACIAN_TEST.
        //
        //  Discussion:
        //
        //    LAPLACIAN_TEST tests the LAPLACIAN library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2013
        //
        //  Author:
        //
        //   John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LAPLACIAN_TEST");
        Console.WriteLine("  Test the LAPLACIAN library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("LAPLACIAN_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests L1DD and similar routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  A full-storage matrix is returned by:");
        Console.WriteLine("  L1DD: Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN: Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND: Neumann/Dirichlet BC;");
        Console.WriteLine("  L1NN: Neumann/Neumann BC;");
        Console.WriteLine("  L1PP: Periodic BC;");

        for (test = 1; test <= 2; test++)
        {
            double h = test switch
            {
                1 => 1.0,
                _ => 1.0 / (n + 1)
            };

            Console.WriteLine("");
            Console.WriteLine("  Using spacing H = " + h + "");

            double[] l = L1DD.l1dd(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1DD:");

            l = L1DN.l1dn(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1DN:");

            l = L1ND.l1nd(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1ND:");

            l = L1NN.l1nn(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1NN:");

            l = L1PP.l1pp(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1PP:");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests L1DD_APPLY and similar functions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 9;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  The Laplacian L is applied to data U by:");
        Console.WriteLine("  L1DD_APPLY for Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN_APPLY for Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND_APPLY for Neumann/Dirichlet BC;");
        Console.WriteLine("  L1NN_APPLY for Neumann/Neumann BC;");
        Console.WriteLine("  L1PP_APPLY for Periodic BC;");

        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (i + 1) / (double) (n + 1);
        }

        const double h = 1.0 / (n + 1);

        Console.WriteLine("");
        Console.WriteLine(" Using spacing H = " + h + "");

        double[] u = new double[n];
        for (i = 0; i < n; i++)
        {
            u[i] = x[i] * (1.0 - x[i]);
        }

        typeMethods.r8vec_print(n, u, "  Vector U:");

        double[] lu = L1DD.l1dd_apply(n, h, u);
        typeMethods.r8vec_print(n, lu, "  L1DD(U):");

        lu = L1DN.l1dn_apply(n, h, u);
        typeMethods.r8vec_print(n, lu, "  L1DN(U):");

        lu = L1ND.l1nd_apply(n, h, u);
        typeMethods.r8vec_print(n, lu, "  L1ND(U):");

        lu = L1NN.l1nn_apply(n, h, u);
        typeMethods.r8vec_print(n, lu, "  L1NN(U):");

        lu = L1PP.l1pp_apply(n, h, u);
        typeMethods.r8vec_print(n, lu, "  L1PP(U):");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests L1DD_EIGEN and similar routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Compute eigen information for the Laplacian:");
        Console.WriteLine("  L1DD_EIGEN for Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN_EIGEN for Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND_EIGEN for Neumann/Dirichlet BC;");
        Console.WriteLine("  L1NN_EIGEN for Neumann/Neumann BC;");
        Console.WriteLine("  L1PP_EIGEN for Periodic BC;");

        double[] v = new double[n * n];
        double[] lambda = new double[n];

        for (test = 1; test <= 2; test++)
        {
            double h = test switch
            {
                1 => 1.0,
                _ => 1.0 / (n + 1)
            };

            Console.WriteLine("");
            Console.WriteLine("  Using spacing H = " + h + "");

            double[] a = L1DD.l1dd(n, h);
            L1DD.l1dd_eigen(n, h, ref v, ref lambda);
            typeMethods.r8vec_print(n, lambda, "  L1DD Eigenvalues:");
            typeMethods.r8mat_print(n, n, v, "  L1DD Eigenvectors:");
            double err = Eigen.eigen_error(n, n, a, v, lambda);
            Console.WriteLine("");
            Console.WriteLine("  L1DD eigenerror = " + err + "");

            a = L1DN.l1dn(n, h);
            L1DN.l1dn_eigen(n, h, ref v, ref lambda);
            typeMethods.r8vec_print(n, lambda, "  L1DN Eigenvalues:");
            typeMethods.r8mat_print(n, n, v, "  L1DN Eigenvectors:");
            err = Eigen.eigen_error(n, n, a, v, lambda);
            Console.WriteLine("");
            Console.WriteLine("  L1DN eigenerror = " + err + "");

            a = L1ND.l1nd(n, h);
            L1ND.l1nd_eigen(n, h, ref v, ref lambda);
            typeMethods.r8vec_print(n, lambda, "  L1ND Eigenvalues:");
            typeMethods.r8mat_print(n, n, v, "  L1ND Eigenvectors:");
            err = Eigen.eigen_error(n, n, a, v, lambda);
            Console.WriteLine("");
            Console.WriteLine("  L1ND eigenerror = " + err + "");

            a = L1NN.l1nn(n, h);
            L1NN.l1nn_eigen(n, h, ref v, ref lambda);
            typeMethods.r8vec_print(n, lambda, "  L1NN Eigenvalues:");
            typeMethods.r8mat_print(n, n, v, "  L1NN Eigenvectors:");
            err = Eigen.eigen_error(n, n, a, v, lambda);
            Console.WriteLine("");
            Console.WriteLine("  L1NN eigenerror = " + err + "");

            a = L1PP.l1pp(n, h);
            L1PP.l1pp_eigen(n, h, ref v, ref lambda);
            typeMethods.r8vec_print(n, lambda, "  L1PP Eigenvalues:");
            typeMethods.r8mat_print(n, n, v, "  L1PP Eigenvectors:");
            err = Eigen.eigen_error(n, n, a, v, lambda);
            Console.WriteLine("");
            Console.WriteLine("  L1PP eigenerror = " + err + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests L1DD_INVERSE and similar routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  The inverse of a full-storage matrix is returned by:");
        Console.WriteLine("  L1DD_INVERSE: Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN_INVERSE: Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND_INVERSE: Neumann/Dirichlet BC;");

        for (test = 1; test <= 2; test++)
        {
            double h = test switch
            {
                1 => 1.0,
                _ => 1.0 / (n + 1)
            };

            Console.WriteLine("");
            Console.WriteLine("  Using spacing H = " + h + "");

            double[] l = L1DD.l1dd(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1DD:");
            double[] linv = L1DD.l1dd_inverse(n, h);
            typeMethods.r8mat_print(n, n, linv, "  L1DD_INVERSE:");
            double err = Inverse.inverse_error(n, l, linv);
            Console.WriteLine("");
            Console.WriteLine("  L1DD inverse error = " + err + "");

            l = L1DN.l1dn(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1DN:");
            linv = L1DN.l1dn_inverse(n, h);
            typeMethods.r8mat_print(n, n, linv, "  L1DN_INVERSE:");
            err = Inverse.inverse_error(n, l, linv);
            Console.WriteLine("");
            Console.WriteLine("  L1DN inverse error = " + err + "");

            l = L1ND.l1nd(n, h);
            typeMethods.r8mat_print(n, n, l, "  L1ND:");
            linv = L1ND.l1nd_inverse(n, h);
            typeMethods.r8mat_print(n, n, linv, "  L1ND_INVERSE:");
            err = Inverse.inverse_error(n, l, linv);
            Console.WriteLine("");
            Console.WriteLine("  L1ND inverse error = " + err + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests L1DD_CHOLESKY and similar routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Compute upper Cholesky factors for the Laplacian:");
        Console.WriteLine("  L1DD_CHOLESKY for Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN_CHOLESKY for Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND_CHOLESKY for Neumann/Dirichlet BC;");
        Console.WriteLine("  L1NN_CHOLESKY for Neumann/Neumann BC;");
        Console.WriteLine("  L1PP_CHOLESKY for Periodic BC;");

        for (test = 1; test <= 2; test++)
        {
            double h = test switch
            {
                1 => 1.0,
                _ => 1.0 / (n + 1)
            };

            Console.WriteLine("");
            Console.WriteLine("  Using spacing H = " + h + "");

            double[] a = L1DD.l1dd(n, h);
            double[] c = L1DD.l1dd_cholesky(n, h);
            typeMethods.r8mat_print(n, n, c, "  L1DD Cholesky factor:");
            double err = Cholesky.cholesky_upper_error(n, a, c);
            Console.WriteLine("");
            Console.WriteLine("  L1DD Cholesky error = " + err + "");

            a = L1DN.l1dn(n, h);
            c = L1DN.l1dn_cholesky(n, h);
            typeMethods.r8mat_print(n, n, c, "  L1DN Cholesky factor:");
            err = Cholesky.cholesky_upper_error(n, a, c);
            Console.WriteLine("");
            Console.WriteLine("  L1DN Cholesky error = " + err + "");

            a = L1ND.l1nd(n, h);
            c = L1ND.l1nd_cholesky(n, h);
            typeMethods.r8mat_print(n, n, c, "  L1ND Cholesky factor:");
            err = Cholesky.cholesky_upper_error(n, a, c);
            Console.WriteLine("");
            Console.WriteLine("  L1ND Cholesky error = " + err + "");

            a = L1NN.l1nn(n, h);
            c = L1NN.l1nn_cholesky(n, h);
            typeMethods.r8mat_print(n, n, c, "  L1NN Cholesky factor:");
            err = Cholesky.cholesky_upper_error(n, a, c);
            Console.WriteLine("");
            Console.WriteLine("  L1NN Cholesky error = " + err + "");

            a = L1PP.l1pp(n, h);
            c = L1PP.l1pp_cholesky(n, h);
            typeMethods.r8mat_print(n, n, c, "  L1PP Cholesky factor:");
            err = Cholesky.cholesky_upper_error(n, a, c);
            Console.WriteLine("");
            Console.WriteLine("  L1PP Cholesky error = " + err + "");
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests L1DD_LU and similar routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Compute LU factors for the Laplacian:");
        Console.WriteLine("  L1DD_LU for Dirichlet/Dirichlet BC;");
        Console.WriteLine("  L1DN_LU for Dirichlet/Neumann BC;");
        Console.WriteLine("  L1ND_LU for Neumann/Dirichlet BC;");
        Console.WriteLine("  L1NN_LU for Neumann/Neumann BC;");
        Console.WriteLine("  L1PP_LU for Periodic BC;");

        double[] l = new double[n * n];
        double[] u = new double[n * n];

        for (test = 1; test <= 2; test++)
        {
            double h = test switch
            {
                1 => 1.0,
                _ => 1.0 / (n + 1)
            };

            Console.WriteLine("");
            Console.WriteLine("  Using spacing H = " + h + "");

            double[] a = L1DD.l1dd(n, h);
            L1DD.l1dd_lu(n, h, ref l, ref u);
            typeMethods.r8mat_print(n, n, l, "  L1DD L factor:");
            typeMethods.r8mat_print(n, n, u, "  L1DD U factor:");
            double err = LU.lu_error(n, a, l, u);
            Console.WriteLine("");
            Console.WriteLine("  L1DD LU error = " + err + "");

            a = L1DN.l1dn(n, h);
            L1DN.l1dn_lu(n, h, ref l, ref u);
            typeMethods.r8mat_print(n, n, l, "  L1DN L factor:");
            typeMethods.r8mat_print(n, n, u, "  L1DN U factor:");
            err = LU.lu_error(n, a, l, u);
            Console.WriteLine("");
            Console.WriteLine("  L1DN LU error = " + err + "");

            a = L1ND.l1nd(n, h);
            L1ND.l1nd_lu(n, h, ref l, ref u);
            typeMethods.r8mat_print(n, n, l, "  L1ND L factor:");
            typeMethods.r8mat_print(n, n, u, "  L1ND U factor:");
            err = LU.lu_error(n, a, l, u);
            Console.WriteLine("");
            Console.WriteLine("  L1ND LU error = " + err + "");

            a = L1NN.l1nn(n, h);
            L1NN.l1nn_lu(n, h, ref l, ref u);
            typeMethods.r8mat_print(n, n, l, "  L1NN L factor:");
            typeMethods.r8mat_print(n, n, u, "  L1NN U factor:");
            err = LU.lu_error(n, a, l, u);
            Console.WriteLine("");
            Console.WriteLine("  L1NN LU error = " + err + "");

            a = L1PP.l1pp(n, h);
            L1PP.l1pp_lu(n, h, ref l, ref u);
            typeMethods.r8mat_print(n, n, l, "  L1PP L factor:");
            typeMethods.r8mat_print(n, n, u, "  L1PP U factor:");
            err = LU.lu_error(n, a, l, u);
            Console.WriteLine("");
            Console.WriteLine("  L1PP LU error = " + err + "");
        }
    }
}