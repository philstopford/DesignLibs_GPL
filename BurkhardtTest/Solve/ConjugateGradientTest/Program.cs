using System;
using Burkardt.RandomMatrix;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ConjugateGradientTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CG_TEST.
        //
        //  Discussion:
        //
        //    CG_TEST tests the CG library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CG_TEST");

        Console.WriteLine("  Test the CG library.");

        r8ge_cg_test();
        r83_cg_test();
        r83s_cg_test();
        r83t_cg_test();
        r8pbu_cg_test();
        r8sd_cg_test();
        r8sp_cg_test();

        Console.WriteLine(" ");
        Console.WriteLine("CG_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine(" ");
    }

    private static void r8ge_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_CG_TEST tests R8GE_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        typeMethods.r8NormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("R8GE_CG_TEST");
        Console.WriteLine("  R8GE_CG applies CG to a full storage matrix.");

        int n = 10;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        int seed = 123456789;
        double[] a = PositiveDefiniteSymmetric.pds_random(n, ref data, ref seed);
        //
        //  Choose a random solution.
        //
        seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r8ge_mv(n, n, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r8ge_cg(n, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r8ge_res(n, n, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r83_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_CG_TEST tests R83_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R83_CG_TEST");
        Console.WriteLine("  R83_CG applies CG to an R83 matrix.");

        int n = 10;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r83_dif2(n, n);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r83_mv(n, n, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r83_cg(n, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r83_res(n, n, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r83s_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_CG_TEST tests R83S_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R83S_CG_TEST");
        Console.WriteLine("  R83S_CG applies CG to an R83S matrix.");

        const int n = 10;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r83s_dif2(n, n);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r83s_mv(n, n, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r83s_cg(n, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r83s_res(n, n, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r83t_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_CG_TEST tests R83T_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R83T_CG_TEST");
        Console.WriteLine("  R83T_CG applies CG to an R83T matrix.");

        const int n = 10;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r83t_dif2(n, n);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r83t_mv(n, n, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r83t_cg(n, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r83t_res(n, n, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r8pbu_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_CG_TEST tests R8PBU_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R8PBU_CG_TEST");
        Console.WriteLine("  R8PBU_CG applies CG to an R8PBU matrix.");

        const int n = 10;
        const int mu = 1;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r8pbu_dif2(n, n, mu);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r8pbu_mv(n, n, mu, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r8pbu_cg(n, mu, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r8pbu_res(n, n, mu, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r8sd_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_CG_TEST tests R8SD_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R8SD_CG_TEST");
        Console.WriteLine("  R8SD_CG applies CG to an R8SD matrix.");

        const int n = 10;

        const int ndiag = 2;
        int[] offset = new int[ndiag];
        offset[0] = 0;
        offset[1] = 1;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r8sd_dif2(n, n, ndiag, offset);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r8sd_mv(n, n, ndiag, offset, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r8sd_cg(n, ndiag, offset, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r8sd_res(n, n, ndiag, offset, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    private static void r8sp_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_CG_TEST tests R8SP_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R8SP_CG_TEST");
        Console.WriteLine("  R8SP_CG applies CG to an R8SP matrix.");

        const int n = 10;

        const int nz_num = 3 * n;
        int[] row = new int[nz_num];
        int[] col = new int[nz_num];
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a = typeMethods.r8sp_dif2(n, n, nz_num, row, col);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r8sp_mv(n, n, nz_num, row, col, a, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        typeMethods.r8sp_cg(n, nz_num, row, col, a, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = typeMethods.r8sp_res(n, n, nz_num, row, col, a, x2, b);
        double r_norm = typeMethods.r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = typeMethods.r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }
}