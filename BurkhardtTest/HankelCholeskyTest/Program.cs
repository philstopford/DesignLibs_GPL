using System;
using Burkardt.CholeskyNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HankelCholeskyTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HANKEL_CHOLESKY_TEST tests HANKEL_CHOLESKY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HANKEL_CHOLESKY_TEST");
        Console.WriteLine("  Test the HANKEL_CHOLESKY library.");

        hankel_cholesky_upper_test();

        Console.WriteLine("");
        Console.WriteLine("HANKEL_CHOLESKY_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void hankel_cholesky_upper_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HANKEL_CHOLESKY_UPPER_TEST tests HANKEL_CHOLESKY_UPPER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int flag = 0;
        int i;
        int j;

        int n = 5;

        Console.WriteLine("");
        Console.WriteLine("HANKEL_CHOLESKY_UPPER_TEST");
        Console.WriteLine("  HANKEL_CHOLESKY_UPPER is given a Hankel matrix H and");
        Console.WriteLine("  computes an upper triangular matrix R such that");
        Console.WriteLine("  H = R' * R");
        //
        //  Get a Hankel matrix that is symmetric positive definite.
        //
        int seed = 123456789;
        double[] lii = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        
        double[] liim1 = UniformRNG.r8vec_uniform_01_new(n - 1, ref seed);
        double[] l = HankelCholesky.hankel_spd_cholesky_lower(n, lii, liim1);
        double[] h = typeMethods.r8mat_mmt_new(n, n, n, l, l);
        typeMethods.r8mat_print(n, n, h, "  The Hankel matrix H:");
        //
        //  Compute R using R8MAT_CHOLESKY_FACTOR_UPPER.
        //
        double[] r1 = typeMethods.r8mat_cholesky_factor_upper(n, h, ref flag);
        if (flag != 0)
        {
            Console.WriteLine("");
            Console.WriteLine(" R8MAT_CHOLESKY_FACTOR_UPPER says H is not positive definite.");
        }
        else
        {
            typeMethods.r8mat_print(n, n, r1, "  R computed by R8MAT_CHOLESKY_FACTOR_UPPER:");
        }

        //
        //  Compute R using HANKEL_CHOLESKY.
        //
        double[] hanti = new double[2 * n - 1];
        for (i = 0; i < n; i++)
        {
            hanti[i] = h[i + 0 * n];
        }

        for (j = 1; j < n; j++)
        {
            hanti[n - 1 + j] = h[n - 1 + j * n];
        }

        double[] r2 = HankelCholesky.hankel_cholesky_upper(n, hanti);
        typeMethods.r8mat_print(n, n, r2, "  R computed by HANKEL_CHOLESKY:");
    }
}