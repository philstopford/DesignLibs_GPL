using System;
using Burkardt.CholeskyNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HankelSPDCholeskyTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hankel_spd_test tests hankel_spd.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("hankel_spd_test");
            
        Console.WriteLine("  Test the hankel_spd library.");

        hankel_spd_cholesky_lower_test01();
        hankel_spd_cholesky_lower_test02();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("hankel_spd_test");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void hankel_spd_cholesky_lower_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hankel_spd_cholesky_lower_test01 tests hankel_spd_cholesky_lower.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("hankel_spd_cholesky_lower_test01");
        Console.WriteLine("  hankel_spd_cholesky_lower computes a lower Cholesky");
        Console.WriteLine("  matrix L such that the matrix H = L * L' is a");
        Console.WriteLine("  symmetric positive definite (SPD) Hankel matrix.");

        int n = 5;
        //
        //  Example 1:
        //
        double[] lii = new double[n];
        for (i = 0; i < n; i++)
        {
            lii[i] = 1.0;
        }

        double[] liim1 = new double[n - 1];
        for (i = 0; i < n - 1; i++)
        {
            liim1[i] = 1.0;
        }

        double[] l = HankelSPDCholesky.hankel_spd_cholesky_lower(n, lii, liim1);

        typeMethods.r8mat_print(n, n, l, "  The Cholesky factor L:");

        double[] h = typeMethods.r8mat_mmt_new(n, n, n, l, l);

        typeMethods.r8mat_print(n, n, h, "  The Hankel matrix H = L * L':");

        //
        //  Example 2:
        //
        lii = new double[n];
        for (i = 0; i < n; i++)
        {
            lii[i] = i + 1;
        }

        liim1 = new double[n - 1];
        for (i = 0; i < n - 1; i++)
        {
            liim1[i] = n - 1 - i;
        }

        l = HankelSPDCholesky.hankel_spd_cholesky_lower(n, lii, liim1);

        typeMethods.r8mat_print(n, n, l, "  The Cholesky factor L:");

        h = typeMethods.r8mat_mmt_new(n, n, n, l, l);

        typeMethods.r8mat_print(n, n, h, "  The Hankel matrix H = L * L':");

        //
        //  Example 3:
        //
        int seed = 123456789;
        lii = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        liim1 = UniformRNG.r8vec_uniform_01_new(n - 1, ref seed);

        l = HankelSPDCholesky.hankel_spd_cholesky_lower(n, lii, liim1);

        typeMethods.r8mat_print(n, n, l, "  The Cholesky factor L:");

        h = typeMethods.r8mat_mmt_new(n, n, n, l, l);

        typeMethods.r8mat_print(n, n, h, "  The Hankel matrix H = L * L':");

    }

    private static void hankel_spd_cholesky_lower_test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    hankel_spd_cholesky_lower_test02 tests hankel_spd_cholesky_lower.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int flag = 0;
        int i;

        Console.WriteLine("");
        Console.WriteLine("hankel_spd_cholesky_lower_test02");
        Console.WriteLine("  hankel_spd_cholesky_lower computes a lower Cholesky");
        Console.WriteLine("  matrix L such that the matrix H = L * L' is a");
        Console.WriteLine("  symmetric positive definite (SPD) Hankel matrix.");

        int n = 5;

        double[] lii = new double[n];
        for (i = 0; i < n; i++)
        {
            lii[i] = 1.0;
        }

        double[] liim1 = new double[n - 1];
        for (i = 0; i < n - 1; i++)
        {
            liim1[i] = 1.0;
        }

        double[] l = HankelSPDCholesky.hankel_spd_cholesky_lower(n, lii, liim1);

        typeMethods.r8mat_print(n, n, l, "  The Cholesky factor L:");

        double[] h = typeMethods.r8mat_mmt_new(n, n, n, l, l);

        typeMethods.r8mat_print(n, n, h, "  The Hankel matrix H = L * L':");

        double[] l2 = typeMethods.r8mat_cholesky_factor(n, h, ref flag);

        typeMethods.r8mat_print(n, n, l2, "  The Cholesky factor L2 of H:");

        double[] h2 = typeMethods.r8mat_mmt_new(n, n, n, l2, l2);

        typeMethods.r8mat_print(n, n, h2, "  The Hankel matrix H2 = L2 * L2':");
    }
}