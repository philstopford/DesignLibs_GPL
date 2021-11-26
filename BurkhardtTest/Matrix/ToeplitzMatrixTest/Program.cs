﻿using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace ToeplitzMatrixTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN tests TOEPLITZ_CHOLESKY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOEPLITZ_CHOLESKY_TEST:");
        Console.WriteLine("  Test the TOEPLITZ_CHOLESKY library.");

        t_cholesky_lower_test();
        toep_cholesky_lower_test();
        toeplitz_cholesky_lower_test();
        t_cholesky_upper_test();
        toep_cholesky_upper_test();
        toeplitz_cholesky_upper_test();
        /*
        Terminate.
        */
        Console.WriteLine("");
        Console.WriteLine("TOEPLITZ_CHOLESKY_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void t_cholesky_lower_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_CHOLESKY_LOWER_TEST tests T_CHOLESKY_LOWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 3;
        double[] t =
        {
            1.0, 0.5, -0.375
        };

        Console.WriteLine("");
        Console.WriteLine("T_CHOLESKY_LOWER_TEST");
        Console.WriteLine("  T_CHOLESKY_LOWER computes the lower Cholesky factor L");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined by the first row.");

        typeMethods.r8vec_print(n, t, "  First row of Toeplitz matrix T:");

        double[] l = ToeplitzMatrix.t_cholesky_lower(n, t);
        typeMethods.r8mat_print(n, n, l, "  Computed lower Cholesky factor L:");

        double[] b = typeMethods.r8mat_mmt_new(n, n, n, l, l);
        typeMethods.r8mat_print(n, n, b, "  Product LL':");

    }

    private static void toep_cholesky_lower_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEP_CHOLESKY_LOWER_TEST tests TOEP_CHOLESKY_LOWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] g =
        {
            1.0, 0.0,
            0.5, 0.5,
            -0.375, -0.375
        };
        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("TOEP_CHOLESKY_LOWER_TEST");
        Console.WriteLine("  TOEP_CHOLESKY_LOWER computes the lower Cholesky factor L");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined by a (2,N) array.");

        typeMethods.r8mat_print(2, n, g, "  Compressed Toeplitz matrix G:");

        double[] l = ToeplitzMatrix.toep_cholesky_lower(n, g);
        typeMethods.r8mat_print(n, n, l, "  Computed lower Cholesky factor L:");

        double[] b = typeMethods.r8mat_mmt_new(n, n, n, l, l);
        typeMethods.r8mat_print(n, n, b, "  Product LL':");
    }

    private static void toeplitz_cholesky_lower_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEPLITZ_CHOLESKY_LOWER_TEST tests TOEPLITZ_CHOLESKY_LOWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            1.0, 0.5, -0.375,
            0.5, 1.0, 0.5,
            -0.375, 0.5, 1.0
        };
        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("TOEPLITZ_CHOLESKY_LOWER_TEST");
        Console.WriteLine("  TOEPLITZ_CHOLESKY_LOWER computes the lower Cholesky factor L");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined as an NxN array.");

        typeMethods.r8mat_print(n, n, a, "  Toeplitz matrix A:");

        double[] l = ToeplitzMatrix.toeplitz_cholesky_lower(n, a);
        typeMethods.r8mat_print(n, n, l, "  Computed lower Cholesky factor L:");

        double[] b = typeMethods.r8mat_mmt_new(n, n, n, l, l);
        typeMethods.r8mat_print(n, n, b, "  Product LL':");

    }

    private static void t_cholesky_upper_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_CHOLESKY_UPPER_TEST tests T_CHOLESKY_UPPER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 3;
        double[] t =
        {
            1.0, 0.5, -0.375
        };

        Console.WriteLine("");
        Console.WriteLine("T_CHOLESKY_UPPER_TEST");
        Console.WriteLine("  T_CHOLESKY_UPPER computes the upper Cholesky factor R");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined by the first row.");

        typeMethods.r8vec_print(n, t, "  First row of Toeplitz matrix T:");

        double[] r = ToeplitzMatrix.t_cholesky_upper(n, t);
        typeMethods.r8mat_print(n, n, r, "  Computed upper Cholesky factor R:");

        double[] b = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, b, "  Product R'R:");
    }

    private static void toep_cholesky_upper_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEP_CHOLESKY_UPPER_TEST tests TOEP_CHOLESKY_UPPER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] g =
        {
            1.0, 0.0,
            0.5, 0.5,
            -0.375, -0.375
        };
        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("TOEP_CHOLESKY_UPPER_TEST");
        Console.WriteLine("  TOEP_CHOLESKY_UPPER computes the upper Cholesky factor L");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined by a (2,N) array.");

        typeMethods.r8mat_print(2, n, g, "  Compressed Toeplitz matrix G:");

        double[] r = ToeplitzMatrix.toep_cholesky_upper(n, g);
        typeMethods.r8mat_print(n, n, r, "  Computed upper Cholesky factor R:");

        double[] b = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, b, "  Product R'R:");
    }

    private static void toeplitz_cholesky_upper_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEPLITZ_CHOLESKY_UPPER_TEST tests TOEPLITZ_CHOLESKY_UPPER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            1.0, 0.5, -0.375,
            0.5, 1.0, 0.5,
            -0.375, 0.5, 1.0
        };
        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("TOEPLITZ_CHOLESKY_UPPER_TEST");
        Console.WriteLine("  TOEPLITZ_CHOLESKY_UPPER computes the upper Cholesky factor L");
        Console.WriteLine("  of a positive definites symmetric Toeplitz matrix");
        Console.WriteLine("  defined as an NxN array.");

        typeMethods.r8mat_print(n, n, a, "  Toeplitz matrix A:");

        double[] r = ToeplitzMatrix.toeplitz_cholesky_upper(n, a);
        typeMethods.r8mat_print(n, n, r, "  Computed upper Cholesky factor R:");

        double[] b = typeMethods.r8mat_mtm_new(n, n, n, r, r);
        typeMethods.r8mat_print(n, n, b, "  Product R'R:");
    }
}