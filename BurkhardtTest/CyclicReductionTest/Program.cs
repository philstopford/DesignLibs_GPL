using System;
using Burkardt.Types;

namespace CyclicReductionTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CYCLIC_REDUCTION_TEST.
        //
        //  Discussion:
        //
        //    CYCLIC_REDUCTION_TEST tests CYCLIC_REDUCTION.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CYCLIC_REDUCTION_TEST");
            
        Console.WriteLine("  Test the CYCLIC_REDUCTION library.");

        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("CYCLIC_REDUCTION_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests R83_CR_FA, R83_CR_SLS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const bool debug = true;
        int i;
        int j;
        const int n = 5;
        const int nb = 2;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  R83_CR_FA factors a real tridiagonal matrix;");
        Console.WriteLine("  R83_CR_SLS solves 1 or more systems.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  Demonstrate multiple system solution method.");
        //
        //  Set the matrix.
        //
        double[] a = new double[3 * n];

        a[0 + 0 * 3] = 0.0;
        for (j = 1; j < n; j++)
        {
            a[0 + j * 3] = -1.0;
        }

        for (j = 0; j < n; j++)
        {
            a[1 + j * 3] = 2.0;
        }

        for (j = 0; j < n - 1; j++)
        {
            a[2 + j * 3] = -1.0;
        }

        a[2 + (n - 1) * 3] = 0.0;

        switch (debug)
        {
            case true:
                typeMethods.r83_print(n, a, "  Input matrix:");
                break;
        }

        //
        //  Factor the matrix once.
        //
        double[] a_cr = typeMethods.r83_cr_fa(n, a);

        switch (debug)
        {
            case true:
                typeMethods.r83_print(2 * n + 1, a_cr, "  Cyclic reduction factor information:");
                break;
        }

        //
        //  Solve 2 systems simultaneously.
        //
        double[] b = new double[n * nb];

        for (i = 0; i < n - 1; i++)
        {
            b[i + 0 * n] = 0.0;
        }

        b[n - 1 + 0 * n] = n + 1;

        b[0 + 1 * n] = 1.0;
        for (i = 1; i < n - 1; i++)
        {
            b[i + 1 * n] = 0.0;
        }

        b[n - 1 + 1 * n] = 1.0;
        //
        //  Solve the linear systems.
        //
        double[] x = typeMethods.r83_cr_sls(n, a_cr, nb, b);

        typeMethods.r8mat_print_some(n, nb, x, 1, 1, 10, nb, "  Solutions:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests R83_CR_FA, R83_CR_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const bool debug = true;
        int j;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For a real tridiagonal matrix,");
        Console.WriteLine("  using CYCLIC REDUCTION,");
        Console.WriteLine("  R83_CR_FA factors;");
        Console.WriteLine("  R83_CR_SL solves a system.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  The matrix is NOT symmetric.");
        //
        //  Set the matrix values.
        //
        double[] a = new double[3 * n];

        a[0 + 0 * 3] = 0.0;
        for (j = 2; j <= n; j++)
        {
            a[0 + (j - 1) * 3] = j;
        }

        for (j = 1; j <= n; j++)
        {
            a[1 + (j - 1) * 3] = 4.0 * j;
        }

        for (j = 1; j <= n - 1; j++)
        {
            a[2 + (j - 1) * 3] = j;
        }

        a[2 + (n - 1) * 3] = 0.0;

        switch (debug)
        {
            case true:
                typeMethods.r83_print(n, a, "  The matrix:");
                break;
        }

        //
        //  Set the desired solution.
        //
        double[] x = typeMethods.r8vec_indicator_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = typeMethods.r83_mxv_new(n, a, x);

        switch (debug)
        {
            case true:
                typeMethods.r8vec_print(n, b, "  The right hand side:");
                break;
        }

        //
        //  Factor the matrix.
        //
        double[] a_cr = typeMethods.r83_cr_fa(n, a);

        switch (debug)
        {
            case true:
                typeMethods.r83_print(2 * n + 1, a_cr, "  The factor information:");
                break;
        }

        //
        //  Solve the linear system.
        //
        x = typeMethods.r83_cr_sl(n, a_cr, b);

        typeMethods.r8vec_print(n, x, "  The solution:");
    }
}