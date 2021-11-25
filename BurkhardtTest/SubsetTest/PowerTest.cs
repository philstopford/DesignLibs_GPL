using System;
using Burkardt;
using Burkardt.Sequence;
using Burkardt.Types;

namespace SubsetTestNS;

public static class PowerTest
{
    public static void power_mod_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST109 tests POWER_MOD;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POWER_MOD_TEST");
        Console.WriteLine("  POWER_MOD computes the remainder of a power");
        Console.WriteLine("  of an integer modulo another integer.");

        int a = 7;
        int n = 50;
        int m = 11;

        Console.WriteLine("");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  mod ( A^N, M ) = " + Helpers.power_mod(a, n, m) + "");

        a = 3;
        n = 118;
        m = 119;

        Console.WriteLine("");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  mod ( A^N, M ) = " + Helpers.power_mod(a, n, m) + "");
    }

    public static void power_series1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_SERIES1_TEST tests POWER_SERIES1;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;

        double[] a = new double[N];
        double[] b = new double[N];
        int i;

        Console.WriteLine("");
        Console.WriteLine("POWER_SERIES1_TEST");
        Console.WriteLine("  POWER_SERIES1 composes a power series;");

        double alpha = 7.0;

        a[0] = 1.0;
        for (i = 1; i < N; i++)
        {
            a[i] = 0.0;
        }

        for (i = 0; i < N; i++)
        {
            b[i] = 0.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Power series of G(x) = (1+F(x))**alpha");
        Console.WriteLine("");
        Console.WriteLine("  N = " + N + "");
        Console.WriteLine("  ALPHA = " + alpha + "");

        typeMethods.r8vec_print(N, a, "  Series for F(x):");

        PowerSeries.power_series1(N, alpha, a, ref b);

        typeMethods.r8vec_print(N, b, "  Series for G(X):");
    }

    public static void power_series2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_SERIES2_TEST tests POWER_SERIES2;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 4;

        double[] a = new double[N];
        double[] b = new double[N];
        int i;

        Console.WriteLine("");
        Console.WriteLine("POWER_SERIES2_TEST");
        Console.WriteLine("  POWER_SERIES2 composes a power series;");

        a[0] = -4.0;
        for (i = 1; i < N; i++)
        {
            a[i] = 0.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Power series of G(x) = exp(F(x))-1");
        Console.WriteLine("");
        Console.WriteLine("  N = " + N + "");

        typeMethods.r8vec_print(N, a, "  Series for F(X):");

        PowerSeries.power_series2(N, a, ref b);

        typeMethods.r8vec_print(N, b, "  Series for G(X):");
    }

    public static void power_series3_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_SERIES3_TEST tests POWER_SERIES3;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 4;

        double[] a = new double[N];
        double[] b = new double[N];
        double[] c = new double[N];

        Console.WriteLine("");
        Console.WriteLine("POWER_SERIES3_TEST");
        Console.WriteLine("  POWER_SERIES3 composes two power series;");

        a[0] = 1.0;
        a[1] = 1.0;
        a[2] = 0.0;
        a[3] = 0.0;

        typeMethods.r8vec_print(N, a, "  Series for F(X):");

        b[0] = 1.0;
        b[1] = 1.0;
        b[2] = 0.0;
        b[3] = 0.0;

        typeMethods.r8vec_print(N, b, "  Series for G(X):");

        PowerSeries.power_series3(N, a, b, ref c);

        typeMethods.r8vec_print(N, c, "  Series for H(X) = G(F(X)):");
    }

    public static void power_series4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_SERIES4_TEST tests POWER_SERIES4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;
        //
        //  The order of arguments for POWER_SERIES4 is a shame.
        //
        double[] a = new double[N];
        double[] b = new double[N];
        double[] c = new double[N];
        int i;

        Console.WriteLine("");
        Console.WriteLine("POWER_SERIES4_TEST");
        Console.WriteLine("  POWER_SERIES4 composes a power series;");
        Console.WriteLine("  Given power series for F(X) and G(X), we compute");
        Console.WriteLine("  the power series of H(x) = G(1/F(x)).");

        for (i = 0; i < N; i++)
        {
            a[i] = 1.0 / (i + 1);
        }

        typeMethods.r8vec_print(N, a, "  Series for F(x):");

        b[0] = 1.0;
        for (i = 1; i < N; i++)
        {
            b[i] = 0.0;
        }

        typeMethods.r8vec_print(N, b, "  Series for G(x):");

        PowerSeries.power_series4(N, a, b, ref c);

        typeMethods.r8vec_print(N, c, "  Series for H(x):");
    }

}