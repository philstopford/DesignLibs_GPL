using System;
using System.Globalization;
using Burkardt.ODENS;

namespace StochasticRungeKuttaTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STOCHASTIC_RK_TEST.
        //
        //  Discussion:
        //
        //    STOCHASTIC_RK_TEST tests the STOCHASTIC_RK library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_RK_TEST");
        Console.WriteLine("  Test the STOCHASTIC_RK library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_RK_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests RK1_TI_STEP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double t0 = 0.0;
        const double tn = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  RK1_TI_STEP uses a first order RK method");
        Console.WriteLine("  for a problem whose right hand side does not");
        Console.WriteLine("  depend explicitly on time.");

        const int n = 10;
        double[] x = new double[n + 1];
        const double h = (tn - t0) / n;
        const double q = 1.0;
        int seed = 123456789;

        int i = 0;
        double t = t0;
        x[i] = 0.0;

        Console.WriteLine("");
        Console.WriteLine("         I           T             X");
        Console.WriteLine("");
        Console.WriteLine("  " + i.ToString().PadLeft(8)
                               + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        for (i = 1; i <= n; i++)
        {
            t = ((n - i) * t0
                 + i * tn)
                / n;

            x[i] = RungeKutta.rk1_ti_step(x[i - 1], t, h, q, fi, gi, ref seed);

            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static double fi(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FI is a time invariant deterministic right hand side.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double FI, the value.
        //
    {
        const double value = 1.0;

        return value;
    }

    private static double gi(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GI is a time invariant stochastic right hand side.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double GI, the value.
        //
    {
        const double value = 1.0;

        return value;
    }
}