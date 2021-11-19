using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace LocalMinimumTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    local_min_test() tests local_min().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 June 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;

        Console.WriteLine("");
        Console.WriteLine("local_min_test():");
        Console.WriteLine("  local_min() seeks a local minimizer of a function F(X)");
        Console.WriteLine("  in an interval [A,B].");

        a = 0.0;
        b = 3.141592653589793;
        local_min_example(a, b, g_01,
            "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1");

        a = 0.0;
        b = 1.0;
        local_min_example(a, b, g_02,
            "g_02(x) = x * x + exp ( - x )");

        a = -2.0;
        b = 2.0;
        local_min_example(a, b, g_03,
            "g_03(x) = x^4 + 2x^2 + x + 3");

        a = 0.0001;
        b = 1.0;
        local_min_example(a, b, g_04,
            "g_04(x) = exp ( x ) + 1 / ( 100 x )");

        a = 0.0002;
        b = 2.0;
        local_min_example(a, b, g_05,
            "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)");

        a = 1.8;
        b = 1.9;
        local_min_example(a, b, g_06,
            "g_06(x) = -x*sin(10*pi*x)-1.0");

        a = -1.2;
        b = 2.7;
        local_min_example(a, b, g_07,
            "g_07(x) = max(-2(x-1),8(x-1)) + 25*(x-1)^2");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("local_min_test():");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void local_min_example(double a, double b, Func<double, double> f,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    local_min_example() tests local_min() on one test function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double A, B, the endpoints of the interval.
        //
        //    double F ( double x ), the name of a user-supplied
        //    function, whose local minimum is being sought.
        //
        //    string TITLE, a title for the problem.
        //
    {
        int calls = 0;
        double fa;
        double fb;
        double fx;
        double t;
        double x = 0;

        t = Math.Sqrt(typeMethods.r8_epsilon());

        fx = LocalMinimum.local_min(a, b, t, f, ref x, ref calls);
        fa = f(a);
        fb = f(b);

        Console.WriteLine("");
        Console.WriteLine("  " + title + "");
        Console.WriteLine("");
        Console.WriteLine("           A                 X             B");
        Console.WriteLine("         F(A)              F(X)          F(B)");
        Console.WriteLine("  " + a.ToString().PadLeft(14)
                               + "  " + x.ToString().PadLeft(14)
                               + "  " + b.ToString().PadLeft(14) + "");
        Console.WriteLine("  " + fa
                               + "  " + fx.ToString().PadLeft(14)
                               + "  " + fb.ToString().PadLeft(14) + "");
        Console.WriteLine("  Number of calls to F = " + calls + "");
    }

    private static double g_01(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_01() evaluates (x-2)^2 + 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_01, the value of the function at X.
        //
    {
        double value = 0;

        value = (x - 2.0) * (x - 2.0) + 1.0;

        return value;
    }

    private static double g_02(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_02() evaluates x^2 + exp ( - x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_02, the value of the function at X.
        //
    {
        double value = 0;

        value = x * x + Math.Exp(-x);

        return value;
    }

    private static double g_03(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_03() evaluates x^4+2x^2+x+3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_03, the value of the function at X.
        //
    {
        double value = 0;

        value = ((x * x + 2.0) * x + 1.0) * x + 3.0;

        return value;
    }

    private static double g_04(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_04() evaluates exp(x)+1/(100X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_04, the value of the function at X.
        //
    {
        double value = 0;

        value = Math.Exp(x) + 0.01 / x;

        return value;
    }

    private static double g_05(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_05() evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_05, the value of the function at X.
        //
    {
        double value = 0;

        value = Math.Exp(x) - 2.0 * x + 0.01 / x - 0.000001 / x / x;

        return value;
    }

    private static double g_06(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_06() evaluates - x * sin(10 pi x ) - 1.0;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_06, the value of the function at X.
        //
    {
        double r8_pi = 3.141592653589793;
        double value = 0;

        value = -x * Math.Sin(10.0 * r8_pi * x) - 1.0;

        return value;
    }

    private static double g_07(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    g_07() evaluates max(-2(x-1), 8(x-1)) + 25 (x-1)^2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the evaluation point.
        //
        //  Output:
        //
        //    double G_07, the value of the function at X.
        //
    {
        double value = 0;

        value = Math.Max(-2.0 * (x - 1), 8.0 * (x - 1))
                + 25.0 * Math.Pow(x - 1.0, 2);

        return value;
    }
}