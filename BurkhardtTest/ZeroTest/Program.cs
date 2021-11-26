using System;
using System.Globalization;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace ZeroTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    zero_test() tests zero().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("zero_test():");
        Console.WriteLine("  zero() seeks a root of a function F(X)");
        Console.WriteLine("  in an interval [A,B].");

        double a = 1.0;
        double b = 2.0;
        zero_example(a, b, f_01, "f_01(x) = sin ( x ) - x / 2");

        a = 0.0;
        b = 1.0;
        zero_example(a, b, f_02, "f_02(x) = 2 * x - exp ( - x )");

        a = -1.0;
        b = 0.5;
        zero_example(a, b, f_03, "f_03(x) = x * exp ( - x )");

        a = 0.0001;
        b = 20.0;
        zero_example(a, b, f_04, "f_04(x) = exp ( x ) - 1 / ( 100 * x * x )");

        a = -5.0;
        b = 2.0;
        zero_example(a, b, f_05, "f_05(x) = (x+3) * (x-1) * (x-1)");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("zero_test():");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void zero_example(double a, double b, Func<double, double> f, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    zero_example() tests zero() on one test function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double A, B, the endpoints of the change of sign interval.
        //
        //    double F ( double x ), the name of a user-supplied
        //    function which evaluates the function whose zero is being sought.
        //
        //    string TITLE, a title for the problem.
        //
    {
        int calls = 0;

        double t = typeMethods.r8_epsilon();

        double z = Zero.zero(a, b, t, f, ref calls);
        double fz = f(z);
        double fa = f(a);
        double fb = f(b);

        Console.WriteLine("");
        Console.WriteLine("  " + title + "");
        Console.WriteLine("");
        Console.WriteLine("           A                 Z             B");
        Console.WriteLine("         F(A)              F(Z)          F(B)");
        Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + z.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  " + fa.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + fz.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + fb.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("  Number of calls to F = " + calls + "");
    }

    private static double f_01(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_01 evaluates sin ( x ) - x / 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the point at which F is to be evaluated.
        //
        //  Output:
        //
        //    double F_01, the value of the function at X.
        //
    {
        double value = 0;

        value = Math.Sin(x) - 0.5 * x;

        return value;
    }

    private static double f_02(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_02 evaluates 2*x-exp(-x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the point at which F is to be evaluated.
        //
        //  Output:
        //
        //    double F_02, the value of the function at X.
        //
    {
        double value = 0;

        value = 2.0 * x - Math.Exp(-x);

        return value;
    }

    private static double f_03(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_03 evaluates x*exp(-x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the point at which F is to be evaluated.
        //
        //  Output:
        //
        //    double F_03, the value of the function at X.
        //
    {
        double value = 0;

        value = x * Math.Exp(-x);

        return value;
    }

    private static double f_04(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_04 evaluates exp(x) - 1 / (100*x*x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the point at which F is to be evaluated.
        //
        //  Output:
        //
        //    double F_04, the value of the function at X.
        //
    {
        double value = 0;

        value = Math.Exp(x) - 1.0 / 100.0 / x / x;

        return value;
    }

    private static double f_05(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_05 evaluates (x+3)*(x-1)*(x-1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    double X, the point at which F is to be evaluated.
        //
        //  Output:
        //
        //    double F_05, the value of the function at X.
        //
    {
        double value = 0;

        value = (x + 3.0) * (x - 1.0) * (x - 1.0);

        return value;
    }
}