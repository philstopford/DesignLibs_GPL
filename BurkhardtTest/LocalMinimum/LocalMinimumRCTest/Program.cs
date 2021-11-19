using System;
using Burkardt.SolveNS;

namespace LocalMinimumRCTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    local_min_rc_test() tests local_min_rc().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;

        Console.WriteLine("");
        Console.WriteLine("local_min_rc_test():");
        Console.WriteLine("  local_min_rc() seeks a local minimizer of a function F(X)");
        Console.WriteLine("  in an interval [A,B], using reverse communication.");


        a = 0.0;
        b = 3.141592653589793;
        example_test(a, b, g_01, "g_01(x) = ( x - 2 ) * ( x - 2 ) + 1");

        a = 0.0;
        b = 1.0;
        example_test(a, b, g_02, "g_02(x) = x * x + exp ( - x )");

        a = -2.0;
        b = 2.0;
        example_test(a, b, g_03, "g_03(x) = x^4 + 2x^2 + x + 3");

        a = 0.0001;
        b = 1.0;
        example_test(a, b, g_04, "g_04(x) = exp ( x ) + 1 / ( 100 x )");

        a = 0.0002;
        b = 2.0;
        example_test(a, b, g_05, "g_05(x) = exp ( x ) - 2x + 1/(100x) - 1/(1000000x^2)");

        a = 1.8;
        b = 1.9;
        example_test(a, b, g_06, "g_06(x) = - x sin ( 10 pi x ) - 1");

        a = 0.0;
        b = 2.0;
        example_test(a, b, g_07, "g_07(x) = 2x^4 - 4x^2 + x + 20");

        a = -2.0;
        b = 0.0;
        example_test(a, b, g_07, "g_07(x) = 2x^4 - 4x^2 + x + 20");

        Console.WriteLine("");
        Console.WriteLine("local_min_rc_test():");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void example_test(double a, double b, Func<double, double> f, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    example_test() tests local_min_rc() on one test function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2008
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
        double a2;
        double arg;
        double b2;
        int status;
        int step;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("  " + title + "");
        Console.WriteLine("");
        Console.WriteLine("  Step      X                          F(X)");
        Console.WriteLine("");
        step = 0;

        arg = a;
        value = f(arg);
        Console.WriteLine("  " + step.ToString().PadLeft(4)
                               + "  " + arg.ToString("0.################").PadLeft(24)
                               + "  " + value.ToString("0.################").PadLeft(24) + "");

        arg = b;
        value = f(arg);
        Console.WriteLine("  " + step.ToString().PadLeft(4)
                               + "  " + arg.ToString("0.################").PadLeft(24)
                               + "  " + value.ToString("0.################").PadLeft(24) + "");

        a2 = a;
        b2 = b;
        status = 0;

        LocalMinimum.LocalMinimumData data = new();

        for (;;)
        {
            arg = LocalMinimum.local_min_rc(ref data, ref a2, ref b2, ref status, value);

            if (status < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("example_test(): Fatal error!");
                Console.WriteLine("  LOCAL_MIN_RC returned negative status.");
                break;
            }

            value = f(arg);

            step += 1;
            Console.WriteLine("  " + step.ToString().PadLeft(4)
                                   + "  " + arg.ToString("0.################").PadLeft(24)
                                   + "  " + value.ToString("0.################").PadLeft(24) + "");

            if (50 < step)
            {
                Console.WriteLine("");
                Console.WriteLine("example_test() - Fatal error!");
                Console.WriteLine("  Too many steps!");
                break;
            }

            if (status == 0)
            {
                break;
            }
        }
    }

    private static double g_01(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    G_01 evaluates (x-2)^2 + 1.
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_01, the value of the function at X.
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
        //    G_02 evaluates x^2 + exp ( - x ).
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_02, the value of the function at X.
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
        //    G_03 evaluates x^4+2x^2+x+3.
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_03, the value of the function at X.
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
        //    G_04 evaluates exp(x)+1/(100X)
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_04, the value of the function at X.
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
        //    G_05 evaluates exp(x) - 2x + 1/(100x) - 1/(1000000x^2)
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_05, the value of the function at X.
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
        //    G_06 evaluates - x sin ( 10 pi x ) - 1
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
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_06, the value of the function at X.
        //
    {
        const double r8_pi = 3.141592653589793;
        double value = 0;

        value = -x * Math.Sin(10.0 * r8_pi * x) - 1.0;

        return value;
    }

    private static double g_07(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    G_07 evaluates 2x^4 - 4x^2 + x + 20
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the point at which F is to be evaluated.
        //
        //    Output, double G_07, the value of the function at X.
        //
    {
        double value = 0;

        value = 2.0 * x * x * x * x - 4.0 * x * x + x + 20.0;

        return value;
    }
}