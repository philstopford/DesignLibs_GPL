using System;
using Burkardt.Bisection;

namespace BisectionIntegerTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BISECTION_INTEGER_TEST.
        //
        //  Discussion:
        //
        //    BISECTION_INTEGER_TEST tests the BISECTION_INTEGER library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BISECTION_INTEGER_TEST");
        Console.WriteLine("  Test the BISECTION_INTEGER library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("BISECTION_INTEGER_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BISECTION_INTEGER;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int b = 0;
        int c = 0;
        int fc = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  BISECTION_INTEGER attempts to locate an integer root C");
        Console.WriteLine("  of an equation F(C) = 0.");
        Console.WriteLine("  The user supplies a change of sign interval [A,B].");
        Console.WriteLine("  The function considered here has two real roots");
        Console.WriteLine("  as well as an integer root, so the algorithm can");
        Console.WriteLine("  fail depending on how the change of sign interval is chosen.");

        a = 4;
        b = 100;

        Console.WriteLine("");
        Console.WriteLine("  The initial change of sign interval is:");
        Console.WriteLine("  F(" + a + ") = " + f01(a) + "");
        Console.WriteLine("  F(" + b + ") = " + f01(b) + "");

        Integer.bisection_integer(f01, ref a, ref b, ref c, ref fc);

        switch (fc)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("  An exact root was found at C = " + c + "");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  An exact root was NOT found.");
                Console.WriteLine("  The change of sign interval is now:");
                Console.WriteLine("  F(" + a + ") = " + f01(a) + "");
                Console.WriteLine("  F(" + b + ") = " + f01(b) + "");
                break;
        }

        a = -10;
        b = 15;

        Console.WriteLine("");
        Console.WriteLine("  The initial change of sign interval is:");
        Console.WriteLine("  F(" + a + ") = " + f01(a) + "");
        Console.WriteLine("  F(" + b + ") = " + f01(b) + "");

        Integer.bisection_integer(f01, ref a, ref b, ref c, ref fc);

        switch (fc)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("  An exact root was found at C = " + c + "");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  An exact root was NOT found.");
                Console.WriteLine("  The change of sign interval is now:");
                Console.WriteLine("  F(" + a + ") = " + f01(a) + "");
                Console.WriteLine("  F(" + b + ") = " + f01(b) + "");
                break;
        }
    }

    private static int f01(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F01 is a test function.
        //
        //  Discussion:
        //
        //    The polynomial has roots 1/2, 7/2, and 10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the argument.
        //
        //    Output, int F01, the function value.
        //
    {
        int value = (2 * n - 7) * (2 * n - 1) * (n - 10);

        return value;
    }
}