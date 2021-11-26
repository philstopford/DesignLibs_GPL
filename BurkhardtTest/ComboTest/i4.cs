using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ComboTest;

internal static partial class Program
{
    private static void i4_choose_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHOOSE_TEST tests I4_CHOOSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("I4_CHOOSE_TEST");
        Console.WriteLine("  I4_CHOOSE computes binomial coefficients.");

        for (i = -1; i <= 5; i++)
        {
            int j;
            for (j = -1; j <= 5; j++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + j.ToString().PadLeft(4)
                                       + "  " + typeMethods.i4_choose(i, j).ToString().PadLeft(12) + "");
            }
        }
    }

    private static void i4_factorial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL_TEST tests I4_FACTORIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fx = 0;
        int n = 0;
        int x = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_FACTORIAL_TEST:");
        Console.WriteLine("  I4_FACTORIAL evaluates the factorial function.");
        Console.WriteLine("");
        Console.WriteLine("     X       Exact F       FACTORIAL(X)");
        Console.WriteLine("");
            
        for (;;)
        {
            typeMethods.i4_factorial_values(ref n, ref x, ref fx);

            if (n == 0)
            {
                break;
            }

            if (x <= 0.0)
            {
                continue;
            }

            int fx2 = typeMethods.i4_factorial(x);

            Console.WriteLine("  " + x.ToString().PadLeft(4)
                                   + "  " + fx.ToString().PadLeft(12)
                                   + "  " + fx2.ToString().PadLeft(12) + "");
        }
    }

    private static void i4_fall_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FALL_TEST tests I4_FALL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int f1 = 0;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_FALL_TEST");
        Console.WriteLine("  I4_FALL evaluates the falling factorial function.");
        Console.WriteLine("");
        Console.WriteLine("         M         N     Exact  I4_Fall(M,N)");
        Console.WriteLine("");

        int n_data = 0;

        while (true)
        {
            typeMethods.i4_fall_values(ref n_data, ref m, ref n, ref f1);

            if (n_data == 0)
            {
                break;
            }

            int f2 = typeMethods.i4_fall(m, n);

            Console.WriteLine("  " + m.ToString().PadLeft(8)
                                   + "  " + n.ToString().PadLeft(8)
                                   + "  " + f1.ToString().PadLeft(8)
                                   + "  " + f2.ToString().PadLeft(8) + "");
        }
    }

    private static void i4_uniform_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int a = -100;
        const int b = 200;
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4_UNIFORM_AB_TEST");
        Console.WriteLine("  I4_UNIFORM_AB computes pseudorandom values");
        Console.WriteLine("  in an interval [A,B].");

        Console.WriteLine("");
        Console.WriteLine("  The lower endpoint A = " + a + "");
        Console.WriteLine("  The upper endpoint B = " + b + "");
        Console.WriteLine("  The initial seed is " + seed + "");
        Console.WriteLine("");

        for (i = 1; i <= 20; i++)
        {
            int j = UniformRNG.i4_uniform_ab(a, b, ref seed);

            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + j.ToString().PadLeft(8) + "");
        }
    }
}