using System;
using System.Numerics;
using Burkardt.Sequence;

namespace PolPakTest;

public static class tribonacciTest
{
    public static void tribonacci_direct_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_direct_test() tests tribonacci_direct().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        int t;

        Console.WriteLine("");
        Console.WriteLine("tribonacci_direct_test():");
        Console.WriteLine("  tribonacci_direct() computes the Tribonacci sequence.");
        Console.WriteLine("");

        n = 20;
        for (i = 1; i <= n; i++)
        {
            t = Tribonacci.tribonacci_direct(i);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    public static void tribonacci_recursive_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_recursive_test() tests tribonacci_recursive().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] f;
        int i;
        int n;

        Console.WriteLine("");
        Console.WriteLine("tribonacci_recursive_test():");
        Console.WriteLine("  tribonacci_recursive() computes the Tribonacci sequence.");
        Console.WriteLine("");

        n = 22;

        f = Tribonacci.tribonacci_recursive(n);

        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void tribonacci_roots_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_roots_test() tests tribonacci_roots().
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha = 0;
        Complex beta = new();
        Complex gamma = new();
        double p = 0;
        Complex pc = new();

        Console.WriteLine("");
        Console.WriteLine("tribonacci_roots_test():");
        Console.WriteLine("  tribonacci_roots() computes the Tribonacci roots.");
        Console.WriteLine("");

        Tribonacci.tribonacci_roots(ref alpha, ref beta, ref gamma);

        p = Math.Pow(alpha, 3) - Math.Pow(alpha, 2) - alpha - 1.0;
        Console.WriteLine("  alpha = " + alpha
                                       + ", p(alpha) = " + p + "");

        pc = Complex.Pow(beta, 3) - Complex.Pow(beta, 2) - beta - 1.0;
        Console.WriteLine("  beta = " + beta
                                      + ", p(beta) = " + pc + "");

        pc = Complex.Pow(gamma, 3) - Complex.Pow(gamma, 2) - gamma - 1.0;
        Console.WriteLine("  gamma = " + gamma
                                       + ", p(gamma) = " + pc + "");

    }

}