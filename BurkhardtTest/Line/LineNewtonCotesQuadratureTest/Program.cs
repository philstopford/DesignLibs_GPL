using System;

namespace LineNewtonCotesQuadratureTest;

using NewtonCotesQuadrature = Burkardt.LineNS.NewtonCotesQuadrature;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LINE_NCC_RULE_TEST.
        //
        //  Discussion:
        //
        //    LINE_NCC_RULE_TEST tests the LINE_NCC_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LINE_NCC_RULE_TEST");
        Console.WriteLine("  Test the LINE_NCC_RULE library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("LINE_NCC_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 computes and prints NCC rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int i;
        int n;
        double[] w;
        double w_sum;
        double[] x;

        a = -1.0;
        b = +1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  LINE_NCC_RULE computes the Newton-Cotes (closed) rule");
        Console.WriteLine("  using N equally spaced points for an interval [A,B].");

        for (n = 1; n <= 12; n++)
        {
            x = new double[n];
            w = new double[n];

            NewtonCotesQuadrature.line_ncc_rule(n, a, b, x, ref w);
            Console.WriteLine("");
            Console.WriteLine("  Newton-Cotes (Closed) Rule #" + n + "");
            Console.WriteLine("   I       X(I)            W(I)");
            Console.WriteLine("");
            w_sum = 0.0;
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(14)
                                       + "  " + w[i].ToString().PadLeft(14) + "");
                w_sum += Math.Abs(w[i]);
            }

            Console.WriteLine("        Sum(|W)|) =  " + w_sum.ToString().PadLeft(14) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses NCC rules to estimate the integral of exp(x) from 0 to 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double error;
        double exact;
        int i;
        int n;
        double q;
        double[] w;
        double[] x;

        a = 0.0;
        b = +1.0;
        exact = Math.Exp(b) - Math.Exp(a);

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use a sequence of NCC rules to compute an estimate Q");
        Console.WriteLine("  of the integral:");
        Console.WriteLine("    I = integral ( 0 <= x <= 1 ) exp(x) dx.");
        Console.WriteLine("  The exact value is:");
        Console.WriteLine("    I = " + exact + "");

        Console.WriteLine("");
        Console.WriteLine("   N       Q             |Q-I|");
        Console.WriteLine("");

        for (n = 1; n <= 22; n++)
        {
            x = new double[n];
            w = new double[n];

            NewtonCotesQuadrature.line_ncc_rule(n, a, b, x, ref w);

            q = 0.0;
            for (i = 0; i < n; i++)
            {
                q += w[i] * Math.Exp(x[i]);
            }

            error = Math.Abs(exact - q);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + q.ToString().PadLeft(14)
                                   + "  " + error.ToString().PadLeft(14) + "");
        }
    }
}