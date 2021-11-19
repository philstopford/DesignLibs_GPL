using System;
using System.Globalization;
using Burkardt.CircleNS;

namespace CircleRulesTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CIRCLE_RULE_TEST.
        //
        //  Discussion:
        //
        //    CIRCLE_RULE_TEST tests the CIRCLE_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CIRCLE_RULE:");

        Console.WriteLine("  Test the CIRCLE_RULE library.");

        int nt = 8;
        test01(nt);

        nt = 32;
        test01(nt);

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int nt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CIRCLE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[2];
        int e1;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  CIRCLE_RULE can compute a rule Q(f) for the unit circle");
        Console.WriteLine("  using NT equally spaced angles.");
        Console.WriteLine("  Estimate integrals I(f) where f = x^e(1) * y^e(2)");
        Console.WriteLine("  using " + nt + " points.");
        //
        //  Compute the quadrature rule.
        //
        double[] w = new double[nt];
        double[] t = new double[nt];

        QuadratureRule.circle_rule(nt, ref w, ref t);
        //
        //  Apply it to integrands.
        //
        Console.WriteLine("");
        Console.WriteLine("  E(1)  E(2)    I(f)            Q(f)");
        Console.WriteLine("");
        //
        //  Specify a monomial.
        //
        for (e1 = 0; e1 <= 6; e1 += 2)
        {
            e[0] = e1;

            int e2;
            for (e2 = e1; e2 <= 6; e2 += 2)
            {
                e[1] = e2;

                double q = 0.0;
                int i;
                for (i = 0; i < nt; i++)
                {
                    double x = Math.Cos(t[i]);
                    double y = Math.Sin(t[i]);
                    q += w[i] * Math.Pow(x, e[0]) * Math.Pow(y, e[1]);
                }

                q = 2.0 * Math.PI * q;

                double exact = Integrals.circle01_monomial_integral(e);

                Console.WriteLine("  " + e[0].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + e[1].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }
}