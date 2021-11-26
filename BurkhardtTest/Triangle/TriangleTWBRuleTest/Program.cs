using System;
using System.Globalization;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace TriangleTWBRuleTest;

using TriMonomial = Monomial;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    triangle_twb_rule_test tests triangle_twb_rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("triangle_twb_rule_test");
        Console.WriteLine("  Test triangle_twb_rule.");

        const int degree_max = 5;
        triangle_unit_quad_test(degree_max);

        Console.WriteLine("");
        Console.WriteLine("triangle_twb_rule_test");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void triangle_unit_quad_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    triangle_unit_quad_test tests rules for the unit triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the monomials to check.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("triangle_unit_quad_test");
        Console.WriteLine("  Approximate monomial integrals in triangle with TWB rules.");

        int degree = 0;
        int ex = 0;
        int ey = degree;

        while (true)
        {
            Console.WriteLine("");
            Console.WriteLine("  Monomial: x^" + ex + " y^" + ey + "");

            int strength;
            double q;
            for (strength = 1; strength <= 25; strength++)
            {
                int n = TWBRule.twb_rule_n(strength);
                switch (n)
                {
                    case < 1:
                        continue;
                }

                double[] w = TWBRule.twb_rule_w(strength);
                double[] x = TWBRule.twb_rule_x(strength);
                double[] y = TWBRule.twb_rule_y(strength);
                double[] v = Burkardt.MonomialNS.Monomial.monomial_value_2d(n, ex, ey, x, y);
                q = typeMethods.r8vec_dot_product(n, w, v);
                Console.WriteLine("  " + strength.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            }

            q = TriMonomial.triangle_unit_monomial(ex, ey);
            Console.WriteLine("   Exact          " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (ex < degree)
            {
                ex += 1;
                ey -= 1;
            }
            else if (degree < degree_max)
            {
                degree += 1;
                ex = 0;
                ey = degree;
            }
            else
            {
                break;
            }
        }

    }
}