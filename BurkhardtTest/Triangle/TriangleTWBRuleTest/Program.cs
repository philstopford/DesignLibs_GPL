using System;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace TriangleTWBRuleTest
{
    using TriMonomial = Burkardt.TriangleNS.Monomial;

    class Program
    {
        static void Main(string[] args)
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
            int degree_max;

            Console.WriteLine("");
            Console.WriteLine("triangle_twb_rule_test");
            Console.WriteLine("  Test triangle_twb_rule.");

            degree_max = 5;
            triangle_unit_quad_test(degree_max);

            Console.WriteLine("");
            Console.WriteLine("triangle_twb_rule_test");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void triangle_unit_quad_test(int degree_max)

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
            int degree;
            int ex;
            int ey;
            int n;
            double q;
            int strength;
            double[] v;
            double[] w;
            double[] x;
            double[] y;

            Console.WriteLine("");
            Console.WriteLine("triangle_unit_quad_test");
            Console.WriteLine("  Approximate monomial integrals in triangle with TWB rules.");

            degree = 0;
            ex = 0;
            ey = degree;

            while (true)
            {
                Console.WriteLine("");
                Console.WriteLine("  Monomial: x^" + ex + " y^" + ey + "");

                for (strength = 1; strength <= 25; strength++)
                {
                    n = TWBRule.twb_rule_n(strength);
                    if (n < 1)
                    {
                        continue;
                    }

                    w = TWBRule.twb_rule_w(strength);
                    x = TWBRule.twb_rule_x(strength);
                    y = TWBRule.twb_rule_y(strength);
                    v = Burkardt.MonomialNS.Monomial.monomial_value_2d(n, ex, ey, x, y);
                    q = typeMethods.r8vec_dot_product(n, w, v);
                    Console.WriteLine("  " + strength.ToString().PadLeft(6)
                                           + "  " + n.ToString().PadLeft(6)
                                           + "  " + q.ToString().PadLeft(14) + "");

                }

                q = TriMonomial.triangle_unit_monomial(ex, ey);
                Console.WriteLine("   Exact          " + q.ToString().PadLeft(14) + "");

                if (ex < degree)
                {
                    ex = ex + 1;
                    ey = ey - 1;
                }
                else if (degree < degree_max)
                {
                    degree = degree + 1;
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
}