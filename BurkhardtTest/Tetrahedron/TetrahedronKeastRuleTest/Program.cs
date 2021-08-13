using System;
using Burkardt;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetrahedronKeastRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TETRAHEDRON_KEAST_RULE_TEST.
            //
            //  Discussion:
            //
            //    TETRAHEDRON_KEAST_RULE_TEST tests the TETRAHEDRON_KEAST_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_KEAST_RULE_TEST:");
            Console.WriteLine("  Test the TETRAHEDRON_KEAST_RULE library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();

            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_KEAST_RULE_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests KEAST_RULE_NUM, KEAST_DEGREE, KEAST_ORDER_NUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree;
            int order_num;
            int rule;
            int rule_num;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  KEAST_RULE_NUM returns the number of rules;");
            Console.WriteLine("  KEAST_DEGREE returns the degree of a rule;");
            Console.WriteLine("  KEAST_ORDER_NUM returns the order of a rule.");

            rule_num = KeastRule.keast_rule_num();

            Console.WriteLine("");
            Console.WriteLine("  Number of available rules = " + rule_num + "");
            Console.WriteLine("");
            Console.WriteLine("      Rule    Degree     Order");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                order_num = KeastRule.keast_order_num(rule);
                degree = KeastRule.keast_degree(rule);
                Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                       + "  " + degree.ToString().PadLeft(8)
                                       + "  " + order_num.ToString().PadLeft(8) + "");
            }

        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests KEAST_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int order;
            int order_num;
            int rule;
            int rule_num;
            double[] wtab;
            double wtab_sum;
            double[] xyztab;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  KEAST_RULE returns the points and weights");
            Console.WriteLine("  of a Keast rule for the triangle.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply check that the weights");
            Console.WriteLine("  sum to 1.");

            rule_num = KeastRule.keast_rule_num();

            Console.WriteLine("");
            Console.WriteLine("  Number of available rules = " + rule_num + "");

            Console.WriteLine("");
            Console.WriteLine("      Rule      Order     Sum of weights");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                order_num = KeastRule.keast_order_num(rule);

                xyztab = new double[3 * order_num];
                wtab = new double[order_num];

                KeastRule.keast_rule(rule, order_num, ref xyztab, ref wtab);

                wtab_sum = 0.0;
                for (order = 0; order < order_num; order++)
                {
                    wtab_sum = wtab_sum + wtab[order];
                }

                Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                       + "  " + order_num.ToString().PadLeft(8)
                                       + "  " + wtab_sum.ToString().PadLeft(14) + "");
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests KEAST_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int rule;
            int rule_num;
            int suborder;
            int suborder_num;
            double[] suborder_w;
            double[] suborder_xyzz;
            double xyzz_sum;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  KEAST_RULE returns the points and weights");
            Console.WriteLine("  of a Keast rule for the triangle.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply check that, for each");
            Console.WriteLine("  quadrature point, the barycentric coordinates");
            Console.WriteLine("  add up to 1.");

            rule_num = KeastRule.keast_rule_num();

            Console.WriteLine("");
            Console.WriteLine("      Rule    Suborder    Sum of coordinates");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                suborder_num = KeastRule.keast_suborder_num(rule);

                suborder_xyzz = new double[4 * suborder_num];
                suborder_w = new double[suborder_num];

                KeastRule.keast_subrule(rule, suborder_num, ref suborder_xyzz, ref suborder_w);

                Console.WriteLine("");
                Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                       + "  " + suborder_num.ToString().PadLeft(8) + "");

                for (suborder = 0; suborder < suborder_num; suborder++)
                {
                    xyzz_sum = suborder_xyzz[0 + suborder * 4]
                               + suborder_xyzz[1 + suborder * 4]
                               + suborder_xyzz[2 + suborder * 4]
                               + suborder_xyzz[3 + suborder * 4];
                    Console.WriteLine("                    "
                                      + "  " + xyzz_sum.ToString("0.################").PadLeft(25) + "");
                }
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 demonstrates TETRAHEDRON_REFERENCE_TO_PHYSICAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int NODE_NUM = 4;

            int node;
            double[] node_xyz =
            {
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0
            };
            double[] node_xyz2 =
            {
                1.0, 2.0, 3.0,
                2.0, 2.0, 3.0,
                1.0, 3.0, 3.0,
                1.0, 2.0, 9.0
            };
            int order;
            int order_num;
            int rule;
            double volume;
            double volume2;
            double[] w;
            double[] xyz;
            double[] xyz2;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  TETRAHEDRON_REFERENCE_TO_PHYSICAL transforms a rule");
            Console.WriteLine("  on the unit (reference) tetrahedron to a rule on ");
            Console.WriteLine("  an arbitrary (physical) tetrahedron.");

            rule = 2;

            order_num = KeastRule.keast_order_num(rule);

            xyz = new double[3 * order_num];
            xyz2 = new double[3 * order_num];
            w = new double[order_num];

            KeastRule.keast_rule(rule, order_num, ref xyz, ref w);
            //
            //  Here is the reference tetrahedron, and its rule.
            //
            Console.WriteLine("");
            Console.WriteLine("  The reference tetrahedron:");
            Console.WriteLine("");

            for (node = 0; node < NODE_NUM; node++)
            {
                Console.WriteLine("  " + (node + 1).ToString().PadLeft(8)
                                       + "  " + node_xyz[0 + node * 3].ToString().PadLeft(14)
                                       + "  " + node_xyz[1 + node * 3].ToString().PadLeft(14)
                                       + "  " + node_xyz[2 + node * 3].ToString().PadLeft(14) + "");
            }

            volume = Tetrahedron.tetrahedron_volume(node_xyz);

            Console.WriteLine("");
            Console.WriteLine("  Rule " + rule + " for reference tetrahedron");
            Console.WriteLine("  with volume = " + volume + "");
            Console.WriteLine("");
            Console.WriteLine("                X               Y               Z             W");
            Console.WriteLine("");

            for (order = 0; order < order_num; order++)
            {
                Console.WriteLine("  " + order.ToString().PadLeft(8)
                                       + "  " + xyz[0 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz[1 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz[2 + order * 3].ToString().PadLeft(14)
                                       + "  " + w[order].ToString().PadLeft(14) + "");
            }

            //
            //  Transform the rule.
            //
            Tetrahedron.reference_to_physical_t4(node_xyz2, order_num, xyz, ref xyz2);
            //
            //  Here is the physical tetrahedron, and its transformed rule.
            //
            Console.WriteLine("");
            Console.WriteLine("  The physical tetrahedron:");
            Console.WriteLine("");

            for (node = 0; node < NODE_NUM; node++)
            {
                Console.WriteLine("  " + (node + 1).ToString().PadLeft(8)
                                       + "  " + node_xyz2[0 + node * 3].ToString().PadLeft(14)
                                       + "  " + node_xyz2[1 + node * 3].ToString().PadLeft(14)
                                       + "  " + node_xyz2[2 + node * 3].ToString().PadLeft(14) + "");
            }

            volume2 = Tetrahedron.tetrahedron_volume(node_xyz2);

            Console.WriteLine("");
            Console.WriteLine("  Rule " + rule + " for physical tetrahedron");
            Console.WriteLine("  with volume = " + volume2 + "");
            Console.WriteLine("");
            Console.WriteLine("                X               Y               Z             W");
            Console.WriteLine("");

            for (order = 0; order < order_num; order++)
            {
                Console.WriteLine("  " + order.ToString().PadLeft(8)
                                       + "  " + xyz2[0 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz2[1 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz2[2 + order * 3].ToString().PadLeft(14)
                                       + "  " + w[order].ToString().PadLeft(14) + "");
            }
        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests KEAST_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree;
            int degree_max = 3;
            int dim_num = 3;
            int[] expon = new int[3];
            int h;
            double[] mono;
            bool more;
            int order_num;
            double quad;
            int rule;
            int rule_num;
            int t;
            double[] w;
            double[] xyz;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Demonstrate the KEAST rules on a sequence of");
            Console.WriteLine("  monomial integrands X^A Y^B Z^C");
            Console.WriteLine("  on the unit tetrahedron.");

            rule_num = KeastRule.keast_rule_num();

            Console.WriteLine("");
            Console.WriteLine("      Rule     Order     Quad");
            Console.WriteLine("");

            for (degree = 0; degree <= degree_max; degree++)
            {
                more = false;
                h = 0;
                t = 0;

                for (;;)
                {
                    Comp.comp_next(degree, dim_num, ref expon, ref more, ref h, ref t);

                    Console.WriteLine("");
                    Console.WriteLine("  F(X,Y,Z)"
                                      + " = X^" + expon[0]
                                      + " * Y^" + expon[1]
                                      + " * Z^" + expon[2] + "");
                    Console.WriteLine("");

                    for (rule = 1; rule <= rule_num; rule++)
                    {
                        order_num = KeastRule.keast_order_num(rule);

                        xyz = new double[3 * order_num];
                        w = new double[order_num];

                        KeastRule.keast_rule(rule, order_num, ref xyz, ref w);

                        mono = Monomial.monomial_value(dim_num, order_num, xyz, expon);

                        quad = typeMethods.r8vec_dot(order_num, w, mono);

                        Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                               + "  " + order_num.ToString().PadLeft(8)
                                               + "  " + quad.ToString().PadLeft(14) + "");
                    }

                    if (!more)
                    {
                        break;
                    }
                }
            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests KEAST_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree;
            int order;
            int order_num;
            int rule;
            double[] w;
            double[] xyz;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  KEAST_RULE returns the points and weights");
            Console.WriteLine("  of a Keast rule for the triangle.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply print a rule.");

            rule = 10;
            degree = KeastRule.keast_degree(rule);
            order_num = KeastRule.keast_order_num(rule);

            Console.WriteLine("");
            Console.WriteLine("  Rule =   " + rule + "");
            Console.WriteLine("  Degree = " + degree + "");
            Console.WriteLine("  Order =  " + order_num + "");

            Console.WriteLine("");
            Console.WriteLine("         I      W               X               Y               Z");
            Console.WriteLine("");

            xyz = new double[3 * order_num];
            w = new double[order_num];

            KeastRule.keast_rule(rule, order_num, ref xyz, ref w);

            for (order = 0; order < order_num; order++)
            {
                Console.WriteLine("  " + order.ToString().PadLeft(8)
                                       + "  " + w[order].ToString().PadLeft(14)
                                       + "  " + xyz[0 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz[1 + order * 3].ToString().PadLeft(14)
                                       + "  " + xyz[2 + order * 3].ToString().PadLeft(14) + "");
            }
        }
    }
}