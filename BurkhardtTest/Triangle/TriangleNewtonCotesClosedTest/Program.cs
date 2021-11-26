using System;
using System.Globalization;
using Burkardt.FEM;
using Burkardt.TriangleNS;

namespace TriangleNewtonCotesClosedTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_NCC_RULE_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_NCC_RULE_TEST tests the TRIANGLE_NCC_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_NCC_RULE_TEST:");
        Console.WriteLine("  Test the TRIANGLE_NCC_RULE library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_NCC_RULE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests TRIANGLE_NCC_RULE_NUM, TRIANGLE_NCC_DEGREE, TRIANGLE_NCC_ORDER_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  TRIANGLE_NCC_RULE_NUM returns the number of rules;");
        Console.WriteLine("  TRIANGLE_NCC_DEGREE returns the degree of a rule;");
        Console.WriteLine("  TRIANGLE_NCC_ORDER_NUM returns the order of a rule.");

        int rule_num = NewtonCotesClosed.triangle_ncc_rule_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of available rules = " + rule_num + "");
        Console.WriteLine("");
        Console.WriteLine("      Rule    Degree     Order");
        Console.WriteLine("");

        for (rule = 1; rule <= rule_num; rule++)
        {
            int order_num = NewtonCotesClosed.triangle_ncc_order_num(rule);
            int degree = NewtonCotesClosed.triangle_ncc_degree(rule);
            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + degree.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + order_num.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests TRIANGLE_NCC_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  TRIANGLE_NCC_RULE returns the points and weights");
        Console.WriteLine("  of an NCC rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply check that the weights");
        Console.WriteLine("  sum to 1.");

        int rule_num = NewtonCotesClosed.triangle_ncc_rule_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of available rules = " + rule_num + "");

        Console.WriteLine("");
        Console.WriteLine("      Rule      Order     Sum of weights");
        Console.WriteLine("");

        for (rule = 1; rule <= rule_num; rule++)
        {
            int order_num = NewtonCotesClosed.triangle_ncc_order_num(rule);

            double[] xytab = new double[2 * order_num];
            double[] wtab = new double[order_num];

            NewtonCotesClosed.triangle_ncc_rule(rule, order_num, ref xytab, ref wtab);

            double wtab_sum = 0.0;
            int order;
            for (order = 0; order < order_num; order++)
            {
                wtab_sum += wtab[order];
            }

            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + order_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + wtab_sum.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests TRIANGLE_NCC_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  TRIANGLE_NCC_RULE returns the points and weights");
        Console.WriteLine("  of an NCC rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply check that, for each");
        Console.WriteLine("  quadrature point, the barycentric coordinates");
        Console.WriteLine("  add up to 1.");

        int rule_num = NewtonCotesClosed.triangle_ncc_rule_num();

        Console.WriteLine("");
        Console.WriteLine("      Rule    Suborder    Sum of coordinates");
        Console.WriteLine("");

        for (rule = 1; rule <= rule_num; rule++)
        {
            int suborder_num = NewtonCotesClosed.triangle_ncc_suborder_num(rule);

            double[] suborder_xyz = new double[3 * suborder_num];
            double[] suborder_w = new double[suborder_num];

            NewtonCotesClosed.triangle_ncc_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);

            Console.WriteLine("");
            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + suborder_num.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            int suborder;
            for (suborder = 0; suborder < suborder_num; suborder++)
            {
                double xyz_sum = suborder_xyz[0 + suborder * 3]
                                 + suborder_xyz[1 + suborder * 3]
                                 + suborder_xyz[2 + suborder * 3];
                Console.WriteLine("                    "
                                  + "  " + xyz_sum.ToString("0.################").PadLeft(25) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests TRIANGLE_NCC_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a;
        const double area = 0.5;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  TRIANGLE_NCC_RULE returns the points and weights of");
        Console.WriteLine("  an NCC rule for the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  This routine uses those rules to estimate the");
        Console.WriteLine("  integral of monomomials in the unit triangle.");

        int rule_num = NewtonCotesClosed.triangle_ncc_rule_num();

        for (a = 0; a <= 10; a++)
        {
            int b;
            for (b = 0; b <= 10 - a; b++)
            {
                //
                //  Multiplying X^A * Y^B by COEF will give us an integrand
                //  whose integral is exactly 1.  This makes the error calculations easy.
                //
                double coef = (a + b + 2) * (double) (a + b + 1);
                int i;
                for (i = 1; i <= b; i++)
                {
                    coef = coef * (a + i) / i;
                }

                Console.WriteLine("");
                Console.WriteLine("  Integrate " + coef
                                                 + " * X^" + a
                                                 + " * Y^" + b + "");
                Console.WriteLine("");
                Console.WriteLine("      Rule       QUAD           ERROR");
                Console.WriteLine("");

                int rule;
                for (rule = 1; rule <= rule_num; rule++)
                {
                    int order_num = NewtonCotesClosed.triangle_ncc_order_num(rule);

                    double[] xytab = new double[2 * order_num];
                    double[] wtab = new double[order_num];

                    NewtonCotesClosed.triangle_ncc_rule(rule, order_num, ref xytab, ref wtab);

                    double quad = 0.0;

                    int order;
                    for (order = 0; order < order_num; order++)
                    {
                        double x = xytab[0 + order * 2];
                        double y = xytab[1 + order * 2];

                        switch (a)
                        {
                            case 0 when b == 0:
                                value = coef;
                                break;
                            case 0 when b != 0:
                                value = coef * Math.Pow(y, b);
                                break;
                            default:
                            {
                                if (a != 0 && b == 0)
                                {
                                    value = coef * Math.Pow(x, a);
                                }
                                else if (a != 0 && b != 0)
                                {
                                    value = coef * Math.Pow(x, a) * Math.Pow(y, b);
                                }

                                break;
                            }
                        }

                        quad += wtab[order] * value;
                    }

                    quad = area * quad;

                    double exact = 1.0;
                    double err = Math.Abs(exact - quad);

                    Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 demonstrates REFERENCE_TO_PHYSICAL_T3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NODE_NUM = 3;

        int node;
        double[] node_xy =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0
        };
        double[] node_xy2 =
        {
            1.0, 2.0,
            1.0, 1.0,
            3.0, 2.0
        };
        int order;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  REFERENCE_TO_PHYSICAL_T3 transforms a rule");
        Console.WriteLine("  on the unit (reference) triangle to a rule on ");
        Console.WriteLine("  an arbitrary (physical) triangle.");

        const int rule = 3;

        int order_num = NewtonCotesClosed.triangle_ncc_order_num(rule);

        double[] xy = new double[2 * order_num];
        double[] xy2 = new double[2 * order_num];
        double[] w = new double[order_num];

        NewtonCotesClosed.triangle_ncc_rule(rule, order_num, ref xy, ref w);
        //
        //  Here is the reference triangle, and its rule.
        //
        Console.WriteLine("");
        Console.WriteLine("  The reference triangle:");
        Console.WriteLine("");

        for (node = 0; node < NODE_NUM; node++)
        {
            Console.WriteLine("  " + node + 1
                              + "  " + node_xy[0 + node * 2]
                              + "  " + node_xy[1 + node * 2] + "");
        }

        double area = Integrals.triangle_area(node_xy);

        Console.WriteLine("");
        Console.WriteLine("  Rule " + rule + " for reference triangle");
        Console.WriteLine("  with area = " + area + "");
        Console.WriteLine("");
        Console.WriteLine("                X               Y               W");
        Console.WriteLine("");

        for (order = 0; order < order_num; order++)
        {
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xy[0 + order * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + xy[1 + order * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + w[order].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Transform the rule.
        //
        Reference.reference_to_physical_t3(node_xy2, order_num, xy, ref xy2);
        //
        //  Here is the physical triangle, and its transformed rule.
        //
        Console.WriteLine("");
        Console.WriteLine("  The physical triangle:");
        Console.WriteLine("");

        for (node = 0; node < NODE_NUM; node++)
        {
            Console.WriteLine("  " + (node + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + node_xy2[0 + node * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + node_xy2[1 + node * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        double area2 = Integrals.triangle_area(node_xy2);

        Console.WriteLine("");
        Console.WriteLine("  Rule " + rule + " for physical triangle");
        Console.WriteLine("  with area = " + area2 + "");
        Console.WriteLine("");
        Console.WriteLine("                X               Y               W");
        Console.WriteLine("");

        for (order = 0; order < order_num; order++)
        {
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xy2[0 + order * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + xy2[1 + order * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + w[order].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}