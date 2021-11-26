using System;
using System.Globalization;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace TriangleLynessRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_LYNESS_RULE_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_LYNESS_RULE_TEST tests the TRIANGLE_LYNESS_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_LYNESS_RULE_TEST:");
        Console.WriteLine("  Test the TRIANGLE_LYNESS_RULE library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_LYNESS_RULE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LYNESS_RULE_NUM, LYNESS_DEGREE, LYNESS_ORDER_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  LYNESS_RULE_NUM returns the number of rules;");
        Console.WriteLine("  LYNESS_DEGREE returns the degree of a rule;");
        Console.WriteLine("  LYNESS_ORDER_NUM returns the order of a rule.");

        int rule_num = LynessRule.lyness_rule_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of available rules = " + rule_num + "");
        Console.WriteLine("");
        Console.WriteLine("      Rule     Order  Precision");
        Console.WriteLine("");

        for (rule = 0; rule <= rule_num; rule++)
        {
            int order = LynessRule.lyness_order(rule);
            int precision = LynessRule.lyness_precision(rule);
            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + precision.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 performs the weight sum test on Lyness rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  LYNESS_RULE returns the points and weights");
        Console.WriteLine("  of a Lyness rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply check that the weights");
        Console.WriteLine("  sum to 1.");

        int rule_num = LynessRule.lyness_rule_num();

        Console.WriteLine("");
        Console.WriteLine("  Number of available rules = " + rule_num + "");
        Console.WriteLine("");
        Console.WriteLine("      Rule    Sum of weights");
        Console.WriteLine("");

        for (rule = 0; rule <= rule_num; rule++)
        {
            int order = LynessRule.lyness_order(rule);

            double[] x = new double[2 * order];
            double[] w = new double[order];

            LynessRule.lyness_rule(rule, order, ref w, ref x);

            double w_sum = typeMethods.r8vec_sum(order, w);

            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + w_sum.ToString(CultureInfo.InvariantCulture).PadLeft(25) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 performs the barycentric coordinate sum test on Lyness rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  LYNESS_RULE returns the points and weights");
        Console.WriteLine("  of a Lyness rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply check that, for each");
        Console.WriteLine("  quadrature point, the barycentric coordinates");
        Console.WriteLine("  sum to 1.");

        int rule_num = LynessRule.lyness_rule_num();

        Console.WriteLine("");
        Console.WriteLine("      Rule   Suborder    Sum of coordinates");

        for (rule = 0; rule < rule_num; rule++)
        {
            int suborder_num = LynessRule.lyness_suborder_num(rule);

            double[] sub_w = new double[suborder_num];
            double[] sub_xyz = new double[3 * suborder_num];

            LynessRule.lyness_subrule(rule, suborder_num, ref sub_xyz, ref sub_w);

            Console.WriteLine("");
            Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + suborder_num.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            int suborder;
            for (suborder = 0; suborder < suborder_num; suborder++)
            {
                double xyz_sum = typeMethods.r8vec_sum(3, sub_xyz, +3 * suborder);
                Console.WriteLine("                           " + xyz_sum + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 prints a rule generated by LYNESS_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  LYNESS_RULE returns the points and weights");
        Console.WriteLine("  of a Lyness rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply print a rule.");

        int rule = 18;
        int order = LynessRule.lyness_order(rule);
        int precision = LynessRule.lyness_precision(rule);

        Console.WriteLine("");
        Console.WriteLine("  Rule =      " + rule + "");
        Console.WriteLine("  Order =     " + order + "");
        Console.WriteLine("  Precision = " + precision + "");

        double[] w = new double[order];
        double[] x = new double[2 * order];

        LynessRule.lyness_rule(rule, order, ref w, ref x);

        Console.WriteLine("");
        Console.WriteLine("     I      W                   X                   Y");
        Console.WriteLine("");

        for (j = 0; j < order; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + w[j].ToString("0.################").PadLeft(24)
                                   + "  " + x[0 + j * 2].ToString("0.################").PadLeft(24)
                                   + "  " + x[1 + j * 2].ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 writes a rule created by LYNESS_RULE to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] r =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  LYNESS_RULE returns the points and weights");
        Console.WriteLine("  of a Lyness rule for the triangle.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, we simply print a rule.");

        const int rule = 18;
        int order = LynessRule.lyness_order(rule);
        int precision = LynessRule.lyness_precision(rule);

        Console.WriteLine("");
        Console.WriteLine("  Rule =      " + rule + "");
        Console.WriteLine("  Order =     " + order + "");
        Console.WriteLine("  Precision = " + precision + "");

        double[] w = new double[order];
        double[] x = new double[2 * order];

        LynessRule.lyness_rule(rule, order, ref w, ref x);

        typeMethods.r8mat_write("lyness_18_r.txt", 2, 3, r);
        Console.WriteLine("");
        Console.WriteLine("  Wrote the region file \"lyness_18_r.txt\".");
        typeMethods.r8mat_write("lyness_18_w.txt", 1, order, w);
        Console.WriteLine("  Wrote the weight file \"lyness_18_w.txt\".");
        typeMethods.r8mat_write("lyness_18_x.txt", 2, order, x);
        Console.WriteLine("  Wrote the point file \"lyness_18_x.txt\".");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests the Lyness rules for exact integration of monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree;
        const int degree_max = 10;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  LYNESS_RULE returns the points and weights of");
        Console.WriteLine("  a Lyness rule for the unit triangle.");
        Console.WriteLine("");
        Console.WriteLine("  This routine uses those rules to estimate the");
        Console.WriteLine("  integral of monomomials in the unit triangle.");

        int rule_num = LynessRule.lyness_rule_num();

        double area = 0.5;

        for (degree = 0; degree <= degree_max; degree++)
        {
            int a;
            for (a = 0; a <= degree; a++)
            {
                int b = degree - a;
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
                Console.WriteLine("  Integrate " + coef + " * X^" + a + " * Y^" + b + "");
                Console.WriteLine("");
                Console.WriteLine("      Rule       QUAD           ERROR");
                Console.WriteLine("");

                int rule;
                for (rule = 0; rule <= rule_num; rule++)
                {
                    int order = LynessRule.lyness_order(rule);

                    double[] w = new double[order];
                    double[] x = new double[2 * order];

                    LynessRule.lyness_rule(rule, order, ref w, ref x);

                    double quad = 0.0;

                    int j;
                    for (j = 0; j < order; j++)
                    {
                        switch (a)
                        {
                            case 0 when b == 0:
                                value = coef;
                                break;
                            case 0 when b != 0:
                                value = coef * Math.Pow(x[1 + j * 2], b);
                                break;
                            default:
                            {
                                if (a != 0 && b == 0)
                                {
                                    value = coef * Math.Pow(x[0 + j * 2], a);
                                }
                                else if (a != 0 && b != 0)
                                {
                                    value = coef * Math.Pow(x[0 + j * 2], a) * Math.Pow(x[1 + j * 2], b);
                                }

                                break;
                            }
                        }

                        quad += w[j] * value;
                    }

                    quad = area * quad;

                    double exact = 1.0;
                    double err = Math.Abs(exact - quad);

                    Console.WriteLine("  " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + quad.ToString("0.######").PadLeft(14)
                                           + "  " + err.ToString("0.##").PadLeft(10) + "");

                }
            }
        }
    }
}