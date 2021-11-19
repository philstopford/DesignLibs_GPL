using System;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.Types;
using QuadratureRule = Burkardt.Wedge.QuadratureRule;

namespace WedgeRulesTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WEDGE_FELIPPA_RULE_TEST.
        //
        //  Discussion:
        //
        //    FELIPPA_TEST tests the FELIPPA library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int degree_max;

        Console.WriteLine("");
        Console.WriteLine("WEDGE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Test the WEDGE_FELIPPA_RULE library.");

        degree_max = 4;

        test01(degree_max);
        test02(degree_max);

        test03();

        Console.WriteLine("");
        Console.WriteLine("WEDGE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests WEDGE_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        int alpha;
        int beta;
        int[] expon = new int[3];
        int gamma;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  For the unit wedge,");
        Console.WriteLine("  WEDGE_INTEGRAL returns the exact value of the");
        Console.WriteLine("  integral of X^ALPHA Y^BETA Z^GAMMA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + QuadratureRule.wedge_volume() + "");
        Console.WriteLine("");
        Console.WriteLine("   ALPHA      BETA     GAMMA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            expon[0] = alpha;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                expon[1] = beta;
                for (gamma = 0; gamma <= degree_max - alpha - beta; gamma++)
                {
                    expon[2] = gamma;
                    value = QuadratureRule.wedge_integral(expon);
                    Console.WriteLine(expon[0].ToString().PadLeft(8) + "  "
                                                                     + expon[1].ToString().PadLeft(8) + "  "
                                                                     + expon[2].ToString().PadLeft(8) + "  "
                                                                     + value.ToString().PadLeft(14) + "");
                }
            }
        }
    }

    private static void test02(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests the rules for the unit wedge.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        int dim_num = 3;
        int[] expon = new int[3];
        int h = 0;
        int line_order;
        int[] line_order_array =  {
                1, 2, 2, 3, 2, 3, 4
            }
            ;
        bool more;
        int order;
        double quad;
        int t = 0;
        int test;
        int test_num = 7;
        int triangle_order;
        int[] triangle_order_array =  {
                1, 3, -3, 6, -6, 7, 12
            }
            ;
        double[] v;
        double[] w;
        double[] xyz;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For the unit wedge,");
        Console.WriteLine("  we approximate monomial integrals with WEDG_UNIT_RULE.");

        more = false;

        SubCompData data = new();
            
        for (;;)
        {
            SubComp.subcomp_next(ref data, degree_max, dim_num, ref expon, ref more, ref h, ref t);

            if (expon[2] % 2 == 1)
            {
                if (!more)
                {
                    break;
                }

                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("  Monomial exponents:   "
                              + expon[0].ToString().PadLeft(2) + "  "
                              + expon[1].ToString().PadLeft(2) + "  "
                              + expon[2].ToString().PadLeft(2) + "");
            Console.WriteLine("");

            for (test = 0; test < test_num; test++)
            {
                line_order = line_order_array[test];
                triangle_order = triangle_order_array[test];

                order = line_order * Math.Abs(triangle_order);

                w = new double[order];
                xyz = new double[dim_num * order];
                QuadratureRule.wedge_rule(line_order, triangle_order, ref w, ref xyz);
                v = Monomial.monomial_value(dim_num, order, expon, xyz);
                quad = QuadratureRule.wedge_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine(triangle_order.ToString().PadLeft(6) + "  "
                                                                       + line_order.ToString().PadLeft(6) + "  "
                                                                       + quad.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            quad = QuadratureRule.wedge_integral(expon);
            Console.WriteLine("   Exact        "
                              + quad.ToString().PadLeft(14) + "");

            if (!more)
            {
                break;
            }

        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 writes out some rules for the unit wedge.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
    {
        int dim_num = 3;
        int line_order;
        int[] line_order_array =  {
                1, 2, 2, 3, 2, 3, 4
            }
            ;
        int order;
        int rule;
        int rule_num = 7;
        int triangle_order;
        int[] triangle_order_array =  {
                1, 3, -3, 6, -6, 7, 12
            }
            ;
        double[] w;
        string w_filename = "";
        double[] x;
        string x_filename = "";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For the unit wedge,");
        Console.WriteLine("  write some rules to a file");
        Console.WriteLine("");
        Console.WriteLine("   Rule  Trig    Line   Total  W_File X_File");
        Console.WriteLine("         Order   Order  Order");
        Console.WriteLine("");

        for (rule = 0; rule < rule_num; rule++)
        {
            switch (rule)
            {
                case 0:
                    w_filename = "wedge_felippa_1x1_w.txt";
                    x_filename = "wedge_felippa_1x1_x.txt";
                    break;
                case 1:
                    w_filename = "wedge_felippa_3x2_w.txt";
                    x_filename = "wedge_felippa_3x2_x.txt";
                    break;
                case 2:
                    w_filename = "wedge_felippa_3bx2_w.txt";
                    x_filename = "wedge_felippa_3bx2_x.txt";
                    break;
                case 3:
                    w_filename = "wedge_felippa_6x3_w.txt";
                    x_filename = "wedge_felippa_6x3_x.txt";
                    break;
                case 4:
                    w_filename = "wedge_felippa_6bx2_w.txt";
                    x_filename = "wedge_felippa_6bx2_x.txt";
                    break;
                case 5:
                    w_filename = "wedge_felippa_7x3_w.txt";
                    x_filename = "wedge_felippa_7x3_x.txt";
                    break;
                case 6:
                    w_filename = "wedge_felippa_12x4_w.txt";
                    x_filename = "wedge_felippa_12x4_x.txt";
                    break;
            }

            line_order = line_order_array[rule];
            triangle_order = triangle_order_array[rule];

            order = line_order * Math.Abs(triangle_order);

            w = new double[order];
            x = new double[dim_num * order];
            QuadratureRule.wedge_rule(line_order, triangle_order, ref w, ref x);
            typeMethods.r8mat_write(w_filename, 1, order, w);
            typeMethods.r8mat_write(x_filename, dim_num, order, x);
            Console.WriteLine(rule.ToString().PadLeft(6) + "  "
                                                         + triangle_order.ToString().PadLeft(6) + "  "
                                                         + line_order.ToString().PadLeft(6) + "  "
                                                         + order.ToString().PadLeft(6) + "  "
                                                         + w_filename + "  "
                                                         + x_filename + "");

        }
    }
}