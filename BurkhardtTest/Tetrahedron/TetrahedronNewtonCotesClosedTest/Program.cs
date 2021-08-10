using System;
using Burkardt.TetrahedronNS;

namespace TetrahedronNewtonCotesClosedTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TETRAHEDRON_NCC_RULE_TEST.
            //
            //  Discussion:
            //
            //    TETRAHEDRON_NCC_RULE_TEST tests the TETRAHEDRON_NCC_RULE library.
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
            Console.WriteLine("TETRAHEDRON_NCC_RULE_TEST:");
            Console.WriteLine("  Test the TETRAHEDRON_NCC_RULE library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();

            Console.WriteLine("");
            Console.WriteLine("TETRAHEDRON_NCC_RULE_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests TETRAHEDRON_NCC_RULE_NUM, TETRAHEDRON_NCC_DEGREE, TETRAHEDRON_NCC_ORDER_NUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 January 2007
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
            Console.WriteLine("  TETRAHEDRON_NCC_RULE_NUM returns the number of rules;");
            Console.WriteLine("  TETRAHEDRON_NCC_DEGREE returns the degree of a rule;");
            Console.WriteLine("  TETRAHEDRON_NCC_ORDER_NUM returns the order of a rule.");

            rule_num = NewtonCotesClosed.tetrahedron_ncc_rule_num();

            Console.WriteLine("");
            Console.WriteLine("  Number of available rules = " + rule_num + "");
            Console.WriteLine("");
            Console.WriteLine("      Rule    Degree     Order");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                order_num = NewtonCotesClosed.tetrahedron_ncc_order_num(rule);
                degree = NewtonCotesClosed.tetrahedron_ncc_degree(rule);
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
            //    TEST02 tests TETRAHEDRON_NCC_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num = 3;
            int order;
            int order_num;
            int rule;
            int rule_num;
            double[] wtab;
            double wtab_sum;
            double[] xyztab;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  TETRAHEDRON_NCC_RULE returns the points and weights");
            Console.WriteLine("  of an NCC rule for the tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply check that the weights");
            Console.WriteLine("  sum to 1.");

            rule_num = NewtonCotesClosed.tetrahedron_ncc_rule_num();

            Console.WriteLine("");
            Console.WriteLine("  Number of available rules = " + rule_num + "");

            Console.WriteLine("");
            Console.WriteLine("      Rule      Sum of weights");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                order_num = NewtonCotesClosed.tetrahedron_ncc_order_num(rule);

                xyztab = new double[dim_num * order_num];
                wtab = new double[order_num];

                NewtonCotesClosed.tetrahedron_ncc_rule(rule, order_num, ref xyztab, ref wtab);

                wtab_sum = 0.0;
                for (order = 0; order < order_num; order++)
                {
                    wtab_sum = wtab_sum + wtab[order];
                }

                Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                       + "  " + wtab_sum.ToString().PadLeft(14) + "");
            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests TETRAHEDRON_NCC_RULE.
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
            int rule_num;
            int suborder;
            int suborder_num;
            double[] suborder_w;
            double[] suborder_xyz;
            double xyz_sum;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  TETRAHEDRON_NCC_RULE returns the points and weights");
            Console.WriteLine("  of an NCC rule for the tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply check that, for each");
            Console.WriteLine("  quadrature point, the barycentric coordinates");
            Console.WriteLine("  add up to 1.");

            rule_num = NewtonCotesClosed.tetrahedron_ncc_rule_num();

            Console.WriteLine("");
            Console.WriteLine("      Rule    Suborder    Sum of coordinates");
            Console.WriteLine("");

            for (rule = 1; rule <= rule_num; rule++)
            {
                suborder_num = NewtonCotesClosed.tetrahedron_ncc_suborder_num(rule);

                suborder_xyz = new double[4 * suborder_num];
                suborder_w = new double[suborder_num];

                NewtonCotesClosed.tetrahedron_ncc_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);

                Console.WriteLine("");
                Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                       + "  " + suborder_num.ToString().PadLeft(8) + "");

                for (suborder = 0; suborder < suborder_num; suborder++)
                {
                    xyz_sum = suborder_xyz[0 + suborder * 4]
                              + suborder_xyz[1 + suborder * 4]
                              + suborder_xyz[2 + suborder * 4]
                              + suborder_xyz[3 + suborder * 4];
                    Console.WriteLine("                    "
                                      + "  " + xyz_sum.ToString("0.################").PadLeft(25) + "");
                }
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests TETRAHEDRON_NCC_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int a;
            int b;
            int c;
            double coef;
            int dim_num = 3;
            double err;
            double exact;
            int i;
            int order;
            int order_num;
            double quad;
            int rule;
            int rule_num;
            double value;
            double volume;
            double[] wtab;
            double x;
            double[] xyztab;
            double y;
            double z;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  TETRAHEDRON_NCC_RULE returns the points and weights of");
            Console.WriteLine("  an NCC rule for the unit tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  This routine uses those rules to estimate the");
            Console.WriteLine("  integral of monomomials in the unit tetrahedron.");

            rule_num = NewtonCotesClosed.tetrahedron_ncc_rule_num();

            volume = 1.0 / 6.0;

            for (a = 0; a <= 6; a++)
            {
                for (b = 0; b <= 6 - a; b++)
                {
                    for (c = 0; c <= 6 - a - b; c++)
                    {
                        //
                        //  Multiplying X**A * Y**B * Z**C by COEF will give us an integrand
                        //  whose integral is exactly 1.  This makes the error calculations easy.
                        //
                        coef = 1.0;

                        //      for ( i = 1; i <= a; i++ )
                        //      {
                        //        coef = coef * i / i;
                        //      }
                        for (i = a + 1; i <= a + b; i++)
                        {
                            coef = coef * (double)(i) / (double)(i - a);
                        }

                        for (i = a + b + 1; i <= a + b + c; i++)
                        {
                            coef = coef * (double)(i) / (double)(i - a - b);
                        }

                        for (i = a + b + c + 1; i <= a + b + c + 3; i++)
                        {
                            coef = coef * (double)(i);
                        }

                        Console.WriteLine("");
                        Console.WriteLine("  Integrate " + coef
                                                         + " * X^" + a
                                                         + " * Y^" + b
                                                         + " * Z^" + c + "");
                        Console.WriteLine("");
                        Console.WriteLine("      Rule       QUAD           ERROR");
                        Console.WriteLine("");

                        for (rule = 1; rule <= rule_num; rule++)
                        {
                            order_num = NewtonCotesClosed.tetrahedron_ncc_order_num(rule);

                            xyztab = new double[dim_num * order_num];
                            wtab = new double[order_num];

                            NewtonCotesClosed.tetrahedron_ncc_rule(rule, order_num, ref xyztab, ref wtab);

                            quad = 0.0;

                            for (order = 0; order < order_num; order++)
                            {
                                x = xyztab[0 + order * dim_num];
                                y = xyztab[1 + order * dim_num];
                                z = xyztab[2 + order * dim_num];
                                //
                                //  Some tedious calculations to avoid 0**0 complaints.
                                //
                                value = coef;

                                if (a != 0)
                                {
                                    value = value * Math.Pow(x, a);
                                }

                                if (b != 0)
                                {
                                    value = value * Math.Pow(y, b);
                                }

                                if (c != 0)
                                {
                                    value = value * Math.Pow(z, c);
                                }

                                quad = quad + wtab[order] * value;
                            }

                            quad = volume * quad;

                            exact = 1.0;
                            err = Math.Abs(exact - quad);

                            Console.WriteLine("  " + rule.ToString().PadLeft(8)
                                                   + "  " + quad.ToString().PadLeft(14)
                                                   + "  " + err.ToString().PadLeft(14) + "");

                        }
                    }
                }
            }
        }

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 demonstrates REFERENCE_TO_PHYSICAL_T4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
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
                4.0, 5.0, 1.0,
                6.0, 5.0, 1.0,
                4.0, 8.0, 1.0,
                4.0, 5.0, 5.0
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
            Console.WriteLine("TEST06");
            Console.WriteLine("  REFERENCE_TO_PHYSICAL_T4 transforms a rule");
            Console.WriteLine("  on the unit (reference) tetrahedron to a rule on ");
            Console.WriteLine("  an arbitrary (physical) tetrahedron.");

            rule = 3;

            order_num = NewtonCotesClosed.tetrahedron_ncc_order_num(rule);

            xyz = new double[DIM_NUM * order_num];
            xyz2 = new double[DIM_NUM * order_num];
            w = new double[order_num];

            NewtonCotesClosed.tetrahedron_ncc_rule(rule, order_num, ref xyz, ref w);
            //
            //  Here is the reference tetrahedron, and its rule.
            //
            Console.WriteLine("");
            Console.WriteLine("  The reference tetrahedron:");
            Console.WriteLine("");

            for (node = 0; node < NODE_NUM; node++)
            {
                Console.WriteLine("  " + (node + 1).ToString().PadLeft(8)
                                       + "  " + node_xyz[0 + node * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + node_xyz[1 + node * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + node_xyz[2 + node * DIM_NUM].ToString().PadLeft(14) + "");
            }

            volume = Tetrahedron.tetrahedron_volume(node_xyz);

            Console.WriteLine("");
            Console.WriteLine("  Rule " + rule + " for reference tetrahedron");
            Console.WriteLine("  with volume = " + volume + "");
            Console.WriteLine("");
            Console.WriteLine("                X               Y               Z               W");
            Console.WriteLine("");

            for (order = 0; order < order_num; order++)
            {
                Console.WriteLine("  " + order.ToString().PadLeft(8)
                                       + "  " + xyz[0 + order * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + xyz[1 + order * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + xyz[2 + order * DIM_NUM].ToString().PadLeft(14)
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
                                       + "  " + node_xyz2[0 + node * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + node_xyz2[1 + node * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + node_xyz2[2 + node * DIM_NUM].ToString().PadLeft(14) + "");
            }

            volume2 = Tetrahedron.tetrahedron_volume(node_xyz2);

            Console.WriteLine("");
            Console.WriteLine("  Rule " + rule + " for physical tetrahedron");
            Console.WriteLine("  with volume = " + volume2 + "");
            Console.WriteLine("");
            Console.WriteLine("                X               Y               Z               W");
            Console.WriteLine("");

            for (order = 0; order < order_num; order++)
            {
                Console.WriteLine("  " + order.ToString().PadLeft(8)
                                       + "  " + xyz2[0 + order * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + xyz2[1 + order * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + xyz2[2 + order * DIM_NUM].ToString().PadLeft(14)
                                       + "  " + w[order].ToString().PadLeft(14) + "");
            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests TETRAHEDRON_NCC_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num = 3;
            int o;
            int order_num;
            int rule;
            int s;
            int[] suborder;
            int suborder_num;
            double[] suborder_w;
            double[] suborder_xyz;
            double[] w;
            double[] xyz;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  TETRAHEDRON_NCC_RULE returns the points and weights");
            Console.WriteLine("  of an NCC rule for the tetrahedron.");

            rule = 4;

            Console.WriteLine("");
            Console.WriteLine("  In this test, we simply print rule " + rule + "");

            suborder_num = NewtonCotesClosed.tetrahedron_ncc_suborder_num(rule);

            suborder = NewtonCotesClosed.tetrahedron_ncc_suborder(rule, suborder_num);

            suborder_w = new double[suborder_num];
            suborder_xyz = new double[4 * suborder_num];

            NewtonCotesClosed.tetrahedron_ncc_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);

            Console.WriteLine("");
            Console.WriteLine("  The compressed rule:");
            Console.WriteLine("");
            Console.WriteLine("  Number of suborders = " + suborder_num + "");
            Console.WriteLine("");
            Console.WriteLine("     S   Sub     Weight     Xsi1      Xsi2      Xsi3      Xsi4");
            Console.WriteLine("");

            for (s = 0; s < suborder_num; s++)
            {
                Console.WriteLine("  " + (s + 1).ToString().PadLeft(4)
                                       + "  " + suborder[s].ToString().PadLeft(4)
                                       + "  " + suborder_w[s].ToString().PadLeft(8)
                                       + "  " + suborder_xyz[0 + s * 4].ToString().PadLeft(8)
                                       + "  " + suborder_xyz[1 + s * 4].ToString().PadLeft(8)
                                       + "  " + suborder_xyz[2 + s * 4].ToString().PadLeft(8)
                                       + "  " + suborder_xyz[3 + s * 4].ToString("0.####").PadLeft(8) + "");
            }

            order_num = NewtonCotesClosed.tetrahedron_ncc_order_num(rule);

            xyz = new double[dim_num * order_num];
            w = new double[order_num];

            NewtonCotesClosed.tetrahedron_ncc_rule(rule, order_num, ref xyz, ref w);

            Console.WriteLine("");
            Console.WriteLine("  The full rule:");
            Console.WriteLine("");
            Console.WriteLine("  Order = " + order_num + "");
            Console.WriteLine("");
            Console.WriteLine("     O    Weight        X           Y           Z");
            Console.WriteLine("");

            for (o = 0; o < order_num; o++)
            {
                Console.WriteLine("  " + (o + 1).ToString().PadLeft(4)
                                       + "  " + w[o].ToString().PadLeft(8)
                                       + "  " + xyz[0 + o * 3].ToString().PadLeft(8)
                                       + "  " + xyz[1 + o * 3].ToString().PadLeft(8)
                                       + "  " + xyz[2 + o * 3].ToString("0.####").PadLeft(8) + "");
            }
        }
    }
}