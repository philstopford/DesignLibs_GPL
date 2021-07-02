using System;
using Burkardt;
using Burkardt.Types;
using QuadratureRule = Burkardt.TriangleNS.QuadratureRule;

namespace TriangleRulesTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGLE_FELIPPA_RULE_TEST.
            //
            //  Discussion:
            //
            //    TRIANGLE_FELIPPA_RULE_TEST tests the TRIANGLE_FELIPPA_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree_max;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_FELIPPA_RULE_TEST");
            Console.WriteLine("  Test the TRIANGLE_FELIPPA_RULE library.");

            degree_max = 4;
            triangle_unit_monomial_test(degree_max);

            degree_max = 7;
            triangle_unit_quad_test(degree_max);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_FELIPPA_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void triangle_unit_monomial_test(int degree_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_UNIT_MONOMIAL_TEST tests TRIANGLE_UNIT_MONOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 April 2009
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
            int[] expon = new int[2];
            double value;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_UNIT_MONOMIAL_TEST");
            Console.WriteLine("  For the unit triangle,");
            Console.WriteLine("  TRIANGLE_UNIT_MONOMIAL returns the exact value of the");
            Console.WriteLine("  integral of X^ALPHA Y^BETA");
            Console.WriteLine("");
            Console.WriteLine("  Volume = " + QuadratureRule.triangle_unit_volume() + "");
            Console.WriteLine("");
            Console.WriteLine("     ALPHA      BETA      INTEGRAL");
            Console.WriteLine("");

            for (alpha = 0; alpha <= degree_max; alpha++)
            {
                expon[0] = alpha;
                for (beta = 0; beta <= degree_max - alpha; beta++)
                {
                    expon[1] = beta;

                    value = QuadratureRule.triangle_unit_monomial(expon);

                    Console.WriteLine("  " + expon[0].ToString().PadLeft(8)
                        + "  " + expon[1].ToString().PadLeft(8)
                        + "  " + value.ToString().PadLeft(14) + "");
                }
            }
        }

        static void triangle_unit_quad_test(int degree_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_UNIT_QUAD_TEST tests the rules for the unit triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 April 2008
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
            int DIM_NUM = 2;

            int dim;
            int dim_num = DIM_NUM;
            int[] expon = new int[DIM_NUM];
            int h = 0;
            bool more;
            int order;
            double quad;
            int t = 0;
            double[] v;
            double[] w;
            double[] xy;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_UNIT_QUAD_TEST");
            Console.WriteLine("  For the unit triangle,");
            Console.WriteLine("  we approximate monomial integrals with:");
            Console.WriteLine("  TRIANGLE_UNIT_O01,");
            Console.WriteLine("  TRIANGLE_UNIT_O03,");
            Console.WriteLine("  TRIANGLE_UNIT_O03b,");
            Console.WriteLine("  TRIANGLE_UNIT_O06,");
            Console.WriteLine("  TRIANGLE_UNIT_O06b,");
            Console.WriteLine("  TRIANGLE_UNIT_O07,");
            Console.WriteLine("  TRIANGLE_UNIT_O12,");

            more = false;

            SubCompData data = new SubCompData();
            
            for (;;)
            {
                SubComp.subcomp_next(ref data, degree_max, dim_num, ref expon, ref more, ref h, ref t);

                Console.WriteLine("");
                string cout = "  Monomial exponents: ";
                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + expon[dim].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                order = 1;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o01(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 3;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o03(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 3;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o03b(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 6;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o06(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 6;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o06b(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 7;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o07(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                order = 12;
                w = new double[order];
                xy = new double[dim_num * order];
                QuadratureRule.triangle_unit_o12(ref w, ref xy);
                v = Monomial.monomial_value(dim_num, order, expon, xy);
                quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order.ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");

                Console.WriteLine("");
                quad = QuadratureRule.triangle_unit_monomial(expon);
                Console.WriteLine("  " + " Exact"
                    + "  " + quad.ToString().PadLeft(14) + "");

                if (!more)
                {
                    break;
                }
            }
        }
    }
}