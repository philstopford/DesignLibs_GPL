using System;
using Burkardt.Square;

namespace SquareMinimumRuleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SQUARE_MINIMAL_RULE_TEST.
            //
            //  Discussion:
            //
            //    SQUARE_MINIMAL_RULE_TEST tests the SQUARE_MINIMAL_RULE library.
            //    
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree;

            Console.WriteLine("");
            Console.WriteLine("SQUARE_MINIMAL_RULE_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the SQUARE_MINIMAL_RULE library.");

            degree = 8;
            square_minimal_rule_print_test(degree);

            square_minimal_rule_order_test();

            square_minimal_rule_error_max_test();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("SQUARE_MINIMAL_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void square_minimal_rule_print_test(int degree)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_PRINT_TEST tests SQUARE_MINIMAL_RULE_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double area;
            double d;
            int j;
            int order;
            double[] xyw;

            Console.WriteLine("");
            Console.WriteLine("SQUARE_MINIMAL_RULE_PRINT_TEST");
            Console.WriteLine("  SQUARE_MINIMAL_RULE_PRINT prints a quadrature rule");
            Console.WriteLine("  for the symmetric unit square.");
            Console.WriteLine("  Minimal quadrature rule for a square.");
            Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
            //
            //  Retrieve and print a symmetric quadrature rule.
            //
            xyw = MinimalRule.square_minimal_rule(degree);

            order = MinimalRule.square_minimal_rule_order(degree);

            Console.WriteLine("");
            Console.WriteLine("  Number of nodes N = " + order + "");

            Console.WriteLine("");
            Console.WriteLine("     J          X               Y               W");
            Console.WriteLine("");
            for (j = 0; j < order; j++)
            {
                Console.WriteLine("  " + j.ToString().PadLeft(4)
                                       + "  " + xyw[0 + j * 3].ToString().PadLeft(14)
                                       + "  " + xyw[1 + j * 3].ToString().PadLeft(14)
                                       + "  " + xyw[2 + j * 3].ToString().PadLeft(14) + "");
            }

            d = 0.0;
            for (j = 0; j < order; j++)
            {
                d = d + xyw[2 + j * 3];
            }

            Console.WriteLine("");
            Console.WriteLine("   Sum  " + d.ToString().PadLeft(14) + "");
            area = Integrals.squaresym_area();
            Console.WriteLine("  Area  " + area.ToString().PadLeft(14) + "");

        }

        static void square_minimal_rule_order_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_ORDER_TEST tests SQUARE_MINIMAL_RULE_ORDER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    22 February 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Mattia Festa, Alvise Sommariva,
            //    Computing almost minimal formulas on the square,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 17, Number 236, November 2012, pages 4296-4302.
            //
        {
            int degree;
            int degree_max;
            int order;

            Console.WriteLine("");
            Console.WriteLine("SQUARE_MINIMAL_RULE_ORDER_TEST");
            Console.WriteLine("  Print the order (number of points) for each");
            Console.WriteLine("  minimal square rule.");

            degree_max = MinimalRule.square_minimal_rule_degree_max();

            Console.WriteLine("");
            Console.WriteLine(" Degree  Order");
            Console.WriteLine("");
            for (degree = 0; degree <= degree_max; degree++)
            {
                order = MinimalRule.square_minimal_rule_order(degree);
                Console.WriteLine("  " + degree
                                       + "  " + order + "");
            }

        }

        static void square_minimal_rule_error_max_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_ERROR_MAX_TEST tests SQUARE_MINIMAL_RULE_ERROR_MAX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    22 February 2018
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Mattia Festa, Alvise Sommariva,
            //    Computing almost minimal formulas on the square,
            //    Journal of Computational and Applied Mathematics,
            //    Volume 17, Number 236, November 2012, pages 4296-4302.
            //
        {
            int degree;
            int degree_max;
            double error_max;
            int m_num;

            Console.WriteLine("");
            Console.WriteLine("SQUARE_MINIMAL_RULE_ERROR_MAX_TEST");
            Console.WriteLine("  SQUARE_MINIMAL_RULE_ERROR_MAX computes the maximum");
            Console.WriteLine("  error for a rule that should be exact for all monomials");
            Console.WriteLine("  up to a given value of DEGREE.");

            degree_max = MinimalRule.square_minimal_rule_degree_max();

            Console.WriteLine("");
            Console.WriteLine(" Degree  Monomials  Error Max");
            Console.WriteLine("");

            for (degree = 0; degree <= degree_max; degree++)
            {
                error_max = MinimalRule.square_minimal_rule_error_max(degree);
                m_num = ((degree + 1) * (degree + 2)) / 2;
                Console.WriteLine("   " + degree.ToString().PadLeft(4)
                                        + "       " + m_num.ToString().PadLeft(4)
                                        + "  " + error_max + "");
            }

        }
    }
}