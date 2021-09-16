using System;

namespace Burkardt.Square
{
    public static partial class MinimalRule
    {
        public static double[] square_minimal_rule(int degree)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE returns a minimal rule for the square.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    20 February 2018
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
            //  Parameters:
            //
            //    Input, int DEGREE, the degree, between 0 and 55.
            //
            //    Output, double *SQUARE_MINIMAL_RULE[3*O], the rule.
            //
        {
            double[] xw;

            if (degree == 0)
            {
                xw = smr00();
            }
            else if (degree == 1)
            {
                xw = smr01();
            }
            else if (degree == 2)
            {
                xw = smr02();
            }
            else if (degree == 3)
            {
                xw = smr03();
            }
            else if (degree == 4)
            {
                xw = smr04();
            }
            else if (degree == 5)
            {
                xw = smr05();
            }
            else if (degree == 6)
            {
                xw = smr06();
            }
            else if (degree == 7)
            {
                xw = smr07();
            }
            else if (degree == 8)
            {
                xw = smr08();
            }
            else if (degree == 9)
            {
                xw = smr09();
            }
            else if (degree == 10)
            {
                xw = smr10();
            }
            else if (degree == 11)
            {
                xw = smr11();
            }
            else if (degree == 12)
            {
                xw = smr12();
            }
            else if (degree == 13)
            {
                xw = smr13();
            }
            else if (degree == 14)
            {
                xw = smr14();
            }
            else if (degree == 15)
            {
                xw = smr15();
            }
            else if (degree == 16)
            {
                xw = smr16();
            }
            else if (degree == 17)
            {
                xw = smr17();
            }
            else if (degree == 18)
            {
                xw = smr18();
            }
            else if (degree == 19)
            {
                xw = smr19();
            }
            else if (degree == 20)
            {
                xw = smr20();
            }
            else if (degree == 21)
            {
                xw = smr21();
            }
            else if (degree == 22)
            {
                xw = smr22();
            }
            else if (degree == 23)
            {
                xw = smr23();
            }
            else if (degree == 24)
            {
                xw = smr24();
            }
            else if (degree == 25)
            {
                xw = smr25();
            }
            else if (degree == 26)
            {
                xw = smr26();
            }
            else if (degree == 27)
            {
                xw = smr27();
            }
            else if (degree == 28)
            {
                xw = smr28();
            }
            else if (degree == 29)
            {
                xw = smr29();
            }
            else if (degree == 30)
            {
                xw = smr30();
            }
            else if (degree == 31)
            {
                xw = smr31();
            }
            else if (degree == 32)
            {
                xw = smr32();
            }
            else if (degree == 33)
            {
                xw = smr33();
            }
            else if (degree == 34)
            {
                xw = smr34();
            }
            else if (degree == 35)
            {
                xw = smr35();
            }
            else if (degree == 36)
            {
                xw = smr36();
            }
            else if (degree == 37)
            {
                xw = smr37();
            }
            else if (degree == 38)
            {
                xw = smr38();
            }
            else if (degree == 39)
            {
                xw = smr39();
            }
            else if (degree == 40)
            {
                xw = smr40();
            }
            else if (degree == 41)
            {
                xw = smr41();
            }
            else if (degree == 42)
            {
                xw = smr42();
            }
            else if (degree == 43)
            {
                xw = smr43();
            }
            else if (degree == 44)
            {
                xw = smr44();
            }
            else if (degree == 45)
            {
                xw = smr45();
            }
            else if (degree == 46)
            {
                xw = smr46();
            }
            else if (degree == 47)
            {
                xw = smr47();
            }
            else if (degree == 48)
            {
                xw = smr48();
            }
            else if (degree == 49)
            {
                xw = smr49();
            }
            else if (degree == 50)
            {
                xw = smr50();
            }
            else if (degree == 51)
            {
                xw = smr51();
            }
            else if (degree == 52)
            {
                xw = smr52();
            }
            else if (degree == 53)
            {
                xw = smr53();
            }
            else if (degree == 54)
            {
                xw = smr54();
            }
            else if (degree == 55)
            {
                xw = smr55();
            }
            else
            {
                xw = null;
                Console.WriteLine("");
                Console.WriteLine("SQUARE_MINIMAL_RULE - Fatal error!");
                Console.WriteLine("  0 <= DEGREE <= 55 is required.");
            }

            return xw;
        }

        public static int square_minimal_rule_degree_max()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_DEGREE_MAX returns the maximum rule degree.
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
            //  Parameters:
            //
            //    Output, integer DEGREE_MAX, the maximum degree of the minimal
            //    rules for the square.
            //
        {
            const int degree_max = 55;

            return degree_max;
        }

        public static double square_minimal_rule_error_max(int degree)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_ERROR_MAX returns the maximum error.
            //
            //  Discussion:
            //
            //    The rule of given DEGREE should theoretically have zero error
            //    for all monomials of degrees 0 <= D <= DEGREE.  This function
            //    checks every such monomial and reports the maximum error.
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
            //  Parameters:
            //
            //    Input, int DEGREE, the desired total polynomial degree exactness
            //    of the quadrature rule.
            //
            //    Output, double SQUARE_MINIMAL_RULE_ERROR_MAX, the maximum error observed 
            //    when using the rule to compute the integrals of all monomials of degree 
            //    between 0 and DEGREE.
            //
        {
            int d;
            int[] e = new int[2];
            double err;
            double error_max;
            double exact;
            int i;
            int j;
            int k;
            int order;
            double s;
            double[] xyw;

            order = square_minimal_rule_order(degree);
            xyw = square_minimal_rule(degree);

            error_max = 0.0;

            for (d = 0; d <= degree; d++)
            {
                for (i = 0; i <= d; i++)
                {
                    j = d - i;
                    e[0] = i;
                    e[1] = j;
                    exact = Integrals.squaresym_monomial_integral(e);
                    s = 0.0;
                    for (k = 0; k < order; k++)
                    {
                        s = s + xyw[2 + k * 3] * Math.Pow(xyw[0 + k * 3], i) * Math.Pow(xyw[1 + k * 3], j);
                    }

                    err = Math.Abs(exact - s);
                    if (error_max < err)
                    {
                        error_max = err;
                    }
                }
            }

            return error_max;
        }

        public static int square_minimal_rule_order ( int degree )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQUARE_MINIMAL_RULE_ORDER returns the order of a minimal square rule.
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
            //  Parameters:
            //
            //    Input, int DEGREE, the degree of the rule,
            //    between 0 and 55.
            //
            //    Output, int SQUARE_MINIMAL_RULE_ORDER, the order of the rule.
            //
        {
            int degree_max;
            int order;
            int[] order_list = {
                1,    1,    3,    4,    6,
                7,   10,   12,   16,   17,   
                22,   24,   31,   33,   40,   
                43,   52,   54,   64,   67,   
                78,   81,   93,   96,  109,  
                113,  127,  132,  146,  153,  
                167,  172,  189,  197,  211,  
                220,  238,  245,  265,  274,  
                296,  303,  326,  331,  353,  
                359,  387,  396,  417,  427,  
                454,  462,  493,  498,  530,  
                536 };

            degree_max = square_minimal_rule_degree_max ( );

            if ( degree < 0 )
            {
                order = 0;
            }
            else if ( degree_max < degree )
            {
                order = 0;
            }
            else
            {
                order = order_list[degree];
            }

            return order;
        }

    }
}