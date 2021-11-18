using System;

namespace Burkardt.Square;

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

        switch (degree)
        {
            case 0:
                xw = smr00();
                break;
            case 1:
                xw = smr01();
                break;
            case 2:
                xw = smr02();
                break;
            case 3:
                xw = smr03();
                break;
            case 4:
                xw = smr04();
                break;
            case 5:
                xw = smr05();
                break;
            case 6:
                xw = smr06();
                break;
            case 7:
                xw = smr07();
                break;
            case 8:
                xw = smr08();
                break;
            case 9:
                xw = smr09();
                break;
            case 10:
                xw = smr10();
                break;
            case 11:
                xw = smr11();
                break;
            case 12:
                xw = smr12();
                break;
            case 13:
                xw = smr13();
                break;
            case 14:
                xw = smr14();
                break;
            case 15:
                xw = smr15();
                break;
            case 16:
                xw = smr16();
                break;
            case 17:
                xw = smr17();
                break;
            case 18:
                xw = smr18();
                break;
            case 19:
                xw = smr19();
                break;
            case 20:
                xw = smr20();
                break;
            case 21:
                xw = smr21();
                break;
            case 22:
                xw = smr22();
                break;
            case 23:
                xw = smr23();
                break;
            case 24:
                xw = smr24();
                break;
            case 25:
                xw = smr25();
                break;
            case 26:
                xw = smr26();
                break;
            case 27:
                xw = smr27();
                break;
            case 28:
                xw = smr28();
                break;
            case 29:
                xw = smr29();
                break;
            case 30:
                xw = smr30();
                break;
            case 31:
                xw = smr31();
                break;
            case 32:
                xw = smr32();
                break;
            case 33:
                xw = smr33();
                break;
            case 34:
                xw = smr34();
                break;
            case 35:
                xw = smr35();
                break;
            case 36:
                xw = smr36();
                break;
            case 37:
                xw = smr37();
                break;
            case 38:
                xw = smr38();
                break;
            case 39:
                xw = smr39();
                break;
            case 40:
                xw = smr40();
                break;
            case 41:
                xw = smr41();
                break;
            case 42:
                xw = smr42();
                break;
            case 43:
                xw = smr43();
                break;
            case 44:
                xw = smr44();
                break;
            case 45:
                xw = smr45();
                break;
            case 46:
                xw = smr46();
                break;
            case 47:
                xw = smr47();
                break;
            case 48:
                xw = smr48();
                break;
            case 49:
                xw = smr49();
                break;
            case 50:
                xw = smr50();
                break;
            case 51:
                xw = smr51();
                break;
            case 52:
                xw = smr52();
                break;
            case 53:
                xw = smr53();
                break;
            case 54:
                xw = smr54();
                break;
            case 55:
                xw = smr55();
                break;
            default:
                xw = null;
                Console.WriteLine("");
                Console.WriteLine("SQUARE_MINIMAL_RULE - Fatal error!");
                Console.WriteLine("  0 <= DEGREE <= 55 is required.");
                break;
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

        int order = square_minimal_rule_order(degree);
        double[] xyw = square_minimal_rule(degree);

        double error_max = 0.0;

        for (d = 0; d <= degree; d++)
        {
            int i;
            for (i = 0; i <= d; i++)
            {
                int j = d - i;
                e[0] = i;
                e[1] = j;
                double exact = Integrals.squaresym_monomial_integral(e);
                double s = 0.0;
                int k;
                for (k = 0; k < order; k++)
                {
                    s += xyw[2 + k * 3] * Math.Pow(xyw[0 + k * 3], i) * Math.Pow(xyw[1 + k * 3], j);
                }

                double err = Math.Abs(exact - s);
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

        int degree_max = square_minimal_rule_degree_max ( );

        int order = degree switch
        {
            < 0 => 0,
            _ => degree_max < degree ? 0 : order_list[degree]
        };

        return order;
    }

}