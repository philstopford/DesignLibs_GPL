using System;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class LevelToOrder
{
    public static void level_to_order_default(int dim_num, int[] level, int[] rule,
            ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_DEFAULT: default growth.
        //
        //  Discussion:
        //
        //    This function uses:
        //
        //    * exponential growth rates for fully nested quadrature rules,
        //      ( "CC", "F2", "GP");
        //
        //    * linear growth rates for other rules.
        //      ( "GL", "GH", "GGH", "LG", "GLG", "GJ", "GW" ).
        //
        //    * slow exponential growth alternative for fully nested rules:
        //      ("CC_SE", "F2_SE", "GP_SE").
        //
        //    * moderate exponential growth alternative for fully nested rules:
        //      ("CC_ME", "F2_ME", "GP_ME").
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the 1D levels.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
        //     2, "F2",  Fejer Type 2, Open Fully Nested rule.
        //     3, "GP",  Gauss Patterson, Open Fully Nested rule.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
        //     7, "LG",  Gauss Laguerre, Open Non Nested rule.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
        //    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
        //    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested rule.
        //    12, "F2_SE", Fejer Type 2 Slow Exponential, Open Fully Nested rule.
        //    13, "GP_SE", Gauss Patterson Slow Exponential, Open Fully Nested rule.
        //    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested rule.
        //    15, "F2_ME", Fejer Type 2 Moderate Exponential, Open Fully Nested rule.
        //    16, "GP_ME", Gauss Patterson Moderate Exponential, Open Fully Nested rule.
        //    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
        //
        //    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level[dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_DEFAULT - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                default:
                    int o;
                    int p;
                    switch (rule[dim])
                    {
                        case 1 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 1:
                            order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                            break;
                        case 2:
                        case 3:
                            order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                            break;
                        case 4:
                        case 5:
                        case 6:
                        case 7:
                        case 8:
                        case 9:
                        case 10:
                            order[dim] = 2 * level[dim] + 1;
                            break;
                        case 11:
                        {
                            switch (level[dim])
                            {
                                case 0:
                                    o = 1;
                                    break;
                                default:
                                {
                                    o = 2;
                                    while (o < 2 * level[dim] + 1)
                                    {
                                        o = 2 * (o - 1) + 1;
                                    }

                                    break;
                                }
                            }

                            order[dim] = o;
                            break;
                        }
                        case 12:
                        {
                            o = 1;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * o + 1;
                            }

                            order[dim] = o;
                            break;
                        }
                        case 13 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 13:
                        {
                            p = 5;
                            o = 3;
                            while (p < 2 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }

                            order[dim] = o;
                            break;
                        }
                        case 14:
                        {
                            switch (level[dim])
                            {
                                case 0:
                                    o = 1;
                                    break;
                                default:
                                {
                                    o = 2;
                                    while (o < 4 * level[dim] + 1)
                                    {
                                        o = 2 * (o - 1) + 1;
                                    }

                                    break;
                                }
                            }

                            order[dim] = o;
                            break;
                        }
                        case 15:
                        {
                            o = 1;
                            while (o < 4 * level[dim] + 1)
                            {
                                o = 2 * o + 1;
                            }

                            order[dim] = o;
                            break;
                        }
                        case 16 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 16:
                        {
                            p = 5;
                            o = 3;
                            while (p < 4 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }

                            order[dim] = o;
                            break;
                        }
                        case 17:
                            order[dim] = 2 * level[dim] + 1;
                            break;
                        default:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_DEFAULT - Fatal error!");
                            Console.WriteLine("  Unexpected value of RULE["
                                              + dim + "] = " + rule[dim] + ".");
                            return;
                    }

                    break;
            }
        }
    }

    public static void level_to_order_exponential(int dim_num, int[] level, int[] rule,
            ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXPONENTIAL: exponential growth.
        //
        //  Discussion:
        //
        //    The user must preallocate space for the output array ORDER.
        //
        //    Closed rules:
        //
        //      O(0) = 1
        //      O(L) = 2^L + 1;
        //
        //      O = 1, 3, 5, 9, 17, 33, ...
        //
        //    Open rules:
        //
        //      O(L) = 2^(L+1) - 1;
        //
        //      O = 1, 3, 7, 15, 31, 63, ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the 1D levels.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
        //     2, "F2",  Fejer Type 2, Open Fully Nested rule.
        //     3, "GP",  Gauss Patterson, Open Fully Nested rule.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
        //     7, "LG",  Gauss Laguerre, Open Non Nested rule.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
        //    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
        //    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested rule.
        //    12, "F2_SE", Fejer Type 2 Slow Exponential, Open Fully Nested rule.
        //    13, "GP_SE", Gauss Patterson Slow Exponential, Open Fully Nested rule.
        //    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested rule.
        //    15, "F2_ME", Fejer Type 2 Moderate Exponential, Open Fully Nested rule.
        //    16, "GP_ME", Gauss Patterson Moderate Exponential, Open Fully Nested rule.
        //    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
        //
        //    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level[dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                default:
                    switch (rule[dim])
                    {
                        case 1 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 1:
                            order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                            break;
                        case 2:
                        case 3:
                        case 4:
                        case 5:
                        case 6:
                        case 7:
                        case 8:
                        case 9:
                        case 10:
                            order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                            break;
                        case 11 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 11:
                            order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                            break;
                        case 12:
                        case 13:
                            order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                            break;
                        case 14 when level[dim] == 0:
                            order[dim] = 1;
                            break;
                        case 14:
                            order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                            break;
                        case 15:
                        case 16:
                            order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                            break;
                        case 17:
                            order[dim] = (int) Math.Pow(2, level[dim] + 1);
                            break;
                        default:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!");
                            Console.WriteLine("  Unexpected value of RULE["
                                              + dim + "] = " + rule[dim] + ".");
                            return;
                    }

                    break;
            }
        }
    }

    public static void level_to_order_exponential_slow(int dim_num, int[] level, int[] rule,
            ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXPONENTIAL_SLOW: slow exponential growth;
        //
        //  Discussion:
        //
        //    We seek a sequence of quadrature rules with two opposing constraints:
        //    * a measured rise in polynomial precision with increasing level;
        //    * a control on the increase in (new) points per level;
        //
        //    Essentially, we are trying to keep some of the advantages of nesting,
        //    while moderating the cost of the explosive growth in order that occurs
        //    due to the repeated order doubling of nesting.
        //
        //    We wish the number of points at a given level L to be "about" 2 * L + 1,
        //    but we also wish the rules to be completely nested.
        //
        //    One way to do this is to start with a nested family of rules, whose
        //    order will tend to grow exponentially (doubling from one to the next),
        //    but simply to REPEAT each rule as many times as possible.  We move to
        //    the next rule only when the desired precision 2 * L + 1 exceeds the
        //    precision of the current rule.
        //
        //    For both the Clenshaw Curtis and Fejer Type 2 rules, the order and
        //    precision are the same if the order is odd.   That is, an 11 point rule
        //    will integrate exactly all polynomials up to and including degree 11.
        //
        //    For Gauss Patterson rules, the relationship between order and precision
        //    is somewhat more complicated.  For that rule, we take the philosophy
        //    that at each level L, we wish to choose the rule of smallest order
        //    so that the precision of 2 * L + 1 is guaranteed.
        //
        //     L    2*L+1  CC Order    F2 Order    GP Order/Precision
        //
        //     0        1         1           1        1/1
        //     1        3         3           3        3/5
        //     2        5         5           7        3/5
        //     3        7         9           7        7/11
        //     4        9         9          15        7/11
        //     5       11        17          15        7/11
        //     6       13        17          15       15/23
        //     7       15        17          15       15/23
        //     8       17        17          31       15/23
        //     9       19        33          31       15/23
        //    10       21        33          31       15/23
        //    11       23        33          31       15/23
        //    12       25        33          31       31/47
        //    13       27        33          31       31/47
        //    14       29        33          31       31/47
        //    15       31        33          31       31/47
        //    16       33        33          63       31/47
        //    17       35        65          63       31/47
        //    18       37        65          63       31/47
        //    19       39        65          63       31/47
        //    20       41        65          63       31/47
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Knut Petras,
        //    Smolyak Cubature of Given Polynomial Degree with Few Nodes
        //    for Increasing Dimension,
        //    Numerische Mathematik,
        //    Volume 93, Number 4, February 2003, pages 729-753.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the 1D levels.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
        //     2, "F2",  Fejer Type 2, Open Fully Nested rule.
        //     3, "GP",  Gauss Patterson, Open Fully Nested rule.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
        //     7, "LG",  Gauss Laguerre, Open Non Nested rule.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
        //    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
        //    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested rule.
        //    12, "F2_SE", Fejer Type 2 Slow Exponential, Open Fully Nested rule.
        //    13, "GP_SE", Gauss Patterson Slow Exponential, Open Fully Nested rule.
        //    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested rule.
        //    15, "F2_ME", Fejer Type 2 Moderate Exponential, Open Fully Nested rule.
        //    16, "GP_ME", Gauss Patterson Moderate Exponential, Open Fully Nested rule.
        //    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
        //
        //    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level[dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL_SLOW - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
            }
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            int o;
            switch (rule[dim])
            {
                case 1:
                case 11:
                case 14:
                case 17:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            o = 2;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 3:
                case 13:
                case 16:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            int p = 5;
                            o = 3;
                            while (p < 2 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                default:
                {
                    o = 1;
                    while (o < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
            }

            order[dim] = o;
        }
    }

    public static void level_to_order_linear(int dim_num, int[] level, int[] rule,
            ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_LINEAR: linear growth.
        //
        //  Discussion:
        //
        //    The user must preallocate space for the output array ORDER.
        //
        //      O(L) = 2 * L + 1;
        //
        //      O = 1, 3, 5, 7, 9, ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the 1D levels.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested rule.
        //     2, "F2",  Fejer Type 2, Open Fully Nested rule.
        //     3, "GP",  Gauss Patterson, Open Fully Nested rule.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested rule.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested rule.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested rule.
        //     7, "LG",  Gauss Laguerre, Open Non Nested rule.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested rule.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested rule.
        //    10, "GW",  Golub Welsch, (presumed) Open Non Nested rule.
        //    11, "CC_SE", Clenshaw Curtis Slow Exponential, Closed Fully Nested rule.
        //    12, "F2_SE", Fejer Type 2 Slow Exponential, Open Fully Nested rule.
        //    13, "GP_SE", Gauss Patterson Slow Exponential, Open Fully Nested rule.
        //    14, "CC_ME", Clenshaw Curtis Moderate Exponential, Closed Fully Nested rule.
        //    15, "F2_ME", Fejer Type 2 Moderate Exponential, Open Fully Nested rule.
        //    16, "GP_ME", Gauss Patterson Moderate Exponential, Open Fully Nested rule.
        //    17, "CCN", Clenshaw Curtis Nested, Linear, Closed Fully Nested rule.
        //
        //    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level[dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_LINEAR - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
            }
        }

        for (dim = 0; dim < dim_num; dim++)
        {
            order[dim] = 2 * level[dim] + 1;
        }

    }

    public static int level_to_order_exp_cc(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXP_CC is used for Clenshaw-Curtis type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be closed (including both endpoints
        //    except for the level 0 rule) and having a precision
        //    behavior typical of Clenshaw Curtis rules, namely, the ORDER-point
        //    rule is exact for polynomials of degree less than ORDER, and if
        //    ORDER is odd, then the exactness includes polynomials of degree ORDER
        //    as well.
        //
        //    LEVEL  ORDER  ORDER  ORDER
        //           G = 0  G = 1  G = 2
        //    -----  -----  -----  -----
        //        0      1      1      1
        //        1      3      5      3
        //        2      5      9      5
        //        3      9     17      9
        //        4      9     17     17
        //        5     17     33     33
        //        6     17     33     65
        //        7     17     33    129
        //        8     17     33    257       
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //    2, full growth.
        //
        //    Output, int LEVEL_TO_ORDER_EXP_CC, the order of the rule.
        //
    {
        int o;
        switch (growth)
        {
            //
            //  Slow exponential growth.
            //
            case 0 when level == 0:
                o = 1;
                break;
            case 0:
            {
                o = 2;
                while (o < 2 * level + 1)
                {
                    o = 2 * (o - 1) + 1;
                }

                break;
            }
            //
            //  Moderate Exponential Growth.
            //
            case 1 when level == 0:
                o = 1;
                break;
            case 1:
            {
                o = 2;
                while (o < 4 * level + 1)
                {
                    o = 2 * (o - 1) + 1;
                }

                break;
            }
            //
            //  Full Exponential Growth.
            //
            case 2 when level == 0:
                o = 1;
                break;
            case 2:
                o = (int) Math.Pow(2, level) + 1;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_CC - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_exp_f2(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXP_F2 is used for Fejer 2 type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint)
        //    and having a precision behavior typical of Fejer Type 2
        //    rules, namely, the ORDER-point rule is exact for polynomials of degree 
        //    less than ORDER, and if ORDER is odd, then the exactness includes 
        //    polynomials of degree ORDER as well.
        // 
        //    LEVEL  ORDER  ORDER  ORDER
        //           G = 0  G = 1  G = 2
        //
        //        0      1      1      1
        //        1      3      7      3
        //        2      7     15      7
        //        3      7     15     15
        //        4     15     31     31
        //        5     15     31     63
        //        6     15     31    127
        //        7     15     31    255
        //        8     31     63    511  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //    2, full growth.
        //
        //    Output, int LEVEL_TO_ORDER_EXP_F2, the order of the rule.
        //
    {
        int o;
        switch (growth)
        {
            //
            //  Slow exponential growth.
            //
            case 0 when level == 0:
                o = 1;
                break;
            case 0:
            {
                o = 1;
                while (o < 2 * level + 1)
                {
                    o = 2 * o + 1;
                }

                break;
            }
            //
            //  Moderate Exponential Growth.
            //
            case 1 when level == 0:
                o = 1;
                break;
            case 1:
            {
                o = 1;
                while (o < 4 * level + 1)
                {
                    o = 2 * o + 1;
                }

                break;
            }
            //
            //  Full Exponential Growth.
            //
            case 2 when level == 0:
                o = 1;
                break;
            case 2:
                o = (int) Math.Pow(2, level + 1) - 1;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_F2 - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_exp_gauss(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXP_GAUSS is used for Gauss type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint),
        //    and having a precision behavior typical of Gauss rules, namely, the 
        //    ORDER-point rule is exact for polynomials of degree less than 2 * ORDER.
        //
        //    LEVEL  ORDER  ORDER  ORDER
        //           G = 0  G = 1  G = 2
        //
        //        0      1      1      1
        //        1      3      3      3
        //        2      3      7      7
        //        3      7      7     15
        //        4      7     15     31
        //        5      7     15     63
        //        6      7     15    127
        //        7     15     15    255
        //        8     15     31    511
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //    2, full growth.
        //
        //    Output, int LEVEL_TO_ORDER_EXP_GAUSS, the order of the rule.
        //
    {
        int o;
        switch (growth)
        {
            //
            //  Slow exponential growth.
            //
            case 0 when level == 0:
                o = 1;
                break;
            case 0:
            {
                o = 1;
                while (2 * o - 1 < 2 * level + 1)
                {
                    o = 2 * o + 1;
                }

                break;
            }
            //
            //  Moderate Exponential Growth.
            //
            case 1 when level == 0:
                o = 1;
                break;
            case 1:
            {
                o = 1;
                while (2 * o - 1 < 4 * level + 1)
                {
                    o = 2 * o + 1;
                }

                break;
            }
            //
            //  Full Exponential Growth.
            //
            case 2 when level == 0:
                o = 1;
                break;
            case 2:
                o = (int) Math.Pow(2, level + 1) - 1;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_GAUSS - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_exp_gp(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXP_GP is used for Gauss-Patterson type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint)
        //    and having a precision behavior typical of Gauss Patterson rules.
        //
        //    Note that there are onlly 9 rules in the family, and so it is possible to
        //    specify input for which the function will fail.
        //
        //    LEVEL  ORDER  ORDER  ORDER
        //           G = 0  G = 1  G = 2
        //
        //        0      1      1      1
        //        1      3      3      3
        //        2      3      7      7
        //        3      7     15     15
        //        4      7     15     31
        //        5      7     15     63
        //        6     15     31    127
        //        7     15     31    255
        //        8     15     31    511
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //    2, full growth.
        //
        //    Output, int LEVEL_TO_ORDER_EXP_GP, the order of the rule.
        //
    {
        int o;
        int p;
        switch (growth)
        {
            //
            //  Slow exponential growth.
            //
            case 0 when level == 0:
                o = 1;
                break;
            case 0:
            {
                p = 5;
                o = 3;
                while (p < 2 * level + 1)
                {
                    p = 2 * p + 1;
                    o = 2 * o + 1;
                    switch (o)
                    {
                        case > 511:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                            Console.WriteLine("  Request for unavailable Patterson rule.");
                            return 1;
                    }
                }

                break;
            }
            //
            //  Moderate Exponential Growth.
            //
            case 1 when level == 0:
                o = 1;
                break;
            case 1:
            {
                p = 5;
                o = 3;
                while (p < 4 * level + 1)
                {
                    p = 2 * p + 1;
                    o = 2 * o + 1;
                    switch (o)
                    {
                        case > 511:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                            Console.WriteLine("  Request for unavailable Patterson rule.");
                            return 1;
                    }
                }

                break;
            }
            //
            //  Full Exponential Growth.
            //
            case 2 when level == 0:
                o = 1;
                break;
            case 2:
            {
                o = (int) Math.Pow(2, level + 1) - 1;
                switch (o)
                {
                    case > 511:
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                        Console.WriteLine("  Request for unavailable Patterson rule.");
                        return 1;
                }

                break;
            }
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_exp_hgk(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_EXP_HGK is used for Hermite Genz-Keister type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint)
        //    and having a precision behavior typical of Hermite Genz-Keister rules.
        //
        //    Note that there are only 6 rules in the family, and so it is possible to
        //    specify input for which the function will fail.
        //
        //    LEVEL  ORDER  ORDER  ORDER
        //           G = 0  G = 1  G = 2
        //
        //        0      1      1      1
        //        1      3      3      3
        //        2      3      9      9
        //        3      9      9     19
        //        4      9     19     35
        //        5      9     19     43
        //        6      9     19     --
        //        7      9     19     --
        //        8     19     35     --
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //    2, full growth.
        //
        //    Output, int LEVEL_TO_ORDER_EXP_HGK, the order of the rule.
        //
    {
        int l;
        int o;
        int[] o_hgk = {1, 3, 9, 19, 35, 43};
        int p;
        int[] p_hgk = {1, 5, 15, 29, 51, 67};
        switch (growth)
        {
            //
            //  Slow exponential growth.
            //
            case 0:
            {
                l = 0;
                p = p_hgk[l];
                o = o_hgk[l];
                while (p < 2 * level + 1)
                {
                    l += 1;
                    switch (l)
                    {
                        case > 5:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                            Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                            return 1;
                        default:
                            p = p_hgk[l];
                            o = o_hgk[l];
                            break;
                    }
                }

                break;
            }
            case 1:
            {
                l = 0;
                p = p_hgk[l];
                o = o_hgk[l];
                while (p < 4 * level + 1)
                {
                    l += 1;
                    switch (l)
                    {
                        case > 5:
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                            Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                            return 1;
                        default:
                            p = p_hgk[l];
                            o = o_hgk[l];
                            break;
                    }
                }

                break;
            }
            case 2:
            {
                l = level;
                l = Math.Max(l, 0);
                switch (l)
                {
                    case > 5:
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                        Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                        return 1;
                }

                o = o_hgk[l];
                break;
            }
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_linear_nn(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_LINEAR_NN is used for non-nested Gauss type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint),
        //    non-nested, and having a precision behavior typical of Gauss rules.
        //
        //    LEVEL  ORDER  ORDER
        //           G = 0  G = 1
        //
        //        0      1      1
        //        1      2      3
        //        2      3      5
        //        3      4      7
        //        4      5      9
        //        5      6     11
        //        6      7     13
        //        7      8     15
        //        8      9     17
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //
        //    Output, int LEVEL_TO_ORDER_LINEAR_NN, the order of the rule.
        //
    {
        int o;
        switch (growth)
        {
            //
            //  Slow linear growth.
            //
            case 0:
                o = level + 1;
                break;
            //
            //  Moderate linear growth.
            //
            case 1:
            case 2:
                o = 2 * level + 1;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_LINEAR_NN - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static int level_to_order_linear_wn(int level, int growth)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_LINEAR_WN is used for weakly-nested Gauss type rules.
        //
        //  Discussion:
        //
        //    Rules of this type are assumed to be open (not including either endpoint),
        //    nested, and having a precision behavior typical of Gauss rules.
        //
        //    We assume the rules are to be generated with an odd number of points,
        //    and that all the rules will share a single point, namely 0.
        //
        //    Note that the "moderate growth" option for this function results in the
        //    same values as the moderate growth option for LEVEL_TO_ORDER_LINEAR_NN.
        //
        //    LEVEL  ORDER  ORDER
        //           G = 0  G = 1
        //
        //        0      1      1
        //        1      3      3
        //        2      3      5
        //        3      5      7
        //        4      5      9
        //        5      7     11
        //        6      7     13
        //        7      9     15
        //        8      9     17
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level of the rule.
        //
        //    Input, int GROWTH, the growth policy:
        //    0, slow growth;
        //    1, moderate growth;
        //
        //    Output, int LEVEL_TO_ORDER_LINEAR_WN, the order of the rule.
        //
    {
        int o;
        switch (growth)
        {
            //
            //  Slow growth.
            //
            case 0:
                o = 2 * ((level + 1) / 2) + 1;
                break;
            case 1:
            case 2:
                o = 2 * level + 1;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_LINEAR_WN - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return 1;
        }

        return o;
    }

    public static void level_growth_to_order(int dim_num, int[] level, int[] rule,
            int[] growth, ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_GROWTH_TO_ORDER: convert Level and Growth to Order.
        //
        //  Discussion:
        //
        //    This function is given level, rule, and growth information
        //    for each dimension of a quadrature rule, and determines the
        //    corresponding order of the rule in each dimension.
        //
        //    This is a revised version of LEVEL_GROWTH_TO_ORDER.
        //
        //    In particular, it revises the interpretation of the RULE vector as 
        //    far as the values 10, 11, and 12 are concerned.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the 1D levels.
        //
        //    Input, int RULE[DIM_NUM], the rule in each dimension.
        //     1, "CC",  Clenshaw Curtis, Closed Fully Nested.
        //     2, "F2",  Fejer Type 2, Open Fully Nested.
        //     3, "GP",  Gauss Patterson, Open Fully Nested.
        //     4, "GL",  Gauss Legendre, Open Weakly Nested.
        //     5, "GH",  Gauss Hermite, Open Weakly Nested.
        //     6, "GGH", Generalized Gauss Hermite, Open Weakly Nested.
        //     7, "LG",  Gauss Laguerre, Open Non Nested.
        //     8, "GLG", Generalized Gauss Laguerre, Open Non Nested.
        //     9, "GJ",  Gauss Jacobi, Open Non Nested.
        //    10, "HGK", Hermite Genz-Keister, Open Fully Nested.
        //    11, "UO",  User supplied Open, presumably Non Nested.
        //    12, "UC",  User supplied Closed, presumably Non Nested.
        //
        //    Input, int GROWTH[DIM_NUM], the desired growth in each dimension.
        //    0, "DF", default growth associated with this quadrature rule;
        //    1, "SL", slow linear, L+1;
        //    2  "SO", slow linear odd, O=1+2((L+1)/2)
        //    3, "ML", moderate linear, 2L+1;
        //    4, "SE", slow exponential;
        //    5, "ME", moderate exponential;
        //    6, "FE", full exponential.
        //
        //    Output, int ORDER[DIM_NUM], the 1D orders (number of points).
        //
    {
        int dim;
        int o = 0;
        int[] o_hgk = {1, 3, 9, 19, 35, 43};
        int[] p_hgk = {1, 5, 15, 29, 51, 67};
        //
        //  Check the input.
        //
        for (dim = 0; dim < dim_num; dim++)
        {
            switch (level[dim])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
            }

            switch (rule[dim])
            {
                case < 1:
                case > 12:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Illegal value of RULE[DIM]!");
                    Console.WriteLine("  RULE[" + dim + "] = " + rule[dim] + "");
                    return;
            }

            switch (growth[dim])
            {
                case < 0:
                case > 6:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Illegal value of GROWTH[DIM]!");
                    Console.WriteLine("  GROWTH[" + dim + "] = " + growth[dim] + "");
                    return;
            }
        }

        //
        //  Compute the order vector.
        //
        for (dim = 0; dim < dim_num; dim++)
        {
            int l;
            int p;
            switch (rule[dim])
            {
                //
                //  CC
                //  Default is Moderate Exponential Growth.
                //
                case 1 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 1 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 1 when growth[dim] == 3:
                    o = 2 * level[dim] + 1;
                    break;
                case 1 when growth[dim] == 4:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            o = 2;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 1 when growth[dim] == 5 || growth[dim] == 0:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            o = 2;
                            while (o < 4 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 1:
                {
                    o = growth[dim] switch
                    {
                        6 when level[dim] == 0 => 1,
                        6 => (int) Math.Pow(2, level[dim]) + 1,
                        _ => o
                    };

                    break;
                }
                //
                //  F2
                //  Default is Moderate Exponential Growth.
                //
                case 2 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 2 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 2 when growth[dim] == 3:
                    o = 2 * level[dim] + 1;
                    break;
                case 2 when growth[dim] == 4:
                {
                    o = 1;
                    while (o < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 2 when growth[dim] == 5 || growth[dim] == 0:
                {
                    o = 1;
                    while (o < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 2:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GP
                //  Default is Moderate Exponential Growth.
                //
                case 3 when growth[dim] == 1:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 1 for rule 3 not available!");
                    return;
                case 3 when growth[dim] == 2:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 2 for rule 3 not available!");
                    return;
                case 3 when growth[dim] == 3:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 3 for rule 3 not available!");
                    return;
                case 3 when growth[dim] == 4:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            p = 5;
                            o = 3;
                            while (p < 2 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 3 when growth[dim] == 5 || growth[dim] == 0:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            p = 5;
                            o = 3;
                            while (p < 4 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 3:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GL
                //  Default is Moderate Linear Growth.
                //
                case 4 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 4 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 4 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 4 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 4 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 4:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GH
                //  Default is Moderate Linear Growth.
                //
                case 5 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 5 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 5 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 5 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 5 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 5:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GGH
                //  Default is Moderate Linear Growth.
                //
                case 6 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 6 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 6 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 6 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 6 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 6:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  LG
                //  Default is Moderate Linear Growth.
                //
                case 7 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 7 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 7 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 7 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 7 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 7:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GLG
                //  Default is Moderate Linear Growth.
                //
                case 8 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 8 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 8 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 8 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 8 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 8:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  GJ
                //  Default is Moderate Linear Growth.
                //
                case 9 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 9 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 9 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 9 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 9 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 9:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  HGK
                //  Default is Moderate Exponential Growth.
                //  Exponential growth is interpreted to mean simply take successive rules.
                //
                case 10 when growth[dim] == 1:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 1 for rule 10 not available!");
                    return;
                case 10 when growth[dim] == 2:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 2 for rule 10 not available!");
                    return;
                case 10 when growth[dim] == 3:
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Growth rate 3 for rule 10 not available!");
                    return;
                case 10 when growth[dim] == 4:
                {
                    l = 0;
                    p = p_hgk[l];
                    o = o_hgk[l];
                    while (p < 2 * level[dim] + 1)
                    {
                        l += 1;
                        switch (l)
                        {
                            case > 5:
                                Console.WriteLine("");
                                Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                                Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                                return;
                            default:
                                p = p_hgk[l];
                                o = o_hgk[l];
                                break;
                        }
                    }

                    break;
                }
                case 10 when growth[dim] == 5 || growth[dim] == 0:
                {
                    l = 0;
                    p = p_hgk[l];
                    o = o_hgk[l];
                    while (p < 4 * level[dim] + 1)
                    {
                        l += 1;
                        switch (l)
                        {
                            case > 5:
                                Console.WriteLine("");
                                Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                                Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                                return;
                            default:
                                p = p_hgk[l];
                                o = o_hgk[l];
                                break;
                        }
                    }

                    break;
                }
                case 10:
                {
                    switch (growth[dim])
                    {
                        case 6:
                        {
                            l = level[dim];
                            l = Math.Max(l, 0);
                            switch (l)
                            {
                                case > 5:
                                    Console.WriteLine("");
                                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                                    Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                                    return;
                            }

                            o = o_hgk[l];
                            break;
                        }
                    }

                    break;
                }
                //
                //  UO
                //  Default is Moderate Linear Growth.
                //  We assume the rule is of OPEN type and that it
                //  has a precision typical of Gauss rules.
                //
                case 11 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 11 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 11 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 11 when growth[dim] == 4:
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 11 when growth[dim] == 5:
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    break;
                }
                case 11:
                {
                    o = growth[dim] switch
                    {
                        6 => (int) Math.Pow(2, level[dim] + 1) - 1,
                        _ => o
                    };

                    break;
                }
                //
                //  UC
                //  Default is Moderate Linear Growth.
                //  We assume the rule is of CLOSED type and that it
                //  has a precision typical of Clenshaw-Curtis rules.
                //
                case 12 when growth[dim] == 1:
                    o = level[dim] + 1;
                    break;
                case 12 when growth[dim] == 2:
                    o = 2 * ((level[dim] + 1) / 2) + 1;
                    break;
                case 12 when growth[dim] == 3 || growth[dim] == 0:
                    o = 2 * level[dim] + 1;
                    break;
                case 12 when growth[dim] == 4:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            o = 2;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 12 when growth[dim] == 5:
                {
                    switch (level[dim])
                    {
                        case 0:
                            o = 1;
                            break;
                        default:
                        {
                            o = 2;
                            while (o < 4 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }

                            break;
                        }
                    }

                    break;
                }
                case 12:
                {
                    o = growth[dim] switch
                    {
                        6 when level[dim] == 0 => 1,
                        6 => (int) Math.Pow(2, level[dim]) + 1,
                        _ => o
                    };

                    break;
                }
            }

            order[dim] = o;
        }
    }

    public static int[] index_level_own(int level, int level_max, int dim_num, int point_num,
            int[] grid_index, int[] grid_base)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_LEVEL_OWN: determine first level at which given index is generated.
        //
        //  Discussion:
        //
        //    We are constructing a sparse grid of OWN points.  The grid
        //    is built up of product grids, with a characteristic LEVEL.  
        //
        //    We are concerned with identifying points in this product grid which
        //    have actually been generated previously, on a lower value of LEVEL.
        //
        //    This routine determines the lowest value of LEVEL at which each of
        //    the input points would be generated.
        //
        //    In 1D, given LEVEL, the number of points is ORDER = 2**(LEVEL+1) + 1,
        //    (except that LEVEL = 0 implies ORDER = 1), the BASE is (ORDER-1)/2, 
        //    and the point INDEX values range from -BASE to +BASE.
        //
        //    The values of INDEX and BASE allow us to determine the abstract
        //    properties of the point.  In particular, if INDEX is 0, the corresponding
        //    abscissa is 0, the special "nested" value we need to take care of.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int LEVEL, the level at which these points were 
        //    generated.  LEVEL_MIN <= LEVEL <= LEVEL_MAX.
        //
        //    Input, int LEVEL_MAX, the maximum level.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points to be tested.
        //
        //    Input, int GRID_INDEX[DIM_NUM*POINT_NUM], the indices of the 
        //    points to be tested.
        //
        //    Input, int GRID_BASE[DIM_NUM], the "base", which is essentially
        //    the denominator of the index.
        //
        //    Output, int INDEX_LEVEL_OWN[POINT_NUM], the value of LEVEL at 
        //    which the point would first be generated.  This will be the same as
        //    the input value of LEVEL, unless the point has an INDEX of 0 and
        //    a corresponding BASE that is NOT zero.
        //
    {
        int point;

        int[] grid_level = new int[point_num];

        int level_min = dim_num switch
        {
            1 => level_max,
            _ => 0
        };

        //
        //  If a point has a DIM-th component whose INDEX is 0, then the 
        //  value of LEVEL at which this point would first be generated is
        //  less than LEVEL, unless the DIM-th component of GRID_BASE is 0.
        //
        for (point = 0; point < point_num; point++)
        {
            grid_level[point] = Math.Max(level, level_min);

            int dim;
            for (dim = 0; dim < dim_num; dim++)
            {
                grid_level[point] = grid_index[dim + point * dim_num] switch
                {
                    0 => Math.Max(grid_level[point] - grid_base[dim], level_min),
                    _ => grid_level[point]
                };
            }
        }

        return grid_level;
    }

    public static int index_to_level_closed(int dim_num, int[] t, int order, int level_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_TO_LEVEL_CLOSED determines the level of a point given its index.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int T[DIM_NUM], the grid indices of a point in a 1D closed rule.
        //    0 <= T[I] <= ORDER.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, int LEVEL_MAX, the level with respect to which the
        //    index applies.
        //
        //    Output, int INDEX_TO_LEVEL_CLOSED, the first level on which
        //    the point associated with the given index will appear.
        //
    {
        int dim;

        int value = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            int s = t[dim];

            s = typeMethods.i4_modp(s, order);

            int level;
            switch (s)
            {
                case 0:
                    level = 0;
                    break;
                default:
                {
                    level = level_max;

                    while (s % 2 == 0)
                    {
                        s /= 2;
                        level -= 1;
                    }

                    break;
                }
            }

            level = level switch
            {
                0 => 1,
                1 => 0,
                _ => level
            };

            value += level;
        }

        return value;
    }

    public static int index_to_level_open(int dim_num, int[] t, int order, int level_max, int tIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INDEX_TO_LEVEL_OPEN determines the level of a point given its index.
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int T[DIM_NUM], the grid index of a point.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, int LEVEL_MAX, the level with respect to which the
        //    index applies.
        //
        //    Output, int INDEX_TO_LEVEL_OPEN, the first level on which
        //    the point associated with the given index will appear.
        //
    {
        int dim;

        int value = 0;

        for (dim = 0; dim < dim_num; dim++)
        {
            int s = t[(dim + tIndex) % t.Length];

            s = typeMethods.i4_modp(s, order);

            int level;
            switch (s)
            {
                case 0:
                    level = 0;
                    break;
                default:
                {
                    level = level_max;

                    while (s % 2 == 0)
                    {
                        s /= 2;
                        level -= 1;
                    }

                    break;
                }
            }

            level = level switch
            {
                0 => 1,
                1 => 0,
                _ => level
            };

            value += level;
        }

        return value;
    }

    public static void level_to_order_closed(int dim_num, int[] level, ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_CLOSED converts a level to an order for closed rules.
        //
        //  Discussion:
        //
        //    Sparse grids can naturally be nested.  A natural scheme is to use
        //    a series of one-dimensional rules arranged in a series of "levels"
        //    whose order roughly doubles with each step.
        //
        //    The arrangement described here works naturally for the Clenshaw Curtis
        //    and Newton Cotes closed rules.  
        //
        //    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
        //    point at the center, and for all values afterwards, we use the 
        //    relationship
        //
        //      ORDER = 2^LEVEL + 1
        //
        //    The following table shows how the growth will occur:
        //
        //    Level    Order
        //
        //    0          1
        //    1          3 =  2 + 1
        //    2          5 =  4 + 1
        //    3          9 =  8 + 1
        //    4         17 = 16 + 1
        //    5         33 = 32 + 1
        //
        //    For the Clenshaw Curtis and Newton Cotes Closed rules, the point growth
        //    is nested.  If we have ORDER points on a particular LEVEL, the next
        //    level includes all these old points, plus ORDER-1 new points, formed
        //    in the gaps between successive pairs of old points.
        //
        //    Level    Order = New + Old
        //
        //    0          1   =  1  +  0
        //    1          3   =  2  +  1
        //    2          5   =  2  +  3
        //    3          9   =  4  +  5
        //    4         17   =  8  +  9
        //    5         33   = 16  + 17
        //
        //    In this routine, we assume that a vector of levels is given,
        //    and the corresponding orders are desired.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the nesting level.
        //
        //    Output, int ORDER[DIM_NUM], the order (number of points) 
        //    of the rule.
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            order[dim] = level[dim] switch
            {
                < 0 => -1,
                0 => 1,
                _ => (int) Math.Pow(2, level[dim]) + 1
            };
        }
    }

    public static void level_to_order_open(int dim_num, int[] level, ref int[] order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_TO_ORDER_OPEN converts a level to an order for open rules.
        //
        //  Discussion:
        //
        //    Sparse grids can naturally be nested.  A natural scheme is to use
        //    a series of one-dimensional rules arranged in a series of "levels"
        //    whose order roughly doubles with each step.
        //
        //    The arrangement described here works naturally for the Fejer Type 1,
        //    Fejer Type 2, Newton Cotes Open, Newton Cotes Half Open,
        //    and Gauss-Patterson rules.  It also can be used, partially, to describe
        //    the growth of Gauss-Legendre rules.
        //
        //    The idea is that we start with LEVEL = 0, ORDER = 1 indicating the single 
        //    point at the center, and for all values afterwards, we use the relationship
        //
        //      ORDER = 2**(LEVEL+1) - 1.
        //
        //    The following table shows how the growth will occur:
        //
        //    Level    Order
        //
        //    0          1
        //    1          3 =  4 - 1
        //    2          7 =  8 - 1
        //    3         15 = 16 - 1
        //    4         31 = 32 - 1
        //    5         63 = 64 - 1
        //
        //    For the Fejer Type 1, Fejer Type 2, Newton Cotes Open, 
        //    Newton Cotes Open Half, and Gauss-Patterson rules, the point growth is
        //    nested.  If we have ORDER points on a particular LEVEL, the next level 
        //    includes all these old points, plus ORDER+1 new points, formed in the 
        //    gaps between successive pairs of old points plus an extra point at each 
        //    end.
        //
        //    Level    Order = New + Old
        //
        //    0          1   =  1  +  0
        //    1          3   =  2  +  1
        //    2          7   =  4  +  3
        //    3         15   =  8  +  7
        //    4         31   = 16  + 15
        //    5         63   = 32  + 31
        //
        //    If we use a series of Gauss-Legendre rules, then there is almost no 
        //    nesting, except that the central point is shared.  If we insist on 
        //    producing a comparable series of such points, then the "nesting" behavior
        //    is as follows:
        //
        //    Level    Order = New + Old
        //
        //    0          1   =  1  +  0
        //    1          3   =  2  +  1
        //    2          7   =  6  +  1
        //    3         15   = 14  +  1
        //    4         31   = 30  +  1
        //    5         63   = 62  +  1
        //
        //    Moreover, if we consider ALL the points used in such a set of "nested" 
        //    Gauss-Legendre rules, then we must sum the "NEW" column, and we see that
        //    we get roughly twice as many points as for the truly nested rules.
        //
        //    In this routine, we assume that a vector of levels is given,
        //    and the corresponding orders are desired.
        //
        //  Modified:
        //
        //    19 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fabio Nobile, Raul Tempone, Clayton Webster,
        //    A Sparse Grid Stochastic Collocation Method for Partial Differential
        //    Equations with Random Input Data,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 46, Number 5, 2008, pages 2309-2345.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int LEVEL[DIM_NUM], the nesting level.
        //
        //    Output, int ORDER[DIM_NUM], the order (number of points) 
        //    of the rule.
        //
    {
        int dim;

        for (dim = 0; dim < dim_num; dim++)
        {
            order[dim] = level[dim] switch
            {
                < 0 => -1,
                0 => 1,
                _ => (int) Math.Pow(2, level[dim] + 1) - 1
            };
        }
    }


}