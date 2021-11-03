using System;

namespace Burkardt.Quadrature
{
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
            int o;
            int p;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (level[dim] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_DEFAULT - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                }
                else if (rule[dim] == 1)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                    }
                }
                else if (rule[dim] == 2)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 3)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 4)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 5)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 6)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 7)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 8)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 9)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 10)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else if (rule[dim] == 11)
                {
                    if (level[dim] == 0)
                    {
                        o = 1;
                    }
                    else
                    {
                        o = 2;
                        while (o < 2 * level[dim] + 1)
                        {
                            o = 2 * (o - 1) + 1;
                        }
                    }

                    order[dim] = o;
                }
                else if (rule[dim] == 12)
                {
                    o = 1;
                    while (o < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    order[dim] = o;
                }
                else if (rule[dim] == 13)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        p = 5;
                        o = 3;
                        while (p < 2 * level[dim] + 1)
                        {
                            p = 2 * p + 1;
                            o = 2 * o + 1;
                        }

                        order[dim] = o;
                    }
                }
                else if (rule[dim] == 14)
                {
                    if (level[dim] == 0)
                    {
                        o = 1;
                    }
                    else
                    {
                        o = 2;
                        while (o < 4 * level[dim] + 1)
                        {
                            o = 2 * (o - 1) + 1;
                        }
                    }

                    order[dim] = o;
                }
                else if (rule[dim] == 15)
                {
                    o = 1;
                    while (o < 4 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }

                    order[dim] = o;
                }
                else if (rule[dim] == 16)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        p = 5;
                        o = 3;
                        while (p < 4 * level[dim] + 1)
                        {
                            p = 2 * p + 1;
                            o = 2 * o + 1;
                        }

                        order[dim] = o;
                    }
                }
                else if (rule[dim] == 17)
                {
                    order[dim] = 2 * level[dim] + 1;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_DEFAULT - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE["
                                      + dim + "] = " + rule[dim] + ".");
                    return;
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
                if (level[dim] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                }
                else if (rule[dim] == 1)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                    }
                }
                else if (rule[dim] == 2)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 3)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 4)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 5)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 6)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 7)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 8)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 9)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 10)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 11)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                    }
                }
                else if (rule[dim] == 12)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 13)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 14)
                {
                    if (level[dim] == 0)
                    {
                        order[dim] = 1;
                    }
                    else
                    {
                        order[dim] = (int) Math.Pow(2, level[dim]) + 1;
                    }
                }
                else if (rule[dim] == 15)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 16)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1) - 1;
                }
                else if (rule[dim] == 17)
                {
                    order[dim] = (int) Math.Pow(2, level[dim] + 1);
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL - Fatal error!");
                    Console.WriteLine("  Unexpected value of RULE["
                                      + dim + "] = " + rule[dim] + ".");
                    return;
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
            int o;
            int p;

            for (dim = 0; dim < dim_num; dim++)
            {
                if (level[dim] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXPONENTIAL_SLOW - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                }
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                if (rule[dim] == 1 || rule[dim] == 11 || rule[dim] == 14 || rule[dim] == 17)
                {
                    if (level[dim] == 0)
                    {
                        o = 1;
                    }
                    else
                    {
                        o = 2;
                        while (o < 2 * level[dim] + 1)
                        {
                            o = 2 * (o - 1) + 1;
                        }
                    }
                }
                else if (rule[dim] == 3 || rule[dim] == 13 || rule[dim] == 16)
                {
                    if (level[dim] == 0)
                    {
                        o = 1;
                    }
                    else
                    {
                        p = 5;
                        o = 3;
                        while (p < 2 * level[dim] + 1)
                        {
                            p = 2 * p + 1;
                            o = 2 * o + 1;
                        }
                    }
                }
                else
                {
                    o = 1;
                    while (o < 2 * level[dim] + 1)
                    {
                        o = 2 * o + 1;
                    }
                }

                order[dim] = o;
            }

            return;
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
                if (level[dim] < 0)
                {
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
            //
            //  Slow exponential growth.
            //
            if (growth == 0)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 2;
                    while (o < 2 * level + 1)
                    {
                        o = 2 * (o - 1) + 1;
                    }
                }
            }
            //
            //  Moderate Exponential Growth.
            //
            else if (growth == 1)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 2;
                    while (o < 4 * level + 1)
                    {
                        o = 2 * (o - 1) + 1;
                    }
                }
            }
            //
            //  Full Exponential Growth.
            //
            else if (growth == 2)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = (int) Math.Pow(2, level) + 1;
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_CC - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow exponential growth.
            //
            if (growth == 0)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 1;
                    while (o < 2 * level + 1)
                    {
                        o = 2 * o + 1;
                    }
                }
            }
            //
            //  Moderate Exponential Growth.
            //
            else if (growth == 1)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 1;
                    while (o < 4 * level + 1)
                    {
                        o = 2 * o + 1;
                    }
                }
            }
            //
            //  Full Exponential Growth.
            //
            else if (growth == 2)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = (int) Math.Pow(2, level + 1) - 1;
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_F2 - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow exponential growth.
            //
            if (growth == 0)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 1;
                    while (2 * o - 1 < 2 * level + 1)
                    {
                        o = 2 * o + 1;
                    }
                }
            }
            //
            //  Moderate Exponential Growth.
            //
            else if (growth == 1)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = 1;
                    while (2 * o - 1 < 4 * level + 1)
                    {
                        o = 2 * o + 1;
                    }
                }
            }
            //
            //  Full Exponential Growth.
            //
            else if (growth == 2)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = (int) Math.Pow(2, level + 1) - 1;
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_GAUSS - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow exponential growth.
            //
            if (growth == 0)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    p = 5;
                    o = 3;
                    while (p < 2 * level + 1)
                    {
                        p = 2 * p + 1;
                        o = 2 * o + 1;
                        if (511 < o)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                            Console.WriteLine("  Request for unavailable Patterson rule.");
                            return (1);
                        }
                    }
                }
            }
            //
            //  Moderate Exponential Growth.
            //
            else if (growth == 1)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    p = 5;
                    o = 3;
                    while (p < 4 * level + 1)
                    {
                        p = 2 * p + 1;
                        o = 2 * o + 1;
                        if (511 < o)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                            Console.WriteLine("  Request for unavailable Patterson rule.");
                            return (1);
                        }
                    }
                }
            }
            //
            //  Full Exponential Growth.
            //
            else if (growth == 2)
            {
                if (level == 0)
                {
                    o = 1;
                }
                else
                {
                    o = (int) Math.Pow(2, level + 1) - 1;
                    if (511 < o)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                        Console.WriteLine("  Request for unavailable Patterson rule.");
                        return (1);
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_GP - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow exponential growth.
            //
            if (growth == 0)
            {
                l = 0;
                p = p_hgk[l];
                o = o_hgk[l];
                while (p < 2 * level + 1)
                {
                    l = l + 1;
                    if (5 < l)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                        Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                        return (1);
                    }

                    p = p_hgk[l];
                    o = o_hgk[l];
                }
            }
            else if (growth == 1)
            {
                l = 0;
                p = p_hgk[l];
                o = o_hgk[l];
                while (p < 4 * level + 1)
                {
                    l = l + 1;
                    if (5 < l)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                        Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                        return (1);
                    }

                    p = p_hgk[l];
                    o = o_hgk[l];
                }
            }
            else if (growth == 2)
            {
                l = level;
                l = (int) Math.Max(l, 0);
                if (5 < l)
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                    Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                    return (1);
                }

                o = o_hgk[l];
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_EXP_HGK - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow linear growth.
            //
            if (growth == 0)
            {
                o = level + 1;
            }
            //
            //  Moderate linear growth.
            //
            else if (growth == 1)
            {
                o = 2 * level + 1;
            }
            else if (growth == 2)
            {
                o = 2 * level + 1;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_LINEAR_NN - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            //
            //  Slow growth.
            //
            if (growth == 0)
            {
                o = 2 * ((level + 1) / 2) + 1;
            }
            else if (growth == 1)
            {
                o = 2 * level + 1;
            }
            else if (growth == 2)
            {
                o = 2 * level + 1;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LEVEL_TO_ORDER_LINEAR_WN - Fatal error!");
                Console.WriteLine("  Illegal value of GROWTH = " + growth + "");
                return (1);
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
            int l;
            int o = 0;
            int[] o_hgk = {1, 3, 9, 19, 35, 43};
            int p;
            int[] p_hgk = {1, 5, 15, 29, 51, 67};
            //
            //  Check the input.
            //
            for (dim = 0; dim < dim_num; dim++)
            {
                if (level[dim] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Negative value of LEVEL[DIM]!");
                    Console.WriteLine("  LEVEL[" + dim + "] = " + level[dim] + "");
                    return;
                }

                if (rule[dim] < 1 || 12 < rule[dim])
                {
                    Console.WriteLine("");
                    Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                    Console.WriteLine("  Illegal value of RULE[DIM]!");
                    Console.WriteLine("  RULE[" + dim + "] = " + rule[dim] + "");
                    return;
                }

                if (growth[dim] < 0 || 6 < growth[dim])
                {
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
                //
                //  CC
                //  Default is Moderate Exponential Growth.
                //
                if (rule[dim] == 1)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = 2;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 5 || growth[dim] == 0)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = 2;
                            while (o < 4 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = (int) Math.Pow(2, level[dim]) + 1;
                        }
                    }
                }
                //
                //  F2
                //  Default is Moderate Exponential Growth.
                //
                else if (rule[dim] == 2)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (o < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5 || growth[dim] == 0)
                    {
                        o = 1;
                        while (o < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GP
                //  Default is Moderate Exponential Growth.
                //
                else if (rule[dim] == 3)
                {
                    if (growth[dim] == 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 1 for rule 3 not available!");
                        return;
                    }
                    else if (growth[dim] == 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 2 for rule 3 not available!");
                        return;
                    }
                    else if (growth[dim] == 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 3 for rule 3 not available!");
                        return;
                    }
                    else if (growth[dim] == 4)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            p = 5;
                            o = 3;
                            while (p < 2 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 5 || growth[dim] == 0)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            p = 5;
                            o = 3;
                            while (p < 4 * level[dim] + 1)
                            {
                                p = 2 * p + 1;
                                o = 2 * o + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GL
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 4)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GH
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 5)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GGH
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 6)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  LG
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 7)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GLG
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 8)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  GJ
                //  Default is Moderate Linear Growth.
                //
                else if (rule[dim] == 9)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  HGK
                //  Default is Moderate Exponential Growth.
                //  Exponential growth is interpreted to mean simply take successive rules.
                //
                else if (rule[dim] == 10)
                {
                    if (growth[dim] == 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 1 for rule 10 not available!");
                        return;
                    }
                    else if (growth[dim] == 2)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 2 for rule 10 not available!");
                        return;
                    }
                    else if (growth[dim] == 3)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                        Console.WriteLine("  Growth rate 3 for rule 10 not available!");
                        return;
                    }
                    else if (growth[dim] == 4)
                    {
                        l = 0;
                        p = p_hgk[l];
                        o = o_hgk[l];
                        while (p < 2 * level[dim] + 1)
                        {
                            l = l + 1;
                            if (5 < l)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                                Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                                return;
                            }

                            p = p_hgk[l];
                            o = o_hgk[l];
                        }
                    }
                    else if (growth[dim] == 5 || growth[dim] == 0)
                    {
                        l = 0;
                        p = p_hgk[l];
                        o = o_hgk[l];
                        while (p < 4 * level[dim] + 1)
                        {
                            l = l + 1;
                            if (5 < l)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                                Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                                return;
                            }

                            p = p_hgk[l];
                            o = o_hgk[l];
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        l = level[dim];
                        l = Math.Max(l, 0);
                        if (5 < l)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("LEVEL_GROWTH_TO_ORDER - Fatal error!");
                            Console.WriteLine("  Hermite Genz-Keister maximum level exceeded.");
                            return;
                        }

                        o = o_hgk[l];
                    }
                }
                //
                //  UO
                //  Default is Moderate Linear Growth.
                //  We assume the rule is of OPEN type and that it
                //  has a precision typical of Gauss rules.
                //
                else if (rule[dim] == 11)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        o = 1;
                        while (2 * o - 1 < 2 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        o = 1;
                        while (2 * o - 1 < 4 * level[dim] + 1)
                        {
                            o = 2 * o + 1;
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        o = (int) Math.Pow(2, level[dim] + 1) - 1;
                    }
                }
                //
                //  UC
                //  Default is Moderate Linear Growth.
                //  We assume the rule is of CLOSED type and that it
                //  has a precision typical of Clenshaw-Curtis rules.
                //
                else if (rule[dim] == 12)
                {
                    if (growth[dim] == 1)
                    {
                        o = level[dim] + 1;
                    }
                    else if (growth[dim] == 2)
                    {
                        o = 2 * ((level[dim] + 1) / 2) + 1;
                    }
                    else if (growth[dim] == 3 || growth[dim] == 0)
                    {
                        o = 2 * level[dim] + 1;
                    }
                    else if (growth[dim] == 4)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = 2;
                            while (o < 2 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 5)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = 2;
                            while (o < 4 * level[dim] + 1)
                            {
                                o = 2 * (o - 1) + 1;
                            }
                        }
                    }
                    else if (growth[dim] == 6)
                    {
                        if (level[dim] == 0)
                        {
                            o = 1;
                        }
                        else
                        {
                            o = (int) Math.Pow(2, level[dim]) + 1;
                        }
                    }
                }

                order[dim] = o;
            }
        }

    }
}