using System;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static partial class NewtonCotesClosed
{
    public static int triangle_ncc_degree(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCC_DEGREE returns the degree of an NCC rule for the triangle.
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
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int TRIANGLE_NCC_DEGREE, the polynomial degree of exactness of
        //    the rule.
        //
    {
        int degree;

        switch (rule)
        {
            case >= 1 and <= 9:
                degree = rule - 1;
                break;
            default:
                degree = -1;
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_NCC_DEGREE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return 1;
        }

        return degree;
    }

    public static int triangle_ncc_order_num(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCC_ORDER_NUM returns the order of an NCC rule for the triangle.
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
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int TRIANGLE_NCC_ORDER_NUM, the order (number of points)
        //    of the rule.
        //
    {
        int order;
        int order_num;
        int[] suborder;
        int suborder_num;

        suborder_num = triangle_ncc_suborder_num(rule);

        suborder = triangle_ncc_suborder(rule, suborder_num);

        order_num = 0;
        for (order = 0; order < suborder_num; order++)
        {
            order_num += suborder[order];
        }

        return order_num;
    }

    public static void triangle_ncc_rule(int rule, int order_num, ref double[] xy, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCC_RULE returns the points and weights of an NCC rule.
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
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int ORDER_NUM, the order (number of points) of the rule.
        //
        //    Output, double XY[2*ORDER_NUM], the points of the rule.
        //
        //    Output, double W[ORDER_NUM], the weights of the rule.
        //
    {
        int k;
        int o;
        int s;
        int[] suborder;
        int suborder_num;
        double[] suborder_w;
        double[] suborder_xyz;
        //
        //  Get the suborder information.
        //
        suborder_num = triangle_ncc_suborder_num(rule);

        suborder_xyz = new double[3 * suborder_num];
        suborder_w = new double[suborder_num];

        suborder = triangle_ncc_suborder(rule, suborder_num);

        triangle_ncc_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);
        //
        //  Expand the suborder information to a full order rule.
        //
        o = 0;

        for (s = 0; s < suborder_num; s++)
        {
            switch (suborder[s])
            {
                case 1:
                    xy[0 + o * 2] = suborder_xyz[0 + s * 3];
                    xy[1 + o * 2] = suborder_xyz[1 + s * 3];
                    w[o] = suborder_w[s];
                    o += 1;
                    break;
                case 3:
                {
                    for (k = 0; k < 3; k++)
                    {
                        xy[0 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        xy[1 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = suborder_w[s];
                        o += 1;
                    }

                    break;
                }
                case 6:
                {
                    for (k = 0; k < 3; k++)
                    {
                        xy[0 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        xy[1 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = suborder_w[s];
                        o += 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        xy[0 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        xy[1 + o * 2] = suborder_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        w[o] = suborder_w[s];
                        o += 1;
                    }

                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGLE_NCC_RULE - Fatal error!");
                    Console.WriteLine("  Illegal SUBORDER(" + s + ") = " + suborder[s] + "");
                    return;
            }
        }
    }

    public static int triangle_ncc_rule_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NCC_RULE_NUM returns the number of NCC rules available.
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
        //  Reference:
        //
        //    Peter Silvester,
        //    Symmetric Quadrature Formulae for Simplexes,
        //    Mathematics of Computation,
        //    Volume 24, Number 109, January 1970, pages 95-100.
        //
        //  Parameters:
        //
        //    Output, int TRIANGLE_NCC_RULE_NUM, the number of rules available.
        //
    {
        int rule_num;

        rule_num = 9;

        return rule_num;
    }
}