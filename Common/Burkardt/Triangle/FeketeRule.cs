using System;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static partial class FeketeRule
{
    public static int fekete_degree(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEKETE_DEGREE returns the degree of a Fekete rule for the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Mark Taylor, Beth Wingate, Rachel Vincent,
        //    An Algorithm for Computing Fekete Points in the Triangle,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 38, Number 5, 2000, pages 1707-1720.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int FEKETE_DEGREE, the polynomial degree of exactness of
        //    the rule.
        //
    {
        int degree;

        switch (rule)
        {
            case 1:
                degree = 3;
                break;
            case 2:
                degree = 6;
                break;
            case 3:
                degree = 9;
                break;
            case 4:
            case 5:
                degree = 12;
                break;
            case 6:
                degree = 15;
                break;
            case 7:
                degree = 18;
                break;
            default:
                degree = -1;
                Console.WriteLine("");
                Console.WriteLine("FEKETE_DEGREE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                break;
        }

        return degree;
    }

    public static int fekete_order_num(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEKETE_ORDER_NUM returns the order of a Fekete rule for the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Mark Taylor, Beth Wingate, Rachel Vincent,
        //    An Algorithm for Computing Fekete Points in the Triangle,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 38, Number 5, 2000, pages 1707-1720.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int FEKETE_ORDER_NUM, the order (number of points) of the rule.
        //
    {
        int order;

        int suborder_num = fekete_suborder_num(rule);

        int[] suborder = fekete_suborder(rule, suborder_num);

        int order_num = 0;
        for (order = 0; order < suborder_num; order++)
        {
            order_num += suborder[order];
        }

        return order_num;
    }

    public static void fekete_rule(int rule, int order_num, ref double[] xy, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEKETE_RULE returns the points and weights of a Fekete rule.
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
        //    Mark Taylor, Beth Wingate, Rachel Vincent,
        //    An Algorithm for Computing Fekete Points in the Triangle,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 38, Number 5, 2000, pages 1707-1720.
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
        int s;
        //
        //  Get the suborder information.
        //
        int suborder_num = fekete_suborder_num(rule);

        double[] suborder_xyz = new double[3 * suborder_num];
        double[] suborder_w = new double[suborder_num];

        int[] suborder = fekete_suborder(rule, suborder_num);

        fekete_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);
        //
        //  Expand the suborder information to a full order rule.
        //
        int o = 0;

        for (s = 0; s < suborder_num; s++)
        {
            int k;
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
                    Console.WriteLine("FEKETE_RULE - Fatal error!");
                    Console.WriteLine("  Illegal SUBORDER(" + s + ") = " + suborder[s] + "");
                    return;
            }
        }
    }

    public static int fekete_rule_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEKETE_RULE_NUM returns the number of Fekete rules available.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Mark Taylor, Beth Wingate, Rachel Vincent,
        //    An Algorithm for Computing Fekete Points in the Triangle,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 38, Number 5, 2000, pages 1707-1720.
        //
        //  Parameters:
        //
        //    Output, int FEKETE_RULE_NUM, the number of rules available.
        //
    {
        const int rule_num = 7;

        return rule_num;
    }
}