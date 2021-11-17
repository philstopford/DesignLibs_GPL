using System;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static partial class LynessRule
{
    public static int lyness_order(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_ORDER returns the order of a Lyness quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int LYNESS_ORDER, the order of the rule.
        //
    {
        int order;

        switch (rule)
        {
            case 0:
                order = 1;
                break;
            case 1:
                order = 3;
                break;
            case 2:
            case 3:
                order = 4;
                break;
            case 4:
                order = 7;
                break;
            case 5:
                order = 6;
                break;
            case 6:
                order = 10;
                break;
            case 7:
                order = 9;
                break;
            case 8:
                order = 7;
                break;
            case 9:
                order = 10;
                break;
            case 10:
                order = 12;
                break;
            case 11:
                order = 16;
                break;
            case 12:
            case 13:
                order = 13;
                break;
            case 14:
            case 15:
                order = 16;
                break;
            case 16:
                order = 21;
                break;
            case 17:
                order = 16;
                break;
            case 18:
                order = 19;
                break;
            case 19:
                order = 22;
                break;
            case 20:
                order = 27;
                break;
            case 21:
                order = 28;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LYNESS_ORDER - Fatal error!");
                Console.WriteLine("  Unrecognized rule index.");
                return 1;
        }

        return order;
    }

    public static int lyness_precision(int rule)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_PRECISION returns the precision of a Lyness quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Output, int LYNESS_PRECISION, the precision of the rule.
        //
    {
        int precision;

        switch (rule)
        {
            case 0:
                precision = 1;
                break;
            case 1:
            case 2:
                precision = 2;
                break;
            case 3:
            case 4:
                precision = 3;
                break;
            case 5:
            case 6:
            case 7:
                precision = 4;
                break;
            case 8:
            case 9:
                precision = 5;
                break;
            case 10:
            case 11:
            case 12:
                precision = 6;
                break;
            case 13:
            case 14:
                precision = 7;
                break;
            case 15:
            case 16:
            case 17:
                precision = 8;
                break;
            case 18:
            case 19:
                precision = 9;
                break;
            case 20:
            case 21:
                precision = 11;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LYNESS_PRECISION - Fatal error!");
                Console.WriteLine("  Unrecognized rule index.");
                return 1;
        }

        return precision;
    }

    public static void lyness_rule(int rule, int order, ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_RULE returns the points and weights of a Lyness quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Input, int RULE, the index of the rule.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Output, double W[ORDER], the weights.
        //
        //    Output, double X[2*ORDER], the points.
        //
    {
        int k;
        int o;
        int s;
        int[] suborder;
        int suborder_num;
        double[] sub_w;
        double[] sub_xyz;
        //
        //  Get the suborder information.
        //
        suborder_num = lyness_suborder_num(rule);

        suborder = lyness_suborder(rule, suborder_num);

        sub_xyz = new double[3 * suborder_num];
        sub_w = new double[suborder_num];

        lyness_subrule(rule, suborder_num, ref sub_xyz, ref sub_w);
        //
        //  Expand the suborder information to a full order rule.
        //
        o = 0;

        for (s = 0; s < suborder_num; s++)
        {
            switch (suborder[s])
            {
                case 1:
                    x[0 + o * 2] = sub_xyz[0 + s * 3];
                    x[1 + o * 2] = sub_xyz[1 + s * 3];
                    w[o] = sub_w[s];
                    o += 1;
                    break;
                case 3:
                {
                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 3.0;
                        o += 1;
                    }

                    break;
                }
                case 6:
                {
                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 6.0;
                        o += 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 6.0;
                        o += 1;
                    }

                    break;
                }
                default:
                    Console.WriteLine("");
                    Console.WriteLine("LYNESS_RULE - Fatal error!");
                    Console.WriteLine("  Illegal SUBORDER[" + s + "] = " + suborder[2] + "");
                    return;
            }
        }
    }

    public static int lyness_rule_num()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LYNESS_RULE_NUM returns the number of Lyness quadrature rules.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    James Lyness, Dennis Jespersen,
        //    Moderate Degree Symmetric Quadrature Rules for the Triangle,
        //    Journal of the Institute of Mathematics and its Applications,
        //    Volume 15, Number 1, February 1975, pages 19-32.
        //
        //  Parameters:
        //
        //    Output, int LYNESS_RULE_NUM, the number of rules.
        //
    {
        int rule_num;

        rule_num = 21;

        return rule_num;
    }
}