using System;
using Burkardt.Types;

namespace Burkardt.TriangleNS
{
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

            if (rule == 0)
            {
                order = 1;
            }
            else if (rule == 1)
            {
                order = 3;
            }
            else if (rule == 2)
            {
                order = 4;
            }
            else if (rule == 3)
            {
                order = 4;
            }
            else if (rule == 4)
            {
                order = 7;
            }
            else if (rule == 5)
            {
                order = 6;
            }
            else if (rule == 6)
            {
                order = 10;
            }
            else if (rule == 7)
            {
                order = 9;
            }
            else if (rule == 8)
            {
                order = 7;
            }
            else if (rule == 9)
            {
                order = 10;
            }
            else if (rule == 10)
            {
                order = 12;
            }
            else if (rule == 11)
            {
                order = 16;
            }
            else if (rule == 12)
            {
                order = 13;
            }
            else if (rule == 13)
            {
                order = 13;
            }
            else if (rule == 14)
            {
                order = 16;
            }
            else if (rule == 15)
            {
                order = 16;
            }
            else if (rule == 16)
            {
                order = 21;
            }
            else if (rule == 17)
            {
                order = 16;
            }
            else if (rule == 18)
            {
                order = 19;
            }
            else if (rule == 19)
            {
                order = 22;
            }
            else if (rule == 20)
            {
                order = 27;
            }
            else if (rule == 21)
            {
                order = 28;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LYNESS_ORDER - Fatal error!");
                Console.WriteLine("  Unrecognized rule index.");
                return (1);
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

            if (rule == 0)
            {
                precision = 1;
            }
            else if (rule == 1)
            {
                precision = 2;
            }
            else if (rule == 2)
            {
                precision = 2;
            }
            else if (rule == 3)
            {
                precision = 3;
            }
            else if (rule == 4)
            {
                precision = 3;
            }
            else if (rule == 5)
            {
                precision = 4;
            }
            else if (rule == 6)
            {
                precision = 4;
            }
            else if (rule == 7)
            {
                precision = 4;
            }
            else if (rule == 8)
            {
                precision = 5;
            }
            else if (rule == 9)
            {
                precision = 5;
            }
            else if (rule == 10)
            {
                precision = 6;
            }
            else if (rule == 11)
            {
                precision = 6;
            }
            else if (rule == 12)
            {
                precision = 6;
            }
            else if (rule == 13)
            {
                precision = 7;
            }
            else if (rule == 14)
            {
                precision = 7;
            }
            else if (rule == 15)
            {
                precision = 8;
            }
            else if (rule == 16)
            {
                precision = 8;
            }
            else if (rule == 17)
            {
                precision = 8;
            }
            else if (rule == 18)
            {
                precision = 9;
            }
            else if (rule == 19)
            {
                precision = 9;
            }
            else if (rule == 20)
            {
                precision = 11;
            }
            else if (rule == 21)
            {
                precision = 11;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LYNESS_PRECISION - Fatal error!");
                Console.WriteLine("  Unrecognized rule index.");
                return (1);
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
                if (suborder[s] == 1)
                {
                    x[0 + o * 2] = sub_xyz[0 + s * 3];
                    x[1 + o * 2] = sub_xyz[1 + s * 3];
                    w[o] = sub_w[s];
                    o = o + 1;
                }
                else if (suborder[s] == 3)
                {
                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 3.0;
                        o = o + 1;
                    }
                }
                else if (suborder[s] == 6)
                {
                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 6.0;
                        o = o + 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        x[0 + o * 2] = sub_xyz[typeMethods.i4_wrap(k + 1, 0, 2) + s * 3];
                        x[1 + o * 2] = sub_xyz[typeMethods.i4_wrap(k, 0, 2) + s * 3];
                        w[o] = sub_w[s] / 6.0;
                        o = o + 1;
                    }
                }
                else
                {
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
}