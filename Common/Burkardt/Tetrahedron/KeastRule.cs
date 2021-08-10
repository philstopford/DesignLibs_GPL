using System;
using Burkardt.Types;

namespace Burkardt.TetrahedronNS
{
    public static class KeastRule
    {
        public static int keast_degree(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_DEGREE returns the degree of a Keast rule for the triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Output, int KEAST_DEGREE, the polynomial degree of exactness of
            //    the rule.
            //
        {
            int degree;

            if (rule == 1)
            {
                degree = 0;
            }
            else if (rule == 2)
            {
                degree = 1;
            }
            else if (rule == 3)
            {
                degree = 2;
            }
            else if (rule == 4)
            {
                degree = 3;
            }
            else if (rule == 5)
            {
                degree = 4;
            }
            else if (rule == 6)
            {
                degree = 4;
            }
            else if (rule == 7)
            {
                degree = 5;
            }
            else if (rule == 8)
            {
                degree = 6;
            }
            else if (rule == 9)
            {
                degree = 7;
            }
            else if (rule == 10)
            {
                degree = 8;
            }
            else
            {
                degree = -1;
                Console.WriteLine("");
                Console.WriteLine("KEAST_DEGREE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return (1);
            }

            return degree;
        }

        public static int keast_order_num(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_ORDER_NUM returns the order of a Keast rule for the triangle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Output, int KEAST_ORDER_NUM, the order (number of points) of the rule.
            //
        {
            int order;
            int order_num;
            int[] suborder;
            int suborder_num;

            suborder_num = keast_suborder_num(rule);

            suborder = keast_suborder(rule, suborder_num);

            order_num = 0;
            for (order = 0; order < suborder_num; order++)
            {
                order_num = order_num + suborder[order];
            }

            return order_num;
        }

        public static void keast_rule(int rule, int order_num, ref double[] xyz, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_RULE returns the points and weights of a Keast rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Input, int ORDER_NUM, the order (number of points) of the rule.
            //
            //    Output, double XYZ[3*ORDER_NUM], the points of the rule.
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
            double[] suborder_xyzz;
            //
            //  Get the suborder information.
            //
            suborder_num = keast_suborder_num(rule);

            suborder_xyzz = new double[4 * suborder_num];
            suborder_w = new double[suborder_num];

            suborder = keast_suborder(rule, suborder_num);

            keast_subrule(rule, suborder_num, ref suborder_xyzz, ref suborder_w);
            //
            //  Expand the suborder information to a full order rule.
            //
            o = 0;

            for (s = 0; s < suborder_num; s++)
            {
                if (suborder[s] == 1)
                {
                    xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyzz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyzz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                //
                //  For SUBORDER = 4, we list the coordinates of the generator as
                //
                //    A,B,B,B
                //
                //  and we generate
                //
                //    A, B, B = (1,2,3)
                //    B, B, B = (2,3,4)
                //    B, B, A = (3,4,1)
                //    B, A, B = (4,1,2)
                //
                else if (suborder[s] == 4)
                {
                    for (k = 0; k < 4; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k, 0, 3) + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 1, 0, 3) + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 2, 0, 3) + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }
                }
                //
                //  For SUBORDER = 6, we list the coordinates of the generator as
                //
                //    A,A,B,B
                //
                //  and we generate
                //
                //    B, A, A = (4,1,2)
                //    A, B, A = (1,4,2)
                //    A, A, B = (1,2,4)
                //
                //    A, B, B = (1,3,4)
                //    B, A, B = (4,2,3)
                //    B, B, A = (4,3,1)
                //
                else if (suborder[s] == 6)
                {
                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[k + o * 3] = suborder_xyzz[2 + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[2 + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[2 + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[2 + s * 4];
                        xyz[k + o * 3] = suborder_xyzz[0 + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }
                }
                //
                //  For SUBORDER = 12, we list the coordinates of the generator as
                //
                //    A,A,B,C
                //
                //  and we generate
                //
                //    B, A, A
                //    A, B, A
                //    A, A, B
                //
                //    C, A, A
                //    A, C, A
                //    A, A, C
                //
                //    A, B, C
                //    B, C, A
                //    C, A, B
                //    A, C, B
                //    C, B, A
                //    B, A, C
                //
                else if (suborder[s] == 12)
                {
                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[k + o * 3] = suborder_xyzz[2 + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
                        xyz[k + o * 3] = suborder_xyzz[3 + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 1, 1, 3) + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 2, 1, 3) + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 3, 1, 3) + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }

                    for (k = 0; k < 3; k++)
                    {
                        xyz[0 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 1, 1, 3) + s * 4];
                        xyz[1 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 3, 1, 3) + s * 4];
                        xyz[2 + o * 3] = suborder_xyzz[typeMethods.i4_wrap(k + 2, 1, 3) + s * 4];
                        w[o] = suborder_w[s];
                        o = o + 1;
                    }
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("KEAST_RULE - Fatal error!");
                    Console.WriteLine("  Illegal SUBORDER(" + s + ") = " + suborder[s] + "");
                }
            }
        }

        public static int keast_rule_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_RULE_NUM returns the number of Keast rules available.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Output, int KEAST_RULE_NUM, the number of rules available.
            //
        {
            int rule_num;

            rule_num = 10;

            return rule_num;
        }

        public static int[] keast_suborder(int rule, int suborder_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBORDER returns the suborders for a Keast rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int KEAST_SUBORDER[SUBORDER_NUM], the suborders of the rule.
            //
        {
            int[] suborder;

            suborder = new int[suborder_num];

            if (rule == 1)
            {
                suborder[0] = 1;
            }
            else if (rule == 2)
            {
                suborder[0] = 4;
            }
            else if (rule == 3)
            {
                suborder[0] = 1;
                suborder[1] = 4;
            }
            else if (rule == 4)
            {
                suborder[0] = 4;
                suborder[1] = 6;
            }
            else if (rule == 5)
            {
                suborder[0] = 1;
                suborder[1] = 4;
                suborder[2] = 6;
            }
            else if (rule == 6)
            {
                suborder[0] = 6;
                suborder[1] = 4;
                suborder[2] = 4;
            }
            else if (rule == 7)
            {
                suborder[0] = 1;
                suborder[1] = 4;
                suborder[2] = 4;
                suborder[3] = 6;
            }
            else if (rule == 8)
            {
                suborder[0] = 4;
                suborder[1] = 4;
                suborder[2] = 4;
                suborder[3] = 12;
            }
            else if (rule == 9)
            {
                suborder[0] = 1;
                suborder[1] = 4;
                suborder[2] = 4;
                suborder[3] = 4;
                suborder[4] = 6;
                suborder[5] = 12;
            }
            else if (rule == 10)
            {
                suborder[0] = 1;
                suborder[1] = 4;
                suborder[2] = 4;
                suborder[3] = 6;
                suborder[4] = 6;
                suborder[5] = 12;
                suborder[6] = 12;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("KEAST_SUBORDER - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return (null);
            }

            return suborder;
        }

        public static int keast_suborder_num(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBORDER_NUM returns the number of suborders for a Keast rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Output, int KEAST_SUBORDER_NUM, the number of suborders of the rule.
            //
        {
            int suborder_num;

            if (rule == 1)
            {
                suborder_num = 1;
            }
            else if (rule == 2)
            {
                suborder_num = 1;
            }
            else if (rule == 3)
            {
                suborder_num = 2;
            }
            else if (rule == 4)
            {
                suborder_num = 2;
            }
            else if (rule == 5)
            {
                suborder_num = 3;
            }
            else if (rule == 6)
            {
                suborder_num = 3;
            }
            else if (rule == 7)
            {
                suborder_num = 4;
            }
            else if (rule == 8)
            {
                suborder_num = 4;
            }
            else if (rule == 9)
            {
                suborder_num = 6;
            }
            else if (rule == 10)
            {
                suborder_num = 7;
            }
            else
            {
                suborder_num = -1;
                Console.WriteLine("");
                Console.WriteLine("KEAST_SUBORDER_NUM - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return (1);
            }

            return suborder_num;
        }

        public static void keast_subrule(int rule, int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE returns a compressed Keast rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int RULE, the index of the rule.
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;

            if (rule == 1)
            {
                keast_subrule_01(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 2)
            {
                keast_subrule_02(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 3)
            {
                keast_subrule_03(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 4)
            {
                keast_subrule_04(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 5)
            {
                keast_subrule_05(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 6)
            {
                keast_subrule_06(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 7)
            {
                keast_subrule_07(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 8)
            {
                keast_subrule_08(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 9)
            {
                keast_subrule_09(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else if (rule == 10)
            {
                keast_subrule_10(suborder_num, ref suborder_xyzz, ref suborder_w);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("KEAST_SUBRULE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return;
            }

            //
            //  Renormalize the weights so they sum to 1.
            //
            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = 6.0 * suborder_w[s];
            }
        }

        public static void keast_subrule_01(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_01 returns a compressed Keast rule 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348..
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_01 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000
            };
            double[] suborder_w_rule_01 =
            {
                0.166666666666666667
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_01[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_01[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_01[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_01[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_01[s];
            }
        }

        public static void keast_subrule_02(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_02 returns a compressed Keast rule 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_02 =
            {
                0.585410196624968500, 0.138196601125010500,
                0.138196601125010500, 0.138196601125010500
            };
            double[] suborder_w_rule_02 =
            {
                0.0416666666666666667
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_02[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_02[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_02[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_02[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_02[s];
            }
        }

        public static void keast_subrule_03(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_03 returns a compressed Keast rule 3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_03 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000,
                0.500000000000000000, 0.166666666666666667,
                0.166666666666666667, 0.166666666666666667
            };
            double[] suborder_w_rule_03 =
            {
                -0.133333333333333333,
                0.075000000000000000
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_03[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_03[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_03[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_03[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_03[s];
            }
        }

        public static void keast_subrule_04(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_04 returns a compressed Keast rule 4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_04 =
            {
                0.568430584196844400, 0.143856471934385200,
                0.143856471934385200, 0.143856471934385200,
                0.500000000000000000, 0.500000000000000000,
                0.000000000000000000, 0.000000000000000000
            };
            double[] suborder_w_rule_04 =
            {
                0.0362941783134009000,
                0.00358165890217718333
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_04[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_04[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_04[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_04[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_04[s];
            }

        }

        public static void keast_subrule_05(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_05 returns a compressed Keast rule 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_05 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000,
                0.785714285714285714, 0.0714285714285714285,
                0.0714285714285714285, 0.0714285714285714285,
                0.399403576166799219, 0.399403576166799219,
                0.100596423833200785, 0.100596423833200785
            };
            double[] suborder_w_rule_05 =
            {
                -0.0131555555555555556,
                0.00762222222222222222,
                0.0248888888888888889
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_05[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_05[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_05[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_05[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_05[s];
            }

        }

        public static void keast_subrule_06(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_06 returns a compressed Keast rule 6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_06 =
            {
                0.500000000000000000, 0.500000000000000000,
                0.000000000000000000, 0.000000000000000000,
                0.698419704324386603, 0.100526765225204467,
                0.100526765225204467, 0.100526765225204467,
                0.0568813795204234229, 0.314372873493192195,
                0.314372873493192195, 0.314372873493192195
            };
            double[] suborder_w_rule_06 =
            {
                0.00317460317460317450,
                0.0147649707904967828,
                0.0221397911142651221
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_06[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_06[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_06[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_06[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_06[s];
            }

        }

        public static void keast_subrule_07(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_07 returns a compressed Keast rule 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_07 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000,
                0.00000000000000000, 0.333333333333333333,
                0.333333333333333333, 0.333333333333333333,
                0.727272727272727273, 0.0909090909090909091,
                0.0909090909090909091, 0.0909090909090909091,
                0.0665501535736642813, 0.0665501535736642813,
                0.433449846426335728, 0.433449846426335728
            };
            double[] suborder_w_rule_07 =
            {
                0.0302836780970891856,
                0.00602678571428571597,
                0.0116452490860289742,
                0.0109491415613864534
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_07[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_07[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_07[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_07[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_07[s];
            }

        }

        public static void keast_subrule_08(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_08 returns a compressed Keast rule 8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_08 =
            {
                0.356191386222544953, 0.214602871259151684,
                0.214602871259151684, 0.214602871259151684,
                0.877978124396165982, 0.0406739585346113397,
                0.0406739585346113397, 0.0406739585346113397,
                0.0329863295731730594, 0.322337890142275646,
                0.322337890142275646, 0.322337890142275646,
                0.0636610018750175299, 0.0636610018750175299,
                0.269672331458315867, 0.603005664791649076
            };
            double[] suborder_w_rule_08 =
            {
                0.00665379170969464506,
                0.00167953517588677620,
                0.00922619692394239843,
                0.00803571428571428248
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_08[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_08[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_08[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_08[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_08[s];
            }

        }

        public static void keast_subrule_09(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_08 returns a compressed Keast rule 8.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_09 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000,
                0.765360423009044044, 0.0782131923303186549,
                0.0782131923303186549, 0.0782131923303186549,
                0.634470350008286765, 0.121843216663904411,
                0.121843216663904411, 0.121843216663904411,
                0.00238250666073834549, 0.332539164446420554,
                0.332539164446420554, 0.332539164446420554,
                0.500000000000000000, 0.500000000000000000,
                0.00000000000000000, 0.00000000000000000,
                0.100000000000000000, 0.100000000000000000,
                0.200000000000000000, 0.600000000000000000
            };
            double[] suborder_w_rule_09 =
            {
                0.0182642234661087939,
                0.0105999415244141609,
                -0.0625177401143299494,
                0.00489142526307353653,
                0.000970017636684296702,
                0.0275573192239850917
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_09[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_09[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_09[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_09[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_09[s];
            }

            return;
        }

        public static void keast_subrule_10(int suborder_num, ref double[] suborder_xyzz,
                ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEAST_SUBRULE_10 returns a compressed Keast rule 10.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Patrick Keast,
            //    Moderate Degree Tetrahedral Quadrature Formulas,
            //    Computer Methods in Applied Mechanics and Engineering,
            //    Volume 55, Number 3, May 1986, pages 339-348.
            //
            //  Parameters:
            //
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int s;
            double[] suborder_xyzz_rule_10 =
            {
                0.250000000000000000, 0.250000000000000000,
                0.250000000000000000, 0.250000000000000000,
                0.617587190300082967, 0.127470936566639015,
                0.127470936566639015, 0.127470936566639015,
                0.903763508822103123, 0.0320788303926322960,
                0.0320788303926322960, 0.0320788303926322960,
                0.0497770956432810185, 0.0497770956432810185,
                0.450222904356718978, 0.450222904356718978,
                0.183730447398549945, 0.183730447398549945,
                0.316269552601450060, 0.316269552601450060,
                0.231901089397150906, 0.231901089397150906,
                0.0229177878448171174, 0.513280033360881072,
                0.0379700484718286102, 0.0379700484718286102,
                0.730313427807538396, 0.193746475248804382
            };
            double[] suborder_w_rule_10 =
            {
                -0.0393270066412926145,
                0.00408131605934270525,
                0.000658086773304341943,
                0.00438425882512284693,
                0.0138300638425098166,
                0.00424043742468372453,
                0.00223873973961420164
            };

            for (s = 0; s < suborder_num; s++)
            {
                suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_10[0 + s * 4];
                suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_10[1 + s * 4];
                suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_10[2 + s * 4];
                suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_10[3 + s * 4];
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = suborder_w_rule_10[s];
            }

        }
    }
}