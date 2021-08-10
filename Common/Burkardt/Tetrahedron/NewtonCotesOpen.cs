using System;

namespace Burkardt.TetrahedronNS
{
    public static class NewtonCotesOpen
    {
        public static int tetrahedron_nco_degree(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_DEGREE: degree of an NCO rule for the tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Output, int TETRAHEDRON_NCO_DEGREE, the polynomial degree of exactness of
            //    the rule.
            //
        {
            int degree;

            if (1 <= rule && rule <= 7)
            {
                degree = rule - 1;
            }
            else
            {
                degree = -1;
                Console.WriteLine("");
                Console.WriteLine("TETRAHEDRON_NCO_DEGREE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return 1;
            }

            return degree;
        }

        public static int tetrahedron_nco_order_num(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_ORDER_NUM: order of an NCO rule for the tetrahedron.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Output, int TETRAHEDRON_NCO_ORDER_NUM, the order (number of points) 
            //    of the rule.
            //
        {
            int order;
            int order_num;
            int[] suborder;
            int suborder_num;

            suborder_num = tetrahedron_nco_suborder_num(rule);

            suborder = tetrahedron_nco_suborder(rule, suborder_num);

            order_num = 0;
            for (order = 0; order < suborder_num; order++)
            {
                order_num = order_num + suborder[order];
            }

            return order_num;
        }

        public static void tetrahedron_nco_rule(int rule, int order_num, ref double[] xyz, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_RULE returns the points and weights of an NCO rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Output, double XYZ[3*ORDER_NUM], the points of the rule.
            //
            //    Output, double W[ORDER_NUM], the weights of the rule.
            //
        {
            int o;
            int s;
            int[] suborder;
            int suborder_num;
            double[] suborder_w;
            double[] suborder_xyz;
            //
            //  Get the suborder information.
            //
            suborder_num = tetrahedron_nco_suborder_num(rule);

            suborder_xyz = new double[4 * suborder_num];
            suborder_w = new double[suborder_num];

            suborder = tetrahedron_nco_suborder(rule, suborder_num);

            tetrahedron_nco_subrule(rule, suborder_num, ref suborder_xyz, ref suborder_w);
            //
            //  Expand the suborder information to a full order rule.
            //
            o = 0;

            for (s = 0; s < suborder_num; s++)
            {
                if (suborder[s] == 1)
                {
                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                //
                //  Fourfold symmetry on (A,A,A,B)
                //
                //    123 AAA
                //    124 AAB
                //    142 ABA
                //    412 BAA
                //
                else if (suborder[s] == 4)
                {
                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                //
                //  Sixfold symmetry on (A,A,B,B):
                //
                //    123 (A,A,B)
                //    132 (A,B,A),
                //    134 (A,B,B)
                //    312 (B,A,A)
                //    314 (B,A,B)
                //    341 (B,B,A)
                //
                else if (suborder[s] == 6)
                {
                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                //
                //  Twelvefold symmetry on (A,A,B,C):
                //
                //    123 (A,A,B)
                //    124 (A,A,C)
                //    132 (A,B,A)
                //    134 (A,B,C)
                //    142 (A,C,A)
                //    143 (A,C,B)
                //    312 (B,A,A)
                //    314 (B,A,C)
                //    341 (B,C,A)
                //    412 (C,A,A)
                //    413 (C,A,B)
                //    431 (C,B,A)
                //
                else if (suborder[s] == 12)
                {
                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                //
                //  24 fold symmetry on (A,B,C,D):
                //
                //    123 (A,B,C)
                //    124 (A,B,D)
                //    132 (A,C,B)
                //    134 (A,C,D)
                //    142 (A,D,B)
                //    143 (A,D,C)
                //    213 (B,A,C)
                //    214 (B,A,D)
                //    231 (B,C,A)
                //    234 (B,C,D)
                //    241 (B,D,A)
                //    243 (B,D,C)
                //    312 (C,A,B)
                //    314 (C,A,D)
                //    321 (C,B,A)
                //    324 (C,B,D)
                //    341 (C,D,A)
                //    342 (C,D,B)
                //    412 (D,A,B)
                //    413 (D,A,C)
                //    421 (D,B,A)
                //    423 (D,B,C)
                //    431 (D,C,A)
                //    432 (D,C,B)
                //
                else if (suborder[s] == 24)
                {
                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[4 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[3 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[0 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[1 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[2 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[0 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;

                    xyz[0 + o * 3] = suborder_xyz[3 + s * 4];
                    xyz[1 + o * 3] = suborder_xyz[2 + s * 4];
                    xyz[2 + o * 3] = suborder_xyz[1 + s * 4];
                    w[o] = suborder_w[s];
                    o = o + 1;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("TETRAHEDRON_NCO_RULE - Fatal error!");
                    Console.WriteLine("  Illegal SUBORDER(" + s + ") = " + suborder[s] + "");
                    return;
                }
            }
        }

        public static int tetrahedron_nco_rule_num()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_RULE_NUM returns the number of NCO rules available.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Output, int TETRAHEDRON_NCO_RULE_NUM, the number of rules available.
            //
        {
            int rule_num;

            rule_num = 7;

            return rule_num;
        }

        public static int[] tetrahedron_nco_suborder(int rule, int suborder_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBORDER returns the suborders for an NCO rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int TETRAHEDRON_NCO_SUBORDER[SUBORDER_NUM], 
            //    the suborders of the rule.
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
                suborder[0] = 4;
                suborder[1] = 6;
            }
            else if (rule == 4)
            {
                suborder[0] = 4;
                suborder[1] = 12;
                suborder[2] = 4;
            }
            else if (rule == 5)
            {
                suborder[0] = 4;
                suborder[1] = 12;
                suborder[2] = 6;
                suborder[3] = 12;
                suborder[4] = 1;
            }
            else if (rule == 6)
            {
                suborder[0] = 4;
                suborder[1] = 12;
                suborder[2] = 12;
                suborder[3] = 12;
                suborder[4] = 12;
                suborder[5] = 4;
            }
            else if (rule == 7)
            {
                suborder[0] = 4;
                suborder[1] = 12;
                suborder[2] = 12;
                suborder[3] = 12;
                suborder[4] = 6;
                suborder[5] = 24;
                suborder[6] = 4;
                suborder[7] = 4;
                suborder[8] = 6;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TETRAHEDRON_NCO_SUBORDER - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return null;
            }

            return suborder;
        }

        public static int tetrahedron_nco_suborder_num(int rule)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBORDER_NUM returns the number of suborders for an NCO rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Output, int TETRAHEDRON_NCO_SUBORDER_NUM, the number of suborders 
            //    of the rule.
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
                suborder_num = 3;
            }
            else if (rule == 5)
            {
                suborder_num = 5;
            }
            else if (rule == 6)
            {
                suborder_num = 6;
            }
            else if (rule == 7)
            {
                suborder_num = 9;
            }
            else
            {
                suborder_num = -1;
                Console.WriteLine("");
                Console.WriteLine("TETRAHEDRON_NCO_SUBORDER_NUM - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return 1;
            }

            return suborder_num;
        }

        public static void tetrahedron_nco_subrule(int rule, int suborder_num,
                ref double[] suborder_xyz, ref double[] suborder_w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE returns a compressed NCO rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, double SUBORDER_XYZ[4*SUBORDER_NUM],
            //    the barycentric coordinates of the abscissas.
            //
            //    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
            //
        {
            int i;
            int s;
            int suborder_w_d = 0;
            int[] suborder_w_n;
            int suborder_xyz_d = 0;
            int[] suborder_xyz_n;

            suborder_xyz_n = new int[4 * suborder_num];
            suborder_w_n = new int[suborder_num];

            if (rule == 1)
            {
                tetrahedron_nco_subrule_01(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 2)
            {
                tetrahedron_nco_subrule_02(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 3)
            {
                tetrahedron_nco_subrule_03(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 4)
            {
                tetrahedron_nco_subrule_04(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 5)
            {
                tetrahedron_nco_subrule_05(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 6)
            {
                tetrahedron_nco_subrule_06(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else if (rule == 7)
            {
                tetrahedron_nco_subrule_07(suborder_num, ref suborder_xyz_n, ref suborder_xyz_d,
                    ref suborder_w_n, ref suborder_w_d);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TETRAHEDRON_NCO_SUBRULE - Fatal error!");
                Console.WriteLine("  Illegal RULE = " + rule + "");
                return;
            }

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz[i + s * 4] =
                        (double)(1 + suborder_xyz_n[i + s * 4])
                        / (double)(4 + suborder_xyz_d);
                }
            }

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w[s] = (double)suborder_w_n[s] / (double)suborder_w_d;
            }
        }

        public static void tetrahedron_nco_subrule_01(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_01 returns a compressed NCO rule 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_01 =
            {
                0, 0, 0, 0
            };
            int suborder_xyz_d_01 = 0;
            int[] suborder_w_n_01 = { 1 };
            int suborder_w_d_01 = 1;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_01[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_01;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_01[s];
            }

            suborder_w_d = suborder_w_d_01;

        }

        public static void tetrahedron_nco_subrule_02(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_02 returns a compressed NCO rule 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_02 =
            {
                0, 0, 0, 1
            };
            int suborder_xyz_d_02 = 1;
            int[] suborder_w_n_02 = { 1 };
            int suborder_w_d_02 = 4;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_02[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_02;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_02[s];
            }

            suborder_w_d = suborder_w_d_02;

        }

        public static void tetrahedron_nco_subrule_03(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_03 returns a compressed NCO rule 3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_03 =
            {
                0, 0, 0, 2,
                1, 1, 0, 0
            };
            int suborder_xyz_d_03 = 2;
            int[] suborder_w_n_03 = { 11, -4 };
            int suborder_w_d_03 = 20;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_03[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_03;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_03[s];
            }

            suborder_w_d = suborder_w_d_03;

        }

        public static void tetrahedron_nco_subrule_04(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_04 returns a compressed NCO rule 4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_04 =
            {
                0, 0, 0, 3,
                0, 0, 1, 2,
                1, 1, 1, 0
            };
            int suborder_xyz_d_04 = 3;
            int[] suborder_w_n_04 = { 20, 13, -29 };
            int suborder_w_d_04 = 120;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_04[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_04;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_04[s];
            }

            suborder_w_d = suborder_w_d_04;

        }

        public static void tetrahedron_nco_subrule_05(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_05 returns a compressed NCO rule 5.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_05 =
            {
                0, 0, 0, 4,
                0, 0, 3, 1,
                2, 2, 0, 0,
                1, 1, 0, 2,
                1, 1, 1, 1
            };
            int suborder_xyz_d_05 = 4;
            int[] suborder_w_n_05 = { 79, -68, 142, -12, 2 };
            int suborder_w_d_05 = 210;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_05[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_05;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_05[s];
            }

            suborder_w_d = suborder_w_d_05;

            return;
        }

        public static void tetrahedron_nco_subrule_06(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_06 returns a compressed NCO rule 6.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_06 =
            {
                0, 0, 0, 5,
                0, 0, 4, 1,
                0, 0, 3, 2,
                1, 1, 0, 3,
                2, 2, 1, 0,
                1, 1, 1, 2
            };
            int suborder_xyz_d_06 = 5;
            int[] suborder_w_n_06 = { 277, 97, 223, -713, 505, -53 };
            int suborder_w_d_06 = 2240;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_06[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_06;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_06[s];
            }

            suborder_w_d = suborder_w_d_06;

        }

        public static void tetrahedron_nco_subrule_07(int suborder_num, ref int[] suborder_xyz_n,
                ref int suborder_xyz_d, ref int[] suborder_w_n, ref int suborder_w_d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TETRAHEDRON_NCO_SUBRULE_07 returns a compressed NCO rule 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 January 2007
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
            //    Input, int SUBORDER_NUM, the number of suborders of the rule.
            //
            //    Output, int SUBORDER_XYZ_N[4*SUBORDER_NUM],
            //    the numerators of the barycentric coordinates of the abscissas.
            //
            //    Output, int *SUBORDER_XYZ_D,
            //    the denominator of the barycentric coordinates of the abscissas.
            //
            //    Output, int SUBORDER_W_N[SUBORDER_NUM], 
            //    the numerator of the suborder weights.
            //
            //    Output, int SUBORDER_W_D, 
            //    the denominator of the suborder weights.
            //
        {
            int i;
            int s;
            int[] suborder_xyz_n_07 =
            {
                0, 0, 0, 6,
                0, 0, 5, 1,
                0, 0, 4, 2,
                1, 1, 0, 4,
                3, 3, 0, 0,
                3, 2, 1, 0,
                1, 1, 1, 3,
                2, 2, 2, 0,
                2, 2, 1, 1
            };
            int suborder_xyz_d_07 = 6;
            int[] suborder_w_n_07 = { 430, -587, 1327, 187, -1298, -398, 22, 1537, -38 };
            int suborder_w_d_07 = 1512;

            for (s = 0; s < suborder_num; s++)
            {
                for (i = 0; i < 4; i++)
                {
                    suborder_xyz_n[i + s * 4] = suborder_xyz_n_07[i + s * 4];
                }
            }

            suborder_xyz_d = suborder_xyz_d_07;

            for (s = 0; s < suborder_num; s++)
            {
                suborder_w_n[s] = suborder_w_n_07[s];
            }

            suborder_w_d = suborder_w_d_07;

        }
    }
}