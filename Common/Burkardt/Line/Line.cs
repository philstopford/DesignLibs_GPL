﻿using System;
using Burkardt.Types;

namespace Burkardt.LineNS
{
    public class Line
    {
        public static bool line_exp_is_degenerate_nd(int dim_num, double[] p1, double[] p2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
            //
            //  Discussion:
            //
            //    The explicit form of a line in ND is:
            //
            //      the line through the points P1 and P2.
            //
            //    An explicit line is degenerate if the two defining points are equal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double P1[DIM_NUM], P2[DIM_NUM], two points on the line.
            //
            //    Output, bool LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
            //    is degenerate.
            //
        {
            bool value;

            value = typeMethods.r8vec_eq(dim_num, p1, p2);

            return value;
        }

        public static double[] line_exp_perp_2d(double[] p1, double[] p2, double[] p3,
                ref bool flag, int p1Index = 0, int p2Index = 0, int p3Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_EXP_PERP_2D computes a line perpendicular to a line and through a point.
            //
            //  Discussion:
            //
            //    The explicit form of a line in 2D is:
            //
            //      the line through P1 and P2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], two points on the given line.
            //
            //    Input, double P3[2], a point not on the given line, through which the
            //    perpendicular must pass.
            //
            //    Output, double LINE_EXP_PERP_2D[2], a point on the given line, such that the line
            //    through P3 and P4 is perpendicular to the given line.
            //
            //    Output, bool *FLAG, is TRUE if the point could not be computed.
            //
        {
            double bot;
            double[] p4;
            double t;

            p4 = new double[2];

            bot = Math.Pow(p2[p2Index + 0] - p1[p1Index + 0], 2) + Math.Pow(p2[p2Index + 1] - p1[p1Index + 1], 2);

            if (bot == 0.0)
            {
                p4[0] = typeMethods.r8_huge();
                p4[1] = typeMethods.r8_huge();
                flag = true;
                return p4;
            }

            //
            //  (P3-P1) dot (P2-P1) = Norm(P3-P1) * Norm(P2-P1) * Cos(Theta).
            //
            //  (P3-P1) dot (P2-P1) / Norm(P3-P1)**2 = normalized coordinate T
            //  of the projection of (P3-P1) onto (P2-P1).
            //
            t = ((p1[p1Index + 0] - p3[p3Index + 0]) * (p1[p1Index + 0] - p2[p2Index + 0])
                 + (p1[p1Index + 1] - p3[p3Index + 1]) * (p1[p1Index + 1] - p2[p2Index + 1])) / bot;

            p4[0] = p1[p1Index + 0] + t * (p2[p2Index + 0] - p1[p1Index + 0]);
            p4[1] = p1[p1Index + 1] + t * (p2[p2Index + 1] - p1[p1Index + 1]);

            flag = false;

            return p4;
        }

        public static void line_exp2imp_2d(double[] p1, double[] p2, ref double a, ref double b,
                ref double c, int p1Index = 0, int p2Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
            //
            //  Discussion:
            //
            //    The explicit form of a line in 2D is:
            //
            //      the line through P1 and P2
            //
            //    The implicit form of a line in 2D is:
            //
            //      A * X + B * Y + C = 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], two distinct points on the line. 
            //
            //    Output, double *A, *B, *C, three coefficients which describe
            //    the line that passes through P1 and P2.
            //
        {
            //
            //  Take care of degenerate cases.
            //
            if (typeMethods.r8vec_eq(2, p1, p2, p1Index, p2Index))
            {
                Console.WriteLine("");
                Console.WriteLine("LINE_EXP2IMP_2D - Fatal error!");
                Console.WriteLine("  P1 = P2");
                Console.WriteLine("  P1 = " + p1[p1Index + 0] + " " + p1[p1Index + 1] + "");
                Console.WriteLine("  P2 = " + p2[p2Index + 0] + " " + p2[p2Index + 1] + "");
                return;
            }

            a = p2[p2Index + 1] - p1[p1Index + 1];
            b = p1[p1Index + 0] - p2[p2Index + 0];
            c = p2[p2Index + 0] * p1[p1Index + 1] - p1[p1Index + 0] * p2[p2Index + 1];
        }

        public static bool line_imp_is_degenerate_2d(double a, double b, double c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
            //
            //  Discussion:
            //
            //    The implicit form of a line in 2D is:
            //
            //      A * X + B * Y + C = 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, C, the implicit line parameters.
            //
            //    Output, bool LINE_IMP_IS_DEGENERATE_2D, is true if the
            //    line is degenerate.
            //
        {
            bool value;

            value = (a * a + b * b == 0.0);

            return value;
        }

        public static void lines_exp_int_2d(double[] p1, double[] p2, double[] p3, double[] p4,
                ref int ival, ref double[] p, int p1Index = 0, int p2Index = 0, int p3Index = 0, int p4Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
            //
            //  Discussion:
            //
            //    The explicit form of a line in 2D is:
            //
            //      the line through P1 and P2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], define the first line.
            //
            //    Input, double P3[2], P4[2], define the second line.
            //
            //    Output, int *IVAL, reports on the intersection:
            //    0, no intersection, the lines may be parallel or degenerate.
            //    1, one intersection point, returned in P.
            //    2, infinitely many intersections, the lines are identical.
            //
            //    Output, double P[2], if IVAl = 1, then P contains
            //    the intersection point.  Otherwise, P = 0.
            //
        {
            double a1 = 0.0;
            double a2 = 0.0;
            double b1 = 0.0;
            double b2 = 0.0;
            double c1 = 0.0;
            double c2 = 0.0;
            bool point_1 = false;
            bool point_2 = false;

            ival = 0;
            p[0] = 0.0;
            p[1] = 0.0;
            //
            //  Check whether either line is a point.
            //
            if (typeMethods.r8vec_eq(2, p1, p2, p1Index, p2Index))
            {
                point_1 = true;
            }
            else
            {
                point_1 = false;
            }

            if (typeMethods.r8vec_eq(2, p3, p4, p3Index, p4Index))
            {
                point_2 = true;
            }
            else
            {
                point_2 = false;
            }

            //
            //  Convert the lines to ABC format.
            //
            if (!point_1)
            {
                line_exp2imp_2d(p1, p2, ref a1, ref b1, ref c1, p1Index, p2Index);
            }

            if (!point_2)
            {
                line_exp2imp_2d(p3, p4, ref a2, ref b2, ref c2, p3Index, p4Index);
            }

            //
            //  Search for intersection of the lines.
            //
            if (point_1 && point_2)
            {
                if (typeMethods.r8vec_eq(2, p1, p3, p1Index, p3Index))
                {
                    ival = 1;
                    typeMethods.r8vec_copy(2, p1, ref p, p1Index);
                }
            }
            else if (point_1)
            {
                if (a2 * p1[0] + b2 * p1[1] == c2)
                {
                    ival = 1;
                    typeMethods.r8vec_copy(2, p1, ref p, p1Index);
                }
            }
            else if (point_2)
            {
                if (a1 * p3[0] + b1 * p3[1] == c1)
                {
                    ival = 1;
                    typeMethods.r8vec_copy(2, p3, ref p, p3Index);
                }
            }
            else
            {
                lines_imp_int_2d(a1, b1, c1, a2, b2, c2, ref ival, ref p);
            }
        }

        public static void lines_imp_int_2d(double a1, double b1, double c1, double a2, double b2,
                double c2, ref int ival, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
            //
            //  Discussion:
            //
            //    The implicit form of a line in 2D is:
            //
            //      A * X + B * Y + C = 0
            //
            //    22 May 2004: Thanks to John Asmuth for pointing out that the 
            //    B array was not being deallocated on exit.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A1, B1, C1, define the first line.
            //    At least one of A1 and B1 must be nonzero.
            //
            //    Input, double A2, B2, C2, define the second line.
            //    At least one of A2 and B2 must be nonzero.
            //
            //    Output, int *IVAL, reports on the intersection.
            //    -1, both A1 and B1 were zero.
            //    -2, both A2 and B2 were zero.
            //     0, no intersection, the lines are parallel.
            //     1, one intersection point, returned in P.
            //     2, infinitely many intersections, the lines are identical.
            //
            //    Output, double P[2], if IVAL = 1, then P contains
            //    the intersection point.  Otherwise, P = 0.
            //
        {
            double[] a = new double[2 * 2];
            double[] b;

            p[0] = 0.0;
            p[1] = 0.0;
            //
            //  Refuse to handle degenerate lines.
            //
            if (a1 == 0.0 && b1 == 0.0)
            {
                ival = -1;
                return;
            }
            else if (a2 == 0.0 && b2 == 0.0)
            {
                ival = -2;
                return;
            }

            //
            //  Set up a linear system, and compute its inverse.
            //
            a[0 + 0 * 2] = a1;
            a[0 + 1 * 2] = b1;
            a[1 + 0 * 2] = a2;
            a[1 + 1 * 2] = b2;

            b = typeMethods.r8mat_inverse_2d(a);
            //
            //  If the inverse exists, then the lines intersect.
            //  Multiply the inverse times -C to get the intersection point.
            //
            if (b != null)
            {

                ival = 1;
                p[0] = -b[0 + 0 * 2] * c1 - b[0 + 1 * 2] * c2;
                p[1] = -b[1 + 0 * 2] * c1 - b[1 + 1 * 2] * c2;
            }
            //
            //  If the inverse does not exist, then the lines are parallel
            //  or coincident.  Check for parallelism by seeing if the
            //  C entries are in the same ratio as the A or B entries.
            //
            else
            {

                ival = 0;

                if (a1 == 0.0)
                {
                    if (b2 * c1 == c2 * b1)
                    {
                        ival = 2;
                    }
                }
                else
                {
                    if (a2 * c1 == c2 * a1)
                    {
                        ival = 2;
                    }
                }
            }
        }
    }
}