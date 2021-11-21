using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.CircleNS;

public static class Geometry
{
    public static void circle_arc_point_near_2d(double r, double[] pc, double theta1,
            double theta2, double[] p, double[] pn, ref double dist)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_ARC_POINT_NEAR_2D : nearest point on a circular arc.
        //
        //  Discussion:
        //
        //    A circular arc is defined by the portion of a circle (R,PC)
        //    between two angles (THETA1,THETA2).
        //
        //    A point P on a circular arc satisfies
        //
        //      ( X - PC(1) ) * ( X - PC(1) ) + ( Y - PC(2) ) * ( Y - PC(2) ) = R * R
        //
        //    and
        //
        //      Theta1 <= Theta <= Theta2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles defining the arc,
        //    in radians.  Normally, THETA1 < THETA2.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double PN[2], a point on the circular arc which is
        //    nearest to the point.
        //
        //    Output, double *DIST, the distance to the nearest point.
        //
    {
        const int DIM_NUM = 2;

        int i;
        switch (r)
        {
            //
            //  Special case, the zero circle.
            //
            case 0.0:
            {
                typeMethods.r8vec_copy(DIM_NUM, pc, ref pn);

                dist = 0.0;
                for (i = 0; i < DIM_NUM; i++)
                {
                    dist += Math.Pow(p[i] - pn[i], 2);
                }

                dist = Math.Sqrt(dist);
                return;
            }
        }

        //
        //  Determine the angle made by the point.
        //
        double theta = typeMethods.r8_atan(p[1] - pc[1], p[0] - pc[0]);
        //
        //  If the angle is between THETA1 and THETA2, then you can
        //  simply project the point onto the arc.
        //
        if (typeMethods.r8_modp(theta - theta1, 2.0 * Math.PI) <=
            typeMethods.r8_modp(theta2 - theta1, 2.0 * Math.PI))
        {
            double r2 = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                r2 += Math.Pow(p[i] - pc[i], 2);
            }

            r2 = Math.Sqrt(r2);
            for (i = 0; i < DIM_NUM; i++)
            {
                pn[i] = pc[i] + (p[i] - pc[i]) * r / r2;
            }
        }
        //
        //  Otherwise, if the angle is less than the negative of the
        //  average of THETA1 and THETA2, it's on the side of the arc
        //  where the endpoint associated with THETA2 is closest.
        //
        else if (typeMethods.r8_modp(theta - 0.5 * (theta1 + theta2), 2.0 * Math.PI)
                 <= Math.PI)
        {
            pn[0] = pc[0] + r * Math.Cos(theta2);
            pn[1] = pc[1] + r * Math.Sin(theta2);
        }
        //
        //  Otherwise, the endpoint associated with THETA1 is closest.
        //
        else
        {
            pn[0] = pc[0] + r * Math.Cos(theta1);
            pn[1] = pc[1] + r * Math.Sin(theta1);
        }

        dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[i] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);
    }

    public static double circle_area_2d(double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_AREA_2D computes the area of a circle in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Output, double CIRCLE_AREA_2D, the area of the circle.
        //
    {
        double area = Math.PI * r * r;

        return area;
    }

    public static void circle_dia2imp_2d(double[] p1, double[] p2, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_DIA2IMP_2D converts a diameter to an implicit circle in 2D.
        //
        //  Discussion:
        //
        //    The diameter form of a circle is:
        //
        //      P1 and P2 are endpoints of a diameter.
        //
        //    The implicit form of a circle in 2D is:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], are the X and Y coordinates
        //    of two points which form a diameter of the circle.
        //
        //    Output, double *R, the computed radius of the circle.
        //
        //    Output, double PC[2], the computed center of the circle.
        //
    {
        r = 0.5 * Math.Sqrt(Math.Pow(p1[0] - p2[0], 2)
                            + Math.Pow(p1[1] - p2[1], 2));

        pc[0] = 0.5 * (p1[0] + p2[0]);
        pc[1] = 0.5 * (p1[1] + p2[1]);
    }

    public static int circle_exp_contains_point_2d(double[] p1, double[] p2, double[] p3,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_EXP_CONTAINS_POINT_2D determines if an explicit circle contains a point in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a circle in 2D is:
        //
        //      The circle passing through points P1, P2 and P3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], the coordinates of three
        //    points that lie on a circle.
        //
        //    Input, double P[2], the coordinates of a point, whose position
        //    relative to the circle is desired.
        //
        //    Output, int CIRCLE_EXP_CONTAINS_POINT_2D:
        //   -1, the three points are distinct and noncolinear,
        //      and the point lies inside the circle.
        //    0, the three points are distinct and noncolinear,
        //      and the point lies on the circle.
        //    1, the three points are distinct and noncolinear,
        //      and the point lies outside the circle.
        //    2, the three points are distinct and colinear,
        //      and the point lies on the line.
        //    3, the three points are distinct and colinear,
        //      and the point does not lie on the line.
        //    4, two points are distinct, and the point lies on the line.
        //    5, two points are distinct, and the point does not lie on the line.
        //    6, all three points are equal, and the point is equal to them,
        //    7, all three points are equal, and the point is not equal to them.
        //
    {
        double[] a = new double[4 * 4];
        double det;
        int inside;
        //
        //  P1 = P2?
        //
        if (typeMethods.r8vec_eq(2, p1, p2))
        {
            if (typeMethods.r8vec_eq(2, p1, p3))
            {
                inside = typeMethods.r8vec_eq(2, p1, p) ? 6 : 7;
            }
            else
            {
                det = (p1[0] - p3[0]) * (p[1] - p3[1])
                      - (p[0] - p3[0]) * (p1[1] - p3[1]);

                inside = det switch
                {
                    0.0 => 4,
                    _ => 5
                };
            }

            return inside;
        }

        //
        //  P1 does not equal P2.  Does P1 = P3?
        //
        if (typeMethods.r8vec_eq(2, p1, p3))
        {
            det = (p1[0] - p2[0]) * (p[1] - p2[1])
                  - (p[0] - p2[0]) * (p1[1] - p2[1]);

            inside = det switch
            {
                0.0 => 4,
                _ => 5
            };

            return inside;
        }

        //
        //  The points are distinct.  Are they colinear?
        //
        det = (p1[0] - p2[0]) * (p3[1] - p2[1])
              - (p3[0] - p2[0]) * (p1[1] - p2[1]);

        switch (det)
        {
            case 0.0:
            {
                det = (p1[0] - p2[0]) * (p[1] - p2[1])
                      - (p[0] - p2[0]) * (p1[1] - p2[1]);

                inside = det switch
                {
                    0.0 => 2,
                    _ => 3
                };

                return inside;
            }
        }

        //
        //  The points are distinct and non-colinear.
        //
        //  Compute the determinant
        //
        a[0 + 0 * 4] = p1[0];
        a[1 + 0 * 4] = p2[0];
        a[2 + 0 * 4] = p3[0];
        a[3 + 0 * 4] = p[0];

        a[0 + 1 * 4] = p1[1];
        a[1 + 1 * 4] = p2[1];
        a[2 + 1 * 4] = p3[1];
        a[3 + 1 * 4] = p[1];

        a[0 + 2 * 4] = p1[0] * p1[0] + p1[1] * p1[1];
        a[1 + 2 * 4] = p2[0] * p2[0] + p2[1] * p2[1];
        a[2 + 2 * 4] = p3[0] * p3[0] + p3[1] * p3[1];
        a[3 + 2 * 4] = p[0] * p[0] + p[1] * p[1];

        a[0 + 3 * 4] = 1.0;
        a[1 + 3 * 4] = 1.0;
        a[2 + 3 * 4] = 1.0;
        a[3 + 3 * 4] = 1.0;

        det = typeMethods.r8mat_det_4d(a);

        inside = det switch
        {
            < 0.0 => 1,
            0.0 => 0,
            _ => -1
        };

        return inside;
    }

    public static void circle_exp2imp_2d(double[] p1, double[] p2, double[] p3, ref double r,
            ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
        //
        //  Discussion:
        //
        //    The explicit form of a circle in 2D is:
        //
        //      The circle passing through points P1, P2 and P3.
        //
        //    Points P on an implicit circle in 2D satisfy the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    Any three points define a circle, as long as they don't lie on a straight
        //    line.  (If the points do lie on a straight line, we could stretch the
        //    definition of a circle to allow an infinite radius and a center at
        //    some infinite point.)
        //
        //    Instead of the formulas used here, you can use the linear system
        //    approach in the routine TRIANGLE_OUTCIRCLE_2D.
        //
        //    The diameter of the circle can be found by solving a 2 by 2 linear system.
        //    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
        //    and each forms a right triangle with the diameter.  Hence, the dot product
        //    of P2 - P1 with the diameter is equal to the square of the length
        //    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
        //    diameter vector originating at P1.
        //
        //    If all three points are equal, return a circle of radius 0 and
        //    the obvious center.
        //
        //    If two points are equal, return a circle of radius half the distance
        //    between the two distinct points, and center their average.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry,
        //    Second Edition,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], are the coordinates
        //    of three points that lie on the circle.  These points should be
        //    distinct, and not collinear.
        //
        //    Output, double *R, the radius of the circle.  Normally, R will be positive.
        //    R will be (meaningfully) zero if all three points are
        //    equal.  If two points are equal, R is returned as the distance between
        //    two nonequal points.  R is returned as -1 in the unlikely event that
        //    the points are numerically collinear; philosophically speaking, R
        //    should actually be "infinity" in this case.
        //
        //    Output, double PC[2], the center of the circle.
        //
    {
        const int DIM_NUM = 2;

        //
        //  If all three points are equal, then the
        //  circle of radius 0 and center P1 passes through the points.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2) && typeMethods.r8vec_eq(DIM_NUM, p1, p3))
        {
            r = 0.0;
            typeMethods.r8vec_copy(DIM_NUM, p1, ref pc);
            return;
        }

        //
        //  If exactly two points are equal, then the circle is defined as
        //  having the obvious radius and center.
        //
        if (typeMethods.r8vec_eq(DIM_NUM, p1, p2))
        {
            r = 0.5 * Math.Sqrt((p1[0] - p3[0]) * (p1[0] - p3[0])
                                + (p1[1] - p3[1]) * (p1[1] - p3[1]));
            pc[0] = 0.5 * (p1[0] + p3[0]);
            pc[1] = 0.5 * (p1[1] + p3[1]);
            return;
        }

        if (typeMethods.r8vec_eq(DIM_NUM, p1, p3))
        {
            r = 0.5 * Math.Sqrt((p1[0] - p2[0]) * (p1[0] - p2[0])
                                + (p1[1] - p2[1]) * (p1[1] - p2[1]));
            pc[0] = 0.5 * (p1[0] + p2[0]);
            pc[1] = 0.5 * (p1[1] + p2[1]);
            return;
        }
        if (typeMethods.r8vec_eq(DIM_NUM, p2, p3))
        {
            r = 0.5 * Math.Sqrt((p1[0] - p2[0]) * (p1[0] - p2[0])
                                + (p1[1] - p2[1]) * (p1[1] - p2[1]));
            pc[0] = 0.5 * (p1[0] + p2[0]);
            pc[1] = 0.5 * (p1[1] + p2[1]);
            return;
        }

        double a = p2[0] - p1[0];
        double b = p2[1] - p1[1];
        double c = p3[0] - p1[0];
        double d = p3[1] - p1[1];

        double e = a * (p1[0] + p2[0]) + b * (p1[1] + p2[1]);
        double f = c * (p1[0] + p3[0]) + d * (p1[1] + p3[1]);
        //
        //  Our formula is:
        //
        //    G = a * ( d - b ) - b * ( c - a )
        //
        //  but we get slightly better results using the original data.
        //
        double g = a * (p3[1] - p2[1]) - b * (p3[0] - p2[0]);
        switch (g)
        {
            //
            //  We check for collinearity.  A more useful check would compare the
            //  absolute value of G to a small quantity.
            //
            case 0.0:
                pc[0] = 0.0;
                pc[1] = 0.0;
                r = -1.0;
                return;
        }

        //
        //  The center is halfway along the diameter vector from P1.
        //
        pc[0] = 0.5 * (d * e - b * f) / g;
        pc[1] = 0.5 * (a * f - c * e) / g;
        //
        //  Knowing the center, the radius is now easy to compute.
        //
        r = Math.Sqrt((p1[0] - pc[0]) * (p1[0] - pc[0])
                      + (p1[1] - pc[1]) * (p1[1] - pc[1]));

    }

    public static bool circle_imp_contains_point_2d(double r, double[] pc, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_CONTAINS_POINT_2D determines if an implicit circle contains a point in 2D.
        //
        //  Discussion:
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool CIRCLE_IMP_CONTAINS_POINT_2D, is true if the point
        //    is inside or on the circle, false otherwise.
        //
    {
        return Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2) <= r * r;
    }

    public static void circle_imp_line_par_int_2d(double r, double[] pc, double x0, double y0,
            double f, double g, ref int int_num, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_LINE_PAR_INT_2D: ( implicit circle, parametric line ) intersection in 2D.
        //
        //  Discussion:
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    The parametric form of a line in 2D is:
        //
        //      X = X0 + F * T
        //      Y = Y0 + G * T
        //
        //    For normalization, we choose F*F+G*G = 1 and 0 <= F.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double F, G, X0, Y0, the parametric parameters of the line.
        //
        //    Output, int *INT_NUM, the number of intersecting points found.
        //    INT_NUM will be 0, 1 or 2.
        //
        //    Output, double P[2*2], contains the X and Y coordinates of
        //    the intersecting points.
        //
    {
        double t;

        double root = r * r * (f * f + g * g)
                      - (f * (pc[1] - y0) - g * (pc[0] - x0))
                      * (f * (pc[1] - y0) - g * (pc[0] - x0));

        switch (root)
        {
            case < 0.0:
                int_num = 0;
                break;
            case 0.0:
                int_num = 1;

                t = (f * (pc[0] - x0) + g * (pc[1] - y0)) / (f * f + g * g);
                p[0 + 0 * 2] = x0 + f * t;
                p[1 + 0 * 2] = y0 + g * t;
                break;
            case > 0.0:
                int_num = 2;

                t = (f * (pc[0] - x0) + g * (pc[1] - y0) - Math.Sqrt(root))
                    / (f * f + g * g);

                p[0 + 0 * 2] = x0 + f * t;
                p[1 + 0 * 2] = y0 + g * t;

                t = (f * (pc[0] - x0) + g * (pc[1] - y0) + Math.Sqrt(root))
                    / (f * f + g * g);

                p[0 + 1 * 2] = x0 + f * t;
                p[1 + 1 * 2] = y0 + g * t;
                break;
        }
    }

    public static double circle_imp_point_dist_2d(double r, double[] pc, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINT_DIST_2D: distance ( implicit circle, point ) in 2D.
        //
        //  Discussion:
        //
        //    The distance is zero if the point is on the circle.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double CIRCLE_IMP_POINT_DIST_2D, the distance of the point
        //    to the circle.
        //
    {
        double value = 0;

        value = Math.Sqrt(Math.Abs(Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2)
                                   - r * r));

        return value;
    }

    public static double circle_imp_point_dist_signed_2d(double r, double[] pc, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINT_DIST_SIGNED_2D: signed distance ( implicit circle, point ) in 2D.
        //
        //  Discussion:
        //
        //    The signed distance is zero if the point is on the circle.
        //    The signed distance is negative if the point is inside the circle.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double CIRCLE_IMP_POINT_DIST_SIGNED_2D, the signed distance
        //    of the point to the circle.  If the point is inside the circle,
        //    the signed distance is negative.
        //
    {
        double t = Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2) - r * r;

        double value = typeMethods.r8_sign(t) * Math.Sqrt(Math.Abs(t));

        return value;
    }

    public static double circle_imp_point_near_2d(double r, double[] pc, double[] p,
            ref double[] pn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINT_NEAR_2D: nearest ( implicit circle, point ) in 2D.
        //
        //  Discussion:
        //
        //    This routine finds the distance from a point to an implicitly
        //    defined circle, and returns the point on the circle that is
        //    nearest to the given point.
        //
        //    If the given point is the center of the circle, than any point
        //    on the circle is "the" nearest.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double PN[2], the nearest point on the circle.
        //
        //    Output, double CIRCLE_IMP_POINT_NEAR_2D, the distance of the point to the circle.
        //
    {
        const int DIM_NUM = 2;

        double dist;
        int i;

        if (typeMethods.r8vec_eq(DIM_NUM, p, pc))
        {
            dist = r;
            for (i = 0; i < DIM_NUM; i++)
            {
                pn[i] = pc[i] + r / Math.Sqrt(DIM_NUM);
            }

            return dist;
        }

        double r2 = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            r2 += Math.Pow(p[i] - pc[i], 2);
        }

        r2 = Math.Sqrt(r2);

        dist = Math.Abs(r2 - r);

        for (i = 0; i < DIM_NUM; i++)
        {
            pn[i] = pc[i] + r * (p[i] - pc[i]) / r2;
        }

        return dist;
    }

    public static double[] circle_imp_points_2d(double r, double[] pc, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINTS_2D returns N equally spaced points on an implicit circle in 2D.
        //
        //  Discussion:
        //
        //    The first point is always ( PC[0] + R, PC[1] ), and subsequent points
        //    proceed counter clockwise around the circle.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, int N, the number of points desired.  N must be at least 1.
        //
        //    Output, double CIRCLE_IMP_POINTS_2D[2*N], points on the circle.
        //
    {
        int j;

        double[] p = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            double angle = 2.0 * Math.PI * j / n;
            p[0 + j * 2] = pc[0] + r * Math.Cos(angle);
            p[1 + j * 2] = pc[1] + r * Math.Sin(angle);
        }

        return p;
    }

    public static double[] circle_imp_points_3d(double r, double[] pc, double[] nc, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINTS_3D returns points on an implicit circle in 3D.
        //
        //  Discussion:
        //
        //    Points P on an implicit circle in 3D satisfy the equations:
        //
        //      ( P(1) - PC(1) )^2
        //    + ( P(2) - PC(2) )^2
        //    + ( P(3) - PC(3) )^2 = R^2
        //
        //    and
        //
        //      ( P(1) - PC(1) ) * NC(1)
        //    + ( P(2) - PC(2) ) * NC(2)
        //    + ( P(3) - PC(3) ) * NC(3) = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[3], the center of the circle.
        //
        //    Input, double NC[3], a nonzero vector that is normal to
        //    the plane of the circle.  It is customary, but not necessary,
        //    that this vector have unit norm.
        //
        //    Input, int N, the number of points desired.
        //    N must be at least 1.
        //
        //    Output, double CIRCLE_IMP_POINTS_3D[3*N], the coordinates of points
        //    on the circle.
        //
    {
        const int DIM_NUM = 3;

        int j;
        //
        //  Get two unit vectors N1 and N2 which are orthogonal to each other,
        //  and to NC.
        //
        double[] n1 = new double[DIM_NUM];
        double[] n2 = new double[DIM_NUM];

        Plane.Geometry.plane_normal_basis_3d(pc, nc, ref n1, ref n2);
        //
        //  Rotate R units away from PC in the plane of N1 and N2.
        //
        double[] p = new double[DIM_NUM * n];

        for (j = 0; j < n; j++)
        {
            double theta = 2.0 * Math.PI * j / n;

            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                p[i + j * DIM_NUM] = pc[i] + r * (Math.Cos(theta) * n1[i]
                                                  + Math.Sin(theta) * n2[i]);
            }
        }

        return p;
    }

    public static void circle_imp_points_arc_2d(double r, double[] pc, double theta1,
            double theta2, int n, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINTS_ARC_2D returns N points on an arc of an implicit circle in 2D.
        //
        //  Discussion:
        //
        //    The first point is ( PC[0] + R * Math.Cos ( THETA1 ), PC[1] + R * Math.Sin ( THETA1 ) );
        //    The last point is  ( PC[0] + R * Math.Cos ( THETA2 ), PC[1] + R * Math.Sin ( THETA2 ) );
        //    and the intermediate points are evenly spaced in angle between these,
        //    and in counter clockwise order.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angular coordinates of the first
        //    and last points to be drawn, in radians.
        //
        //    Input, int N, the number of points desired.  N must be at least 1.
        //
        //    Output, double P[2*N], the coordinates of points on the circle.
        //
    {
        int i;
        //
        //  THETA3 is the smallest angle, no less than THETA1, which
        //  coincides with THETA2.
        //
        double theta3 = theta1 + typeMethods.r8_modp(theta2 - theta1, 2.0 * Math.PI);

        for (i = 0; i < n; i++)
        {
            double theta = n switch
            {
                > 1 => ((n - i - 1) * theta1 + i * theta3) / (n - 1),
                _ => 0.5 * (theta1 + theta3)
            };

            p[0 + i * 2] = pc[0] + r * Math.Cos(theta);
            p[1 + i * 2] = pc[1] + r * Math.Sin(theta);
        }

    }

    public static void circle_imp_print_2d(double r, double[] pc, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_PRINT_2D prints an implicit circle in 2D.
        //
        //  Discussion:
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, string TITLE, a title.
        //
    {
        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");
        Console.WriteLine("  Center = (" + pc[0] + ",  " + pc[1] + ")");

    }

    public static void circle_imp_print_3d(double r, double[] pc, double[] nc, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_PRINT_2D prints an implicit circle in 3D.
        //
        //  Discussion:
        //
        //    Points P on an implicit circle in 3D satisfy the equations:
        //
        //      ( P(1) - PC(1) )^2
        //    + ( P(2) - PC(2) )^2
        //    + ( P(3) - PC(3) )^2 = R^2
        //
        //    and
        //
        //      ( P(1) - PC(1) ) * NC(1)
        //    + ( P(2) - PC(2) ) * NC(2)
        //    + ( P(3) - PC(3) ) * NC(3) = 0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[3], the center of the circle.
        //
        //    Input, double NC[3], the normal vector to the circle.
        //
        //    Input, string TITLE, a title.
        //
    {
        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        Console.WriteLine("  Radius = " + r + "");
        Console.WriteLine("  Center = (" + pc[0]
                                         + ", " + pc[1]
                                         + ", " + pc[2] + ")");
        Console.WriteLine("  Normal = (" + nc[0]
                                         + ", " + nc[1]
                                         + ", " + nc[2] + ")");

    }

    public static void circle_imp2exp_2d(double r, double[] pc, double[] p1, double[] p2,
            double[] p3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP2EXP_2D converts a circle from implicit to explicit form in 2D.
        //
        //  Discussion:
        //
        //    Points P on an implicit circle in 2D satisfy the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    The explicit form of a circle in 2D is:
        //
        //      The circle passing through points P1, P2 and P3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry,
        //    Second Edition,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, double R, PC[2], the radius and center of the circle.
        //
        //    Output, double P1[2], P2[2], P3[2], three points on the circle.
        //
    {
        double theta = 0.0;
        p1[0] = pc[0] + r * Math.Cos(theta);
        p1[1] = pc[1] + r * Math.Sin(theta);

        theta = 2.0 * Math.PI / 3.0;
        p2[0] = pc[0] + r * Math.Cos(theta);
        p2[1] = pc[1] + r * Math.Sin(theta);

        theta = 4.0 * Math.PI / 3.0;
        p3[0] = pc[0] + r * Math.Cos(theta);
        p3[1] = pc[1] + r * Math.Sin(theta);
    }

    public static double[] circle_llr2imp_2d(double[] p1, double[] p2, double[] q1, double[] q2,
            double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LLR2IMP_2D converts a circle from LLR to implicit form in 2D.
        //
        //  Discussion:
        //
        //    The LLR form of a circle in 2D is:
        //
        //      The circle of radius R tangent to the lines L1 and L2.
        //
        //    The implicit form of a circle in 2D is:
        //
        //      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
        //
        //    Let S be the scaled distance of a point on L1 from P1 to P2,
        //    and let N1 be a unit normal vector to L1.  Then a point P that is
        //    R units from L1 satisfies:
        //
        //      P = P1 + s * ( P2 - P1 ) + R * N1.
        //
        //    Let t be the scaled distance of a point on L2 from Q1 to Q2,
        //    and let N2 be a unit normal vector to L2.  Then a point Q that is
        //    R units from L2 satisfies:
        //
        //      Q = Q1 + t * ( Q2 - Q1 ) + R * N2.
        //
        //    For the center of the circle, then, we have P = Q, that is
        //
        //      ( P2 - P1 ) * s + ( Q2 - Q1 ) * t = - P1 - Q1 - R * ( N1 + N2 )
        //
        //    This is a linear system for ( s and t ) from which we can compute
        //    the points of tangency, and the center.
        //
        //    Note that we have four choices for the circle based on the use
        //    of plus or minus N1 and plus or minus N2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two points on line 1.
        //
        //    Input, double Q1[2], Q2[2], two points on line 2.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Output, double CIRCLE_LLR2IMP_2D[2*4], the centers of the circles.
        //
    {
        double det = 0;

        double[] a = new double[2 * 2];
        double[] b = new double[2];
        double[] pc = new double[2 * 4];
        //
        //  Compute the normals N1 and N2.
        //
        double[] n1 = LineNS.Geometry.line_exp_normal_2d(p1, p2);

        double[] n2 = LineNS.Geometry.line_exp_normal_2d(q1, q2);
        //
        //  Set the linear system.
        //
        a[0 + 0 * 2] = p2[0] - p1[0];
        a[1 + 0 * 2] = p2[1] - p1[1];
        a[0 + 1 * 2] = -q2[0] + q1[0];
        a[1 + 1 * 2] = -q2[1] + q1[1];
        //
        //  Solve the 4 linear systems, using every combination of
        //  signs on the normal vectors.
        //
        b[0] = -p1[0] + q1[0] + r * n1[0] + r * n2[0];
        b[1] = -p1[1] + q1[1] + r * n1[1] + r * n2[1];

        double[] x = typeMethods.r8mat_solve_2d(a, b, ref det);

        pc[0 + 2 * 0] = p1[0] + (p2[0] - p1[0]) * x[0] - r * n1[0];
        pc[1 + 2 * 0] = p1[1] + (p2[1] - p1[1]) * x[0] - r * n1[1];

        b[0] = -p1[0] + q1[0] + r * n1[0] - r * n2[0];
        b[1] = -p1[1] + q1[1] + r * n1[1] - r * n2[1];

        x = typeMethods.r8mat_solve_2d(a, b, ref det);

        pc[0 + 2 * 1] = p1[0] + (p2[0] - p1[0]) * x[0] - r * n1[0];
        pc[1 + 2 * 1] = p1[1] + (p2[1] - p1[1]) * x[0] - r * n1[1];

        b[0] = -p1[0] + q1[0] - r * n1[0] + r * n2[0];
        b[1] = -p1[1] + q1[1] - r * n1[1] + r * n2[1];

        x = typeMethods.r8mat_solve_2d(a, b, ref det);

        pc[0 + 2 * 2] = p1[0] + (p2[0] - p1[0]) * x[0] + r * n1[0];
        pc[1 + 2 * 2] = p1[1] + (p2[1] - p1[1]) * x[0] + r * n1[1];

        b[0] = -p1[0] + q1[0] - r * n1[0] - r * n2[0];
        b[1] = -p1[1] + q1[1] - r * n1[1] - r * n2[1];

        x = typeMethods.r8mat_solve_2d(a, b, ref det);

        pc[0 + 2 * 3] = p1[0] + (p2[0] - p1[0]) * x[0] + r * n1[0];
        pc[1 + 2 * 3] = p1[1] + (p2[1] - p1[1]) * x[0] + r * n1[1];

        return pc;
    }

    public static double circle_lune_angle_by_height_2d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D computes the angle of a circular lune.
        //
        //  Discussion:
        //
        //    Draw the chord connecting two points on the circumference of a circle.
        //    The region between the chord and the circumference is a "lune".
        //    We wish to know the angle subtended by the lune.
        //
        //    The distance from the center of the circle to the midpoint of the chord
        //    is the "height" H of the lune.  It is natural to expect 0 <= H <= R.
        //    However, if we allow -R <= H < 0 as well, this allows us to include
        //    lunes which involve more than half the circle's area.
        //
        //    If H < -R or R < H, then no lune is formed, and we return a zero angle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double H, the height of the lune.
        //
        //    Output, double ANGLE, the angle of the lune.
        //
    {
        double angle;

        if (-r <= h && h <= r)
        {
            angle = 2.0 * Math.Acos(h / r);
        }
        else
        {
            angle = 0.0;
        }

        return angle;
    }

    public static double circle_lune_area_by_angle_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_AREA_BY_ANGLE_2D returns the area of a circular lune in 2D.
        //
        //  Discussion:
        //
        //    A lune is formed by drawing a circular arc, and joining its endpoints.
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
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles defining the arc,
        //    in radians.  Normally, THETA1 < THETA2.
        //
        //    Output, double CIRCLE_LUNE_AREA_BY_ANGLE_2D, the area of the lune.
        //
    {
        double area_sector = circle_sector_area_2d(r, pc, theta1, theta2);
        double area_triangle = circle_triangle_area_2d(r, pc, theta1, theta2);
        double area = area_sector - area_triangle;

        return area;
    }

    public static double circle_lune_area_by_height_2d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_AREA_BY_HEIGHT_2D returns the area of a circular lune in 2D.
        //
        //  Discussion:
        //
        //    A lune is formed by drawing a circular arc, and joining its endpoints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double H, the height of the lune.
        //
        //    Output, double CIRCLE_LUNE_AREA_BY_HEIGHT_2D, the area of the lune.
        //
    {
        double area;

        if (-r <= h && h <= r)
        {
            area = r * r * Math.Acos(h / r) - h * Math.Sqrt(r * r - h * h);
        }
        else
        {
            area = 0.0;
        }

        return area;
    }

    public static double[] circle_lune_centroid_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_CENTROID_2D returns the centroid of a circular lune in 2D.
        //
        //  Discussion:
        //
        //    A lune is formed by drawing a circular arc, and joining its endpoints.
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles of the first and last points
        //    on the circular arc.
        //
        //    Output, double CIRCLE_LUNE_CENTROID_2D[2], the coordinates of the centroid
        //    of the lune.
        //
    {
        double theta = theta2 - theta1;

        double d = theta switch
        {
            0.0 => r,
            _ => 4.0 * r * Math.Pow(Math.Sin(0.5 * theta), 3) / (3.0 * (theta - Math.Sin(theta)))
        };

        double[] centroid = new double[2];

        centroid[0] = pc[0] + d * Math.Cos(theta);
        centroid[1] = pc[1] + d * Math.Sin(theta);

        return centroid;
    }

    public static double circle_lune_height_by_angle_2d(double r, double angle)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D computes the height of a circular lune.
        //
        //  Discussion:
        //
        //    Draw the chord connecting two points on the circumference of a circle.
        //    The region between the chord and the circumference is a "lune".
        //    The lune subtends a given angle between 0 and 2 pi.
        //
        //    The distance from the center of the circle to the midpoint of the chord
        //    is the "height" H of the lune and we wish to determine this value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double ANGLE, the angle subtended by the lune.
        //
        //    Output, double HEIGHT, the height of the lune
        //
    {
        double height = r * Math.Cos(angle / 2.0);

        return height;
    }

    public static void circle_pppr2imp_3d(double[] p1, double[] p2, double[] p3, double r,
            ref double[] pc, ref double[] normal)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_PPR2IMP_3D converts a circle from PPR to implicit form in 3D.
        //
        //  Discussion:
        //
        //    The PPPR form of a circle in 3D is:
        //
        //      The circle of radius R passing through points P1 and P2,
        //      and lying in the plane of P1, P2 and P3.
        //
        //    The implicit form of a circle in 3D is:
        //
        //      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) )^2 = R^2
        //      and the dot product of P - PC with NORMAL is 0.
        //
        //    There may be zero, one, or two circles that satisfy the
        //    requirements of the PPPR form.
        //
        //    If there is no such circle, then PC(1:2,1) and PC(1:2,2)
        //    are set to the midpoint of (P1,P2).
        //
        //    If there is one circle, PC(1:2,1) and PC(1:2,2) will be equal.
        //
        //    If there are two circles, then PC(1:2,1) is the first center,
        //    and PC(1:2,2) is the second.
        //
        //    This calculation is equivalent to finding the intersections of
        //    spheres of radius R at points P1 and P2, which lie in the plane
        //    defined by P1, P2 and P3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], two points on the circle.
        //
        //    Input, double P3[3], a third point.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Output, double PC[3*2], the centers of the two circles.
        //
        //    Output, double NORMAL[3], the normal to the circles.
        //
    {
        const int DIM_NUM = 3;

        int i;
        //
        //  Compute the distance from P1 to P2.
        //
        double dist = typeMethods.r8vec_distance(DIM_NUM, p1, p2);
        //
        //  If R is smaller than DIST, we don't have a circle.
        //
        if (2.0 * r < dist)
        {
            int j;
            for (j = 0; j < 2; j++)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    pc[i + j * DIM_NUM] = 0.5 * (p1[i] + p2[i]);
                }
            }
        }

        //
        //  H is the distance from the midpoint of (P1,P2) to the center.
        //
        double h = Math.Sqrt((r + 0.5 * dist) * (r - 0.5 * dist));
        //
        //  Define a unit direction V that is normal to P2-P1, and lying
        //  in the plane (P1,P2,P3).
        //
        //  To do this, subtract from P3-P1 the component in the direction P2-P1.
        //
        double[] v = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            v[i] = p3[i] - p1[i];
        }

        double dot = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dot += v[i] * (p2[i] - p1[i]);
        }

        dot /= dist;
        for (i = 0; i < DIM_NUM; i++)
        {
            v[i] -= dot * (p2[i] - p1[i]) / dist;
        }

        double length = typeMethods.r8vec_norm(DIM_NUM, v);
        for (i = 0; i < DIM_NUM; i++)
        {
            v[i] /= length;
        }

        //
        //  We can go with or against the given normal direction.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            pc[i + 0 * DIM_NUM] = 0.5 * (p2[i] + p1[i]) + h * v[i];
            pc[i + 1 * DIM_NUM] = 0.5 * (p2[i] + p1[i]) - h * v[i];
        }

        Plane.Geometry.plane_exp_normal_3d(p1, p2, p3, ref normal);
    }

    public static double[] circle_ppr2imp_2d(double[] p1, double[] p2, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_PPR2IMP_2D converts a circle from PPR to implicit form in 2D.
        //
        //  Discussion:
        //
        //    The PPR form of a circle in 2D is:
        //
        //      The circle of radius R passing through points P1 and P2.
        //
        //    The implicit form of a circle in 2D is:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = R * R
        //
        //    There may be zero, one, or two circles that satisfy the
        //    requirements of the PPR form.
        //
        //    If there is no such circle, then the two "solutions" are set to
        //    the midpoint of P1 and P2.
        //
        //    If there is one circle, then the two solutions will be set equal
        //    to the midpoint of P1 and P2.
        //
        //    If there are two distinct circles, then (PC[0],PC[1]) is the first center,
        //    and (PC[2],PC[3]) is the second.
        //
        //    This calculation is equivalent to finding the intersections of
        //    circles of radius R at points P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], two points on the circle.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Output, double PC[2*2], the centers of the two circles.
        //
    {
        const int DIM_NUM = 2;

        double[] pc = new double[DIM_NUM * 2];
        //
        //  Compute the distance from P1 to P2.
        //
        double dist = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2) + Math.Pow(p2[1] - p1[1], 2));
        //
        //  If R is smaller than DIST, we don't have a circle.
        //
        if (2.0 * r < dist)
        {
            int j;
            for (j = 0; j < 2; j++)
            {
                int i;
                for (i = 0; i < DIM_NUM; i++)
                {
                    pc[i + j * DIM_NUM] = 0.5 * (p1[i] + p2[i]);
                }
            }
        }

        //
        //  H is the distance from the midpoint of (P1,P2) to the center.
        //
        double h = Math.Sqrt((r + 0.5 * dist) * (r - 0.5 * dist));
        //
        //  The center is found by going midway between P1 and P2, and then
        //  H units in the unit perpendicular direction.
        //
        //  We actually have two choices for the normal direction.
        //
        pc[0 + 0 * DIM_NUM] = 0.5 * (p2[0] + p1[0]) + h * (p2[1] - p1[1]) / dist;
        pc[1 + 0 * DIM_NUM] = 0.5 * (p2[1] + p1[1]) - h * (p2[0] - p1[0]) / dist;

        pc[0 + 1 * DIM_NUM] = 0.5 * (p2[0] + p1[0]) - h * (p2[1] - p1[1]) / dist;
        pc[1 + 1 * DIM_NUM] = 0.5 * (p2[1] + p1[1]) + h * (p2[0] - p1[0]) / dist;

        return pc;
    }

    public static double circle_sector_area_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_AREA_2D computes the area of a circular sector in 2D.
        //
        //  Discussion:
        //
        //    A circular sector is formed by a circular arc, and the two straight line
        //    segments that join its ends to the center of the circle.
        //
        //    A circular sector is defined by
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    and
        //
        //      Theta1 <= Theta <= Theta2
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
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles of the first and last points
        //    on the circular arc.
        //
        //    Output, double CIRCLE_SECTOR_AREA_2D, the area of the circle.
        //
    {
        double area = 0.5 * r * r * (theta2 - theta1);

        return area;
    }

    public static double[] circle_sector_centroid_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_CENTROID_2D returns the centroid of a circular sector in 2D.
        //
        //  Discussion:
        //
        //    A circular sector is formed by a circular arc, and the two straight line
        //    segments that join its ends to the center of the circle.
        //
        //    A circular sector is defined by
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    and
        //
        //      Theta1 <= Theta <= Theta2
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
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the coordinates of the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles of the first and last points
        //    on the circular arc.
        //
        //    Output, double CIRCLE_SECTOR_CENTROID_2D[2], the coordinates
        //    of the centroid of the sector.
        //
    {
        const int DIM_NUM = 2;

        double theta = theta2 - theta1;

        double d = theta switch
        {
            0.0 => 2.0 * r / 3.0,
            _ => 4.0 * r * Math.Sin(0.5 * theta) / (3.0 * theta)
        };

        double[] centroid = new double[DIM_NUM];

        centroid[0] = pc[0] + d * Math.Cos(theta);
        centroid[1] = pc[1] + d * Math.Sin(theta);

        return centroid;
    }

    public static bool circle_sector_contains_point_2d(double r, double[] pc, double theta1,
            double theta2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_CONTAINS_POINT_2D : is a point inside a circular sector?
        //
        //  Discussion:
        //
        //    A circular sector is formed by a circular arc, and the two straight line
        //    segments that join its ends to the center of the circle.
        //
        //    A circular sector is defined by
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    and
        //
        //      Theta1 <= Theta <= Theta2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles defining the arc,
        //    in radians.  Normally, THETA1 < THETA2.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, logical CIRCLE_SECTOR_CONTAINS_POINT_2D, is TRUE if the point is
        //    inside or on the circular sector.
        //
    {
        bool inside = false;
        //
        //  Is the point inside the (full) circle?
        //
        if (Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2) <= r * r)
        {
            //
            //  Is the point's angle within the arc's range?
            //  Try to force the angles to lie between 0 and 2 * PI.
            //
            double theta = typeMethods.r8_atan(p[1] - pc[1], p[0] - pc[0]);

            if (typeMethods.r8_modp(theta - theta1, 2.0 * Math.PI) <=
                typeMethods.r8_modp(theta2 - theta1, 2.0 * Math.PI))
            {
                inside = true;
            }
        }

        return inside;
    }

    public static void circle_sector_print_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_PRINT_2D prints a circular sector in 2D.
        //
        //  Discussion:
        //
        //    A circular sector is formed by a circular arc, and the two straight line
        //    segments that join its ends to the center of the circle.
        //
        //    A circular sector is defined by
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    and
        //
        //      Theta1 <= Theta <= Theta2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles defining the arc,
        //    in radians.  Normally, THETA1 < THETA2.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("  Circular sector definition:");
        Console.WriteLine("");
        Console.WriteLine("    Radius = " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        Console.WriteLine("    Center = " + pc[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                          + "  " + pc[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        Console.WriteLine("    Theta  = " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                          + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

    }

    public static double circle_triangle_area_2d(double r, double[] pc, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_TRIANGLE_AREA_2D returns the area of a circle triangle in 2D.
        //
        //  Discussion:
        //
        //    A circle triangle is formed by drawing a circular arc, and considering
        //    the triangle formed by the endpoints of the arc plus the center of
        //    the circle.
        //
        //    Note that for angles greater than PI, the triangle will actually
        //    have NEGATIVE area.
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
        //    Input, double R, the radius of the circle.
        //
        //    Input, double PC[2], the center of the circle.
        //
        //    Input, double THETA1, THETA2, the angles defining the arc,
        //    in radians.  Normally, THETA1 < THETA2.
        //
        //    Output, double AREA, the (signed) area of the triangle.
        //
    {
        double area = 0.5 * r * r * Math.Sin(theta2 - theta1);

        return area;
    }

    public static void circle_triple_angles_2d(double r1, double r2, double r3, ref double angle1,
            ref double angle2, ref double angle3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_TRIPLE_ANGLE_2D returns an angle formed by three circles in 2D.
        //
        //  Discussion:
        //
        //    A circle triple is a set of three tangent circles.  We assume
        //    that no circle is contained in another.
        //
        //    We consider the triangle formed by joining the centers of the circles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Kenneth Stephenson,
        //    Circle Packing, The Theory of Discrete Analytic Functions,
        //    Cambridge, 2005.
        //
        //  Parameters:
        //
        //    Input, double R1, R2, R3, the radii of the circles.
        //
        //    Input, double *ANGLE1, *ANGLE2, *ANGLE3, the angles
        //    in the triangle.
        //
    {
        angle1 = typeMethods.r8_acos(
                     Math.Pow(r1 + r2, 2) + Math.Pow(r1 + r3, 2) - Math.Pow(r2 + r3, 2)) /
                 (2.0 * (r1 + r2) * (r1 + r3));

        angle2 = typeMethods.r8_acos(
                     Math.Pow(r2 + r3, 2) + Math.Pow(r2 + r1, 2) - Math.Pow(r3 + r1, 2)) /
                 (2.0 * (r2 + r3) * (r2 + r1));

        angle3 = typeMethods.r8_acos(
                     Math.Pow(r3 + r1, 2) + Math.Pow(r3 + r2, 2) - Math.Pow(r1 + r2, 2)) /
                 (2.0 * (r3 + r1) * (r3 + r2));

    }

    public static double circles_intersect_area_2d(double r1, double r2, double d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLES_INTERSECT_AREA_2D: area of the intersection of two circles.
        //
        //  Discussion:
        //
        //    Circles of radius R1 and R2 are D units apart.  What is the area of
        //    intersection?
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the radiuses of the circles.
        //    R1 and R2 should be positive.
        //
        //    Input, double D, the distance between the circular centers.
        //    D must be positive, and should be no greater than R1 + R2.
        //
        //    Output, double AREA, the area of the intersection.
        //
    {
        double area;

        if (r1 + r2 < d)
        {
            area = 0.0;
        }
        else
        {
            switch (d)
            {
                case 0.0:
                    area = Math.PI * Math.Pow(Math.Min(r1, r2), 2);
                    break;
                default:
                    double h1 = 0.5 * (d * d + r1 * r1 - r2 * r2) / d;
                    double area1 = circle_lune_area_by_height_2d(r1, h1);
                    double h2 = 0.5 * (d * d - r1 * r1 + r2 * r2) / d;
                    double area2 = circle_lune_area_by_height_2d(r2, h2);
                    area = area1 + area2;
                    break;
            }
        }

        return area;
    }

    public static void circles_intersect_points_2d(double r1, double[] pc1, double r2,
            double[] pc2, ref int int_num, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLES_INTERSECT_POINTS_2D: intersection points of two circles in 2D.
        //
        //  Discussion:
        //
        //    Two circles can intersect in 0, 1, 2 or infinitely many points.
        //
        //    The 0 and 2 intersection cases are numerically robust; the 1 and
        //    infinite intersection cases are numerically fragile.  The routine
        //    uses a tolerance to try to detect the 1 and infinite cases.
        //
        //    An implicit circle in 2D satisfies the equation:
        //
        //      Math.Pow ( P[0] - PC[0], 2 ) + Math.Pow ( P[1] - PC[1], 2 ) = Math.Pow ( R, 2 )
        //
        //    Thanks to Mario Pintaric for pointing out, on 13 March 2006,
        //    a place where (R1-R2) had been mistakenly written as (R1-R1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, the radius of the first circle.
        //
        //    Input, double PC1[2], the coordinates of the center of the first circle.
        //
        //    Input, double R2, the radius of the second circle.
        //
        //    Input, double PC2[2], the coordinates of the center of the second circle.
        //
        //    Output, int *INT_NUM, the number of intersecting points found.
        //    INT_NUM will be 0, 1, 2 or 3.  3 indicates that there are an infinite
        //    number of intersection points.
        //
        //    Output, double P[2*2], if INT_NUM is 1 or 2, the
        //    coordinates of the intersecting points.
        //
    {
        double tol = typeMethods.r8_epsilon();

        p[0 + 0 * 2] = 0.0;
        p[1 + 0 * 2] = 0.0;
        p[0 + 1 * 2] = 0.0;
        p[1 + 1 * 2] = 0.0;
        //
        //  Take care of the case in which the circles have the same center.
        //
        double t1 = (Math.Abs(pc1[0] - pc2[0]) + Math.Abs(pc1[1] - pc2[1])) / 2.0;
        double t2 = (Math.Abs(pc1[0]) + Math.Abs(pc2[0])
                                      + Math.Abs(pc1[1]) + Math.Abs(pc2[1]) + 1.0) / 5.0;

        if (t1 <= tol * t2)
        {
            t1 = Math.Abs(r1 - r2);
            t2 = (Math.Abs(r1) + Math.Abs(r2) + 1.0) / 3.0;

            int_num = t1 <= tol * t2 ? 3 : 0;

            return;
        }

        double distsq = (pc1[0] - pc2[0]) * (pc1[0] - pc2[0])
                        + (pc1[1] - pc2[1]) * (pc1[1] - pc2[1]);

        double root = 2.0 * (r1 * r1 + r2 * r2) * distsq - distsq * distsq
                                                         - (r1 - r2) * (r1 - r2) * (r1 + r2) * (r1 + r2);

        if (root < -tol)
        {
            int_num = 0;
            return;
        }

        double sc1 = (distsq - (r2 * r2 - r1 * r1)) / distsq;

        if (root < tol)
        {
            int_num = 1;
            p[0 + 0 * 2] = pc1[0] + 0.5 * sc1 * (pc2[0] - pc1[0]);
            p[1 + 0 * 2] = pc1[1] + 0.5 * sc1 * (pc2[1] - pc1[1]);
            return;
        }

        double sc2 = Math.Sqrt(root) / distsq;

        int_num = 2;

        p[0 + 0 * 2] = pc1[0] + 0.5 * sc1 * (pc2[0] - pc1[0])
                       - 0.5 * sc2 * (pc2[1] - pc1[1]);
        p[1 + 0 * 2] = pc1[1] + 0.5 * sc1 * (pc2[1] - pc1[1])
                              + 0.5 * sc2 * (pc2[0] - pc1[0]);

        p[0 + 1 * 2] = pc1[0] + 0.5 * sc1 * (pc2[0] - pc1[0])
                              + 0.5 * sc2 * (pc2[1] - pc1[1]);
        p[1 + 1 * 2] = pc1[1] + 0.5 * sc1 * (pc2[1] - pc1[1])
                       - 0.5 * sc2 * (pc2[0] - pc1[0]);

    }
}