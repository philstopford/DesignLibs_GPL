using System;
using Burkardt.Types;

namespace Burkardt.Geometry
{
    public static class Angle
    {
        public static void angle_box_2d(double dist, double[] p1, double[] p2, double[] p3,
                ref double[] p4, ref double[] p5)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_BOX_2D "boxes" an angle defined by three points in 2D.
            //
            //  Discussion:
            //
            //    The routine is given points P1, P2 and P3, determining the two lines:
            //      P1 to P2
            //    and
            //      P2 to P3
            //    and a nonnegative distance
            //      DIST.
            //
            //    The routine returns a pair of "corner" points
            //      P4 and P5
            //    both of which are a distance DIST from both lines, and in fact,
            //    both of which are a distance DIST from P2.
            //
            //                         /  P3
            //                        /   /   /
            //     - - - - - - - - -P4 - / -P6 - - -
            //                      /   /   /
            //    P1---------------/--P2-----------------
            //                    /   /   /
            //     - - - - - - -P7 - / -P5 - - - - -
            //                  /   /   /
            //
            //    In the illustration, P1, P2 and P3 represent
            //    the points defining the lines.
            //
            //    P4 and P5 represent the desired "corner points", which
            //    are on the positive or negative sides of both lines.
            //
            //    The numbers P6 and P7 represent the undesired points, which
            //    are on the positive side of one line and the negative of the other.
            //
            //    Special cases:
            //
            //    if P1 = P2, this is the same as extending the line from
            //    P3 through P2 without a bend.
            //
            //    if P3 = P2, this is the same as extending the line from
            //    P1 through P2 without a bend.
            //
            //    if P1 = P2 = P3 this is an error.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double DIST, the nonnegative distance from P1
            //    to the computed points P4 and P5.
            //
            //    Input, double P1[2], P2[2], P3[2].
            //    P1 and P2 are distinct points that define a line.
            //    P2 and P3 are distinct points that define a line.
            //
            //    Output, double P4[2], P5[2], points which lie DIST units from
            //    the line between P1 and P2, and from the line between P2 and P3.
            //
        {
            int DIM_NUM = 2;

            double stheta;
            double temp;
            double temp1;
            double temp2;
            double[] u = new double[DIM_NUM];
            double[] u1 = new double[DIM_NUM];
            double[] u2 = new double[DIM_NUM];
            //
            //  If DIST = 0, assume the user knows best.
            //
            if (dist == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, p2, ref p4);
                typeMethods.r8vec_copy(DIM_NUM, p2, ref p5);
                return;
            }

            //
            //  Fail if all three points are equal.
            //
            if (typeMethods.r8vec_eq(DIM_NUM, p1, p2) && typeMethods.r8vec_eq(DIM_NUM, p2, p3))
            {
                Console.WriteLine("");
                Console.WriteLine("ANGLE_BOX_2D - Fatal error!");
                Console.WriteLine("  Input points P3 = P2 = P3.");
                typeMethods.r8vec_print(DIM_NUM, p1, "  P1:");
                return;
            }

            //
            //  If P1 = P2, extend the line through the doubled point.
            //
            if (typeMethods.r8vec_eq(DIM_NUM, p1, p2))
            {
                u2[0] = p3[1] - p2[1];
                u2[1] = p2[0] - p3[0];
                temp = typeMethods.r8vec_norm(DIM_NUM, u2);
                u2[0] = u2[0] / temp;
                u2[1] = u2[1] / temp;
                p4[0] = p2[0] + dist * u2[0];
                p4[1] = p2[1] + dist * u2[1];
                p5[0] = p2[0] - dist * u2[0];
                p5[1] = p2[1] - dist * u2[1];
                return;
            }

            //
            //  If P2 = P3, extend the line through the doubled point.
            //
            if (typeMethods.r8vec_eq(DIM_NUM, p2, p3))
            {
                u1[0] = p1[1] - p2[1];
                u1[1] = p2[0] - p1[0];
                temp = typeMethods.r8vec_norm(DIM_NUM, u1);
                u1[0] = u1[0] / temp;
                u1[1] = u1[1] / temp;
                p4[0] = p2[0] + dist * u1[0];
                p4[1] = p2[1] + dist * u1[1];
                p5[0] = p2[0] - dist * u1[0];
                p5[1] = p2[1] - dist * u1[1];
                return;
            }

            //
            //  Now compute the unit normal vectors to each line.
            //  We choose the sign so that the unit normal to line 1 has
            //  a positive dot product with line 2.
            //
            u1[0] = p1[1] - p2[1];
            u1[1] = p2[0] - p1[0];
            temp = typeMethods.r8vec_norm(DIM_NUM, u1);
            u1[0] = u1[0] / temp;
            u1[1] = u1[1] / temp;

            temp1 = u1[0] * (p3[0] - p2[0])
                    + u1[1] * (p3[1] - p2[1]);

            if (temp1 < 0.0)
            {
                u1[0] = -u1[0];
                u1[1] = -u1[1];
            }

            u2[0] = p3[1] - p2[1];
            u2[1] = p2[0] - p3[0];
            temp = typeMethods.r8vec_norm(DIM_NUM, u2);
            u2[0] = u2[0] / temp;
            u2[1] = u2[1] / temp;

            temp2 = u2[0] * (p1[0] - p2[0])
                    + u2[1] * (p1[1] - p2[1]);

            if (temp2 < 0.0)
            {
                u2[0] = -u2[0];
                u2[1] = -u2[1];
            }

            //
            //  Try to catch the case where we can't determine the
            //  sign of U1, because both U1 and -U1 are perpendicular
            //  to (P3-P2), and similarly for U2 and (P1-P2).
            //
            if (temp1 == 0.0 || temp2 == 0.0)
            {
                if (u1[0] * u2[0] + u1[1] * u2[1] < 0.0)
                {
                    u1[0] = -u1[0];
                    u2[0] = -u2[0];
                }
            }

            //
            //  Try to catch a line turning back on itself, evidenced by
            //    Cos(theta) = (P3-P2) dot (P2-P1) / ( norm(P3-P2) * norm(P2-P1) )
            //  being -1, or very close to -1.
            //
            temp = (p3[0] - p2[0]) * (p2[0] - p1[0])
                   + (p3[1] - p2[1]) * (p2[1] - p1[1]);

            temp1 = Math.Sqrt(Math.Pow(p3[0] - p2[0], 2) + Math.Pow(p3[1] - p2[1], 2));
            temp2 = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2) + Math.Pow(p2[1] - p1[1], 2));

            temp = temp / (temp1 * temp2);

            if (temp < -0.99)
            {
                temp = Math.Sqrt(Math.Pow(p2[0] - p1[0], 2) + Math.Pow(p2[1] - p1[1], 2));

                p4[0] = p2[0] + dist * (p2[0] - p1[0]) / temp + dist * u1[0];
                p4[1] = p2[1] + dist * (p2[1] - p1[1]) / temp + dist * u1[1];
                p5[0] = p2[0] + dist * (p2[0] - p1[0]) / temp - dist * u1[0];
                p5[1] = p2[1] + dist * (p2[1] - p1[1]) / temp - dist * u1[1];
                return;
            }

            //
            //  Compute the "average" unit normal vector.
            //
            //  The average of the unit normals could be zero, but only when
            //  the second line has the same direction and opposite sense
            //  of the first, and we've already checked for that case.
            //
            //  Well, check again!  This problem "bit" me in the case where
            //  P1 = P2, which I now treat specially just to guarantee I
            //  avoid this problem!
            //
            if (u1[0] * u2[0] + u1[1] * u2[1] < 0.0)
            {
                u2[0] = -u2[0];
                u2[1] = -u2[1];
            }

            u[0] = 0.5 * (u1[0] + u2[0]);
            u[1] = 0.5 * (u1[1] + u2[1]);

            temp = typeMethods.r8vec_norm(DIM_NUM, u);
            u[0] = u[0] / temp;
            u[1] = u[1] / temp;
            //
            //  You must go DIST/STHETA units along this unit normal to
            //  result in a distance DIST from line1 (and line2).
            //
            stheta = u[0] * u1[0] + u[1] * u1[1];

            p4[0] = p2[0] + dist * u[0] / stheta;
            p4[1] = p2[1] + dist * u[1] / stheta;
            p5[0] = p2[0] - dist * u[0] / stheta;
            p5[1] = p2[1] - dist * u[1] / stheta;
        }

        public static bool angle_contains_ray_2d(double[] p1, double[] p2, double[] p3,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_CONTAINS_RAY_2D determines if an angle contains a ray, in 2D.
            //
            //  Discussion:
            //
            //    The angle is defined by the sequence of points P1, P2, P3.
            //
            //    The ray is defined by the sequence of points P2, P.
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
            //    Input, double P1[2], P2[2], P3[2], the coordinates of points on
            //    the angle.
            //
            //    Input, double P[2], the end point of the ray to be checked.
            //    The ray is assumed to have an origin at P2.
            //
            //    Output, bool ANGLE_CONTAINS_RAY_2D, is true if the ray is inside
            //    the angle or on its boundary, and false otherwise.
            //
        {
            double a1;
            double a2;
            bool value;

            a1 = angle_deg_2d(p1, p2, p);

            a2 = angle_deg_2d(p1, p2, p3);

            if (a1 <= a2)
            {
                value = true;
            }
            else
            {
                value = false;
            }

            return value;
        }

        public static double angle_deg_2d(double[] p1, double[] p2, double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_DEG_2D returns the angle in degrees swept out between two rays in 2D.
            //
            //  Discussion:
            //
            //    Except for the zero angle case, it should be true that
            //
            //    ANGLE_DEG_2D(P1,P2,P3) + ANGLE_DEG_2D(P3,P2,P1) = 360.0
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
            //    Input, double P1[2], P2[2], P3[2], define the rays
            //    P1 - P2 and P3 - P2 which define the angle.
            //
            //    Output, double ANGLE_DEG_2D, the angle swept out by the rays, measured
            //    in degrees.  0 <= ANGLE_DEG_2D < 360.  If either ray has zero length,
            //    then ANGLE_DEG_2D is set to 0.
            //
        {
            int DIM_NUM = 2;

            double angle_rad;
            double[] p = new double[DIM_NUM];
            double r8_pi = 3.141592653589793;
            double value;

            p[0] = (p1[0] - p2[0]) * (p3[0] - p2[0])
                   + (p1[1] - p2[1]) * (p3[1] - p2[1]);

            p[1] = (p1[0] - p2[0]) * (p3[1] - p2[1])
                   - (p1[1] - p2[1]) * (p3[0] - p2[0]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                value = 0.0;
            }
            else
            {
                angle_rad = Math.Atan2(p[1], p[0]);

                if (angle_rad < 0.0)
                {
                    angle_rad = angle_rad + 2.0 * r8_pi;
                }

                value = Helpers.radians_to_degrees(angle_rad);

            }

            return value;
        }

        public static double[] angle_half_2d(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_HALF_2D finds half an angle in 2D.
            //
            //  Discussion:
            //
            //    The original angle is defined by the sequence of points P1, P2 and P3.
            //
            //    The point P4 is calculated so that:
            //
            //      (P1,P2,P4) = (P1,P2,P3) / 2
            //
            //        P1
            //        /
            //       /   P4
            //      /  .
            //     / .
            //    P2--------->P3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 May 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], P3[2], points defining the angle.
            //
            //    Input, double ANGLE_HALF_2D[2], a point P4 defining the half angle.
            //    The vector P4 - P2 will have unit norm.
            //
        {
            int i;
            double norm;
            double[] p4;

            p4 = new double[2];

            norm = Math.Sqrt((p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length]) * (p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length])
                             + (p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length]) * (p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length]));

            for (i = 0; i < 2; i++)
            {
                p4[i] = (p1[(i + p1Index) % p1.Length] - p2[(i + p2Index) % p2.Length]) / norm;
            }

            norm = Math.Sqrt((p3[(0 + p3Index) % p3.Length] - p2[(0 + p2Index) % p2.Length]) * (p3[(0 + p3Index) % p3.Length] - p2[(0 + p2Index) % p2.Length])
                             + (p3[(1 + p3Index) % p3.Length] - p2[(1 + p2Index) % p2.Length]) * (p3[(1 + p3Index) % p3.Length] - p2[(1 + p2Index) % p2.Length]));

            for (i = 0; i < 2; i++)
            {
                p4[i] = p4[i] + (p3[(i + p3Index) % p3.Length] - p2[(i + p2Index) % p2.Length]) / norm;
            }

            for (i = 0; i < 2; i++)
            {
                p4[i] = 0.5 * p4[i];
            }

            norm = typeMethods.r8vec_norm(2, p4);

            for (i = 0; i < 2; i++)
            {
                p4[i] = p2[(i + p2Index) % p2.Length] + p4[i] / norm;
            }

            return p4;
        }

        public static double angle_rad_2d(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
            //
            //  Discussion:
            //
            //      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
            //
            //        P1
            //        /
            //       /
            //      /
            //     /
            //    P2--------->P3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], P3[2], define the rays
            //    P1 - P2 and P3 - P2 which define the angle.
            //
            //    Output, double ANGLE_RAD_2D, the angle between the two rays,
            //    in radians.  This value will always be between 0 and 2*PI.  If either 
            //    ray has zero length, then the angle is returned as zero.
            //
        {
            double[] p = new double[2];
            double value;

            p[0] = (p3[(0 + p3Index) % p3.Length] - p2[(0 + p2Index) % p2.Length]) * (p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length])
                   + (p3[(1 + p3Index) % p3.Length] - p2[(1 + p2Index) % p2.Length]) * (p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length]);


            p[1] = (p3[(0 + p3Index) % p3.Length] - p2[(0 + p2Index) % p2.Length]) * (p1[(1 + p1Index) % p1.Length] - p2[(1 + p2Index) % p2.Length])
                   - (p3[(1 + p3Index) % p3.Length] - p2[(1 + p2Index) % p2.Length]) * (p1[(0 + p1Index) % p1.Length] - p2[(0 + p2Index) % p2.Length]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                value = 0.0;
                return value;
            }

            value = Math.Atan2(p[1], p[0]);

            if (value < 0.0)
            {
                value = value + 2.0 * Math.PI;
            }

            return value;
        }

        public static double angle_rad_3d(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_3D returns the angle between two vectors in 3D.
            //
            //  Discussion:
            //
            //    The routine always computes the SMALLER of the two angles between
            //    two vectors.  Thus, if the vectors make an (exterior) angle of 200
            //    degrees, the (interior) angle of 160 is reported.
            //
            //    X dot Y = Norm(X) * Norm(Y) * Cos ( Angle(X,Y) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[3], P2[3], P3[3], points defining an angle.
            //    The rays are P1 - P2 and P3 - P2.
            //
            //    Output, double ANGLE_RAD_3D, the angle between the two vectors, in radians.
            //    This value will always be between 0 and PI.  If either vector has
            //    zero length, then the angle is returned as zero.
            //
        {
            int DIM_NUM = 3;

            double dot;
            int i;
            double v1norm;
            double v2norm;
            double value;

            v1norm = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                v1norm = v1norm + Math.Pow(p1[(i + p1Index) % p1.Length] - p2[(i + p2Index) % p2.Length], 2);
            }

            v1norm = Math.Sqrt(v1norm);

            if (v1norm == 0.0)
            {
                value = 0.0;
                return value;
            }

            v2norm = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                v2norm = v2norm + Math.Pow(p3[(i + p3Index) % p3.Length] - p2[(i + p2Index) % p2.Length], 2);
            }

            v2norm = Math.Sqrt(v2norm);

            if (v2norm == 0.0)
            {
                value = 0.0;
                return value;
            }

            dot = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                dot = dot + (p1[(i + p1Index) % p1.Length] - p2[(i + p2Index) % p2.Length]) * (p3[(i + p3Index) % p3.Length] - p2[(i + p2Index) % p2.Length]);
            }

            value = typeMethods.r8_acos(dot / (v1norm * v2norm));

            return value;
        }

        public static double angle_rad_nd(int dim_num, double[] vec1, double[] vec2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_ND returns the angle between two vectors in ND.
            //
            //  Discussion:
            //
            //    ANGLE_RAD_ND always computes the SMALLER of the two angles between
            //    two vectors.  Thus, if the vectors make an (exterior) angle of
            //    1.5 PI radians, then the (interior) angle of 0.5 radians is returned.
            //
            //    X dot Y = Norm(X) * Norm(Y) * Cos( Angle(X,Y) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double VEC1[DIM_NUM], VEC2[DIM_NUM], the two vectors to be considered.
            //
            //    Output, double ANGLE_RAD_ND, the angle between the vectors, in radians.
            //    This value will always be between 0 and PI.
            //
        {
            double dot;
            double v1norm;
            double v2norm;
            double value;

            dot = typeMethods.r8vec_dot_product(dim_num, vec1, vec2);

            v1norm = typeMethods.r8vec_norm(dim_num, vec1);
            v2norm = typeMethods.r8vec_norm(dim_num, vec2);

            if (v1norm == 0.0 || v2norm == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Acos(dot / (v1norm * v2norm));
            }

            return value;
        }

        public static double angle_turn_2d(double[] p1, double[] p2, double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_TURN_2D computes a turning angle in 2D.
            //
            //  Discussion:
            //
            //    This routine is most useful when considering the vertices of a
            //    polygonal shape.  We wish to distinguish between angles that "turn
            //    in" to the shape, (between 0 and 180 degrees) and angles that
            //    "turn out" (between 180 and 360 degrees), as we traverse the boundary.
            //
            //    If we compute the interior angle and subtract 180 degrees, we get the
            //    supplementary angle, which has the nice property that it is
            //    negative for "in" angles and positive for "out" angles, and is zero if
            //    the three points actually lie along a line.
            //
            //    Assuming P1, P2 and P3 define an angle, the TURN can be
            //    defined to be either:
            //
            //    * the supplementary angle to the angle formed by P1=P2=P3, or
            //
            //    * the angle between the vector ( P3-P2) and the vector -(P1-P2),
            //      where -(P1-P2) can be understood as the vector that continues
            //      through P2 from the direction P1.
            //
            //    The turning will be zero if P1, P2 and P3 lie along a straight line.
            //
            //    It will be a positive angle if the turn from the previous direction
            //    is counter clockwise, and negative if it is clockwise.
            //
            //    The turn is given in radians, and will lie between -PI and PI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], P2[2], P3[2], the points that form
            //    the angle.
            //
            //    Output, double ANGLE_TURN_2D, the turn angle, between -PI and PI.
            //
        {
            int DIM_NUM = 2;

            double[] p = new double[DIM_NUM];
            double r8_pi = 3.141592653589793;
            double turn;

            p[0] = (p3[0] - p2[0]) * (p1[0] - p2[0])
                   + (p3[1] - p2[1]) * (p1[1] - p2[1]);

            p[1] = (p3[0] - p2[0]) * (p1[1] - p2[1])
                   - (p3[1] - p2[1]) * (p1[0] - p2[0]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                turn = 0.0;
            }
            else
            {
                turn = r8_pi - typeMethods.r8_atan(p[1], p[0]);
            }

            return turn;
        }

        public static double anglei_deg_2d(double[] p1, double[] p2, double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLEI_DEG_2D returns the interior angle in degrees between two rays in 2D.
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
            //    Input, double P1[3], P2[3], P3[3], points defining an angle.
            //    The rays are P1 - P2 and P3 - P2.
            //
            //    Output, double ANGLEI_DEG_2D, the angle swept out by the rays, measured
            //    in degrees.  This value satisfies 0 <= ANGLEI_DEG_2D < 180.0.  If either
            //    ray is of zero length, then ANGLEI_deg_2D is returned as 0.
            //
        {
            int DIM_NUM = 2;

            double[] p = new double[DIM_NUM];
            double r8_pi = 3.141592653589793;
            double value;

            p[0] = (p1[0] - p2[0]) * (p3[0] - p2[0])
                   + (p1[1] - p2[1]) * (p3[1] - p2[1]);

            p[1] = (p1[0] - p2[0]) * (p3[1] - p2[1])
                   - (p1[1] - p2[1]) * (p3[0] - p2[0]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Atan2(p[1], p[0]);

                if (value < 0.0)
                {
                    value = value + 2.0 * r8_pi;
                }

                value = Helpers.radians_to_degrees(value);

                if (180.0 < value)
                {
                    value = 360.0 - value;
                }

            }

            return value;
        }

        public static double anglei_rad_2d(double[] p1, double[] p2, double[] p3)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLEI_RAD_2D returns the interior angle in radians between two rays in 2D.
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
            //    Input, double P1[3], P2[3], P3[3], points defining an angle.
            //    The rays are P1 - P2 and P3 - P2.
            //
            //    Output, double ANGLEI_RAD_2D, the angle swept out by the rays, measured
            //    in degrees.  This value satisfies 0 <= ANGLEI_RAD_2D < PI.  If either
            //    ray is of zero length, then ANGLEI_RAD_2D is returned as 0.
            //
        {
            int DIM_NUM = 2;

            double[] p = new double[DIM_NUM];
            double r8_pi = 3.141592653589793;
            double value;

            p[0] = (p1[0] - p2[0]) * (p3[0] - p2[0])
                   + (p1[1] - p2[1]) * (p3[1] - p2[1]);

            p[1] = (p1[0] - p2[0]) * (p3[1] - p2[1])
                   - (p1[1] - p2[1]) * (p3[0] - p2[0]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Atan2(p[1], p[0]);

                if (value < 0.0)
                {
                    value = value + 2.0 * r8_pi;
                }

                if (r8_pi < value)
                {
                    value = 2.0 * r8_pi - value;
                }

            }

            return value;
        }

    }
}