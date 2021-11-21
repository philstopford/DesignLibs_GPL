﻿using System;
using Burkardt.Types;

namespace Burkardt.Ellipse;

public static class Geometry
{

    public static double ellipse_area1(double[] a, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA1 returns the area of an ellipse defined by a matrix.
        //
        //  Discussion:
        //
        //    The points X in the ellipse are described by a 2 by 2
        //    positive definite symmetric matrix A, and a "radius" R, such that
        //      X' * A * X <= R * R
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2*2], the matrix that describes
        //    the ellipse.  A must be symmetric and positive definite.
        //
        //    Input, double R, the "radius" of the ellipse.
        //
        //    Output, double ELLIPSE_AREA1, the area of the ellipse.
        //
    {
        double value = 0;

        value = r * r * Math.PI / Math.Sqrt(a[0 + 0 * 2] * a[1 + 1 * 2] - a[1 + 0 * 2] * a[0 + 1 * 2]);

        return value;
    }

    public static double ellipse_area2(double a, double b, double c, double d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA2 returns the area of an ellipse defined by an equation.
        //
        //  Discussion:
        //
        //    The ellipse is described by the formula
        //      a x^2 + b xy + c y^2 = d
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, coefficients on the left hand side.
        //
        //    Input, double D, the right hand side.
        //
        //    Output, double ELLIPSE_AREA2, the area of the ellipse.
        //
    {
        double value = 0;

        value = 2.0 * d * d * Math.PI / Math.Sqrt(4.0 * a * c - b * b);

        return value;
    }

    public static double ellipse_area3(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA3 returns the area of an ellipse in 2D.
        //
        //  Discussion:
        //
        //    An ellipse in standard position has a center at the origin, and
        //    axes aligned with the coordinate axes.  Any point P on the ellipse
        //    satisfies
        //
        //      Math.Pow (  P[0] / R1, 2 ) + Math.Pow ( P[1] / R2, 2 ) == 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the "radius" of the ellipse in the major
        //    and minor axis directions.  A circle has these values equal.
        //
        //    Output, double ELLIPSE_AREA3, the area of the ellipse.
        //
    {
        double area = Math.PI * r1 * r2;

        return area;
    }

    public static double ellipse_point_dist_2d(double r1, double r2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINT_DIST_2D finds the distance from a point to an ellipse in 2D.
        //
        //  Discussion:
        //
        //    An ellipse in standard position has a center at the origin, and
        //    axes aligned with the coordinate axes.  Any point P on the ellipse
        //    satisfies
        //
        //      Math.Pow (  P[0] / R1, 2 ) + Math.Pow ( P[1] / R2, 2 ) == 1
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
        //  Reference:
        //
        //    Dianne O'Leary,
        //    Elastoplastic Torsion: Twist and Stress,
        //    Computing in Science and Engineering,
        //    July/August 2004, pages 74-76.
        //    September/October 2004, pages 63-65.
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the ellipse parameters.  Normally,
        //    these are both positive quantities.  Generally, they are also
        //    distinct.
        //
        //    Input, double P[2], the point.
        //
        //    Output, double ELLIPSE_POINT_DIST_2D, the distance to the ellipse.
        //
    {
        const int DIM_NUM = 2;

        int i;

        double[] pn = ellipse_point_near_2d(r1, r2, p);

        double dist = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            dist += Math.Pow(p[i] - pn[i], 2);
        }

        dist = Math.Sqrt(dist);

        return dist;
    }

    public static double[] ellipse_point_near_2d(double r1, double r2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINT_NEAR_2D finds the nearest point on an ellipse in 2D.
        //
        //  Discussion:
        //
        //    An ellipse in standard position has a center at the origin, and
        //    axes aligned with the coordinate axes.  Any point P on the ellipse
        //    satisfies
        //
        //      (  P(1) / R1 )^2 + ( P(2) / R2 )^2 == 1
        //
        //    The nearest point PN on the ellipse has the property that the
        //    line from PN to P is normal to the ellipse.  Points on the ellipse
        //    can be parameterized by T, to have the form
        //
        //      ( R1 * Math.Cos ( T ), R2 * Math.Sin ( T ) ).
        //
        //    The tangent vector to the ellipse has the form
        //
        //      ( -R1 * Math.Sin ( T ), R2 * Math.Cos ( T ) )
        //
        //    At PN, the dot product of this vector with  ( P - PN ) must be
        //    zero:
        //
        //      - R1 * Math.Sin ( T ) * ( X - R1 * Math.Cos ( T ) )
        //      + R2 * Math.Cos ( T ) * ( Y - R2 * Math.Sin ( T ) ) = 0
        //
        //    This nonlinear equation for T can be solved by Newton's method.
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
        //    Input, double R1, R2, the ellipse parameters.  Normally,
        //    these are both positive quantities.  Generally, they are also
        //    distinct.
        //
        //    Input, double P[2], the point.
        //
        //    Output, double ELLIPSE_POINT_NEAR_2D[2], the point on the ellipse which
        //    is closest to P.
        //
    {
        const int DIM_NUM = 2;

        const int iteration_max = 100;
        double t;

        double x = Math.Abs(p[0]);
        double y = Math.Abs(p[1]);

        switch (y)
        {
            case 0.0 when r1 * r1 - r2 * r2 <= r1 * x:
                t = 0.0;
                break;
            default:
            {
                switch (x)
                {
                    case 0.0 when r2 * r2 - r1 * r1 <= r2 * y:
                        t = Math.PI / 2.0;
                        break;
                    default:
                    {
                        y = y switch
                        {
                            0.0 => Math.Sqrt(typeMethods.r8_epsilon()) * Math.Abs(r2),
                            _ => y
                        };

                        x = x switch
                        {
                            0.0 => Math.Sqrt(typeMethods.r8_epsilon()) * Math.Abs(r1),
                            _ => x
                        };

                        //
                        //  Initial parameter T:
                        //
                        t = Math.Atan2(y, x);

                        int iteration = 0;

                        for (;;)
                        {
                            double ct = Math.Cos(t);
                            double st = Math.Sin(t);

                            double f = (x - Math.Abs(r1) * ct) * Math.Abs(r1) * st
                                       - (y - Math.Abs(r2) * st) * Math.Abs(r2) * ct;

                            if (Math.Abs(f) <= 100.0 * typeMethods.r8_epsilon())
                            {
                                break;
                            }

                            if (iteration_max <= iteration)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("ELLIPSE_POINT_NEAR_2D - Warning!");
                                Console.WriteLine("  Reached iteration limit.");
                                Console.WriteLine("  T = " + t + "");
                                Console.WriteLine("  F = " + f + "");
                                break;
                            }

                            iteration += 1;

                            double fp = r1 * r1 * st * st + r2 * r2 * ct * ct
                                                          + (x - Math.Abs(r1) * ct) * Math.Abs(r1) * ct
                                                          + (y - Math.Abs(r2) * st) * Math.Abs(r2) * st;

                            t -= f / fp;
                        }

                        break;
                    }
                }

                break;
            }
        }

        //
        //  From the T value, we get the nearest point.
        //
        double[] pn = new double[DIM_NUM];

        pn[0] = Math.Abs(r1) * Math.Cos(t);
        pn[1] = Math.Abs(r2) * Math.Sin(t);
        //
        //  Take care of case where the point was in another quadrant.
        //
        pn[0] = typeMethods.r8_sign(p[0]) * pn[0];
        pn[1] = typeMethods.r8_sign(p[1]) * pn[1];

        return pn;
    }

    public static void ellipse_points_2d(double[] pc, double r1, double r2, double psi,
            int n, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINTS_2D returns N points on an tilted ellipse in 2D.
        //
        //  Discussion:
        //
        //    An ellipse in standard position has a center at the origin, and
        //    axes aligned with the coordinate axes.  Any point P on the ellipse
        //    satisfies
        //
        //      Math.Pow (  P[0] / R1, 2 ) + Math.Pow ( P[1] / R2, 2 ) == 1
        //
        //    The points are "equally spaced" in the angular sense.  They are
        //    not equally spaced along the perimeter of the ellipse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double PC[2], the coordinates of the center of the ellipse.
        //
        //    Input, double R1, R2, the "radius" of the ellipse in the major
        //    and minor axis directions.  A circle has these values equal.
        //
        //    Input, double PSI, the angle that the major axis of the ellipse
        //    makes with the X axis.  A value of 0.0 means that the major and
        //    minor axes of the ellipse will be the X and Y coordinate axes.
        //
        //    Input, int N, the number of points desired.  N must be at least 1.
        //
        //    Output, double P[2*N], the coordinates of points on the ellipse.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            double theta = 2.0 * Math.PI * i / n;

            p[0 + i * 2] = pc[0] + r1 * Math.Cos(psi) * Math.Cos(theta)
                           - r2 * Math.Sin(psi) * Math.Sin(theta);

            p[1 + i * 2] = pc[1] + r1 * Math.Sin(psi) * Math.Cos(theta)
                                 + r2 * Math.Cos(psi) * Math.Sin(theta);

        }

    }

    public static void ellipse_points_arc_2d(double[] pc, double r1, double r2, double psi,
            double theta1, double theta2, int n, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINTS_ARC_2D returns N points on a tilted elliptical arc in 2D.
        //
        //  Discussion:
        //
        //    An ellipse in standard position has a center at the origin, and
        //    axes aligned with the coordinate axes.  Any point P on the ellipse
        //    satisfies
        //
        //      Math.Pow (  P[0] / R1, 2 ) + Math.Pow ( P[1] / R2, 2 ) == 1
        //
        //    The points are "equally spaced" in the angular sense.  They are
        //    not equally spaced along the perimeter of the ellipse.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double PC[2], the coordinates of the center of the ellipse.
        //
        //    Input, double R1, R2, the "radius" of the ellipse in the major
        //    and minor axis directions.  A circle has these values equal.
        //
        //    Input, double PSI, the angle that the major axis of the ellipse
        //    makes with the X axis.  A value of 0.0 means that the major and
        //    minor axes of the ellipse will be the X and Y coordinate axes.
        //
        //    Input, double THETA1, THETA2, the angular coordinates of the first
        //    and last points to be drawn, in radians.  This angle is measured
        //    with respect to the (possibly tilted) major axis.
        //
        //    Input, int N, the number of points desired.  N must be at least 1.
        //
        //    Output, double P[2*N], the coordinates of points on the ellipse.
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

            p[0 + i * 2] = pc[0] + r1 * Math.Cos(psi) * Math.Cos(theta)
                           - r2 * Math.Sin(psi) * Math.Sin(theta);

            p[1 + i * 2] = pc[1] + r1 * Math.Sin(psi) * Math.Cos(theta)
                                 + r2 * Math.Cos(psi) * Math.Sin(theta);

        }

    }
        
    public static void super_ellipse_points_2d ( double[] pc, double r1, double r2,
            double expo, double psi, int n, ref double[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUPER_ELLIPSE_POINTS_2D returns N points on a tilted superellipse in 2D.
        //
        //  Discussion:
        //
        //    The points are "equally spaced" in the angular sense.  They are
        //    not equally spaced along the perimeter.
        //
        //    The parametric formula of the (untilted) superellipse is:
        //
        //      X = R1 * Math.Cos**EXPO ( THETA )
        //      Y = R2 * Math.Sin**EXPO ( THETA )
        //
        //    An implicit form of the (untilted) superellipse is:
        //
        //      (X/R1)**(2/EXPO) + (Y/R2)**(2/EXPO) = 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Martin Gardner,
        //    The Mathematical Carnival,
        //    Knopf, 1975, pages 240-254.
        //
        //  Parameters:
        //
        //    Input, double PC[2], the coordinates of the center of the superellipse.
        //
        //    Input, double R1, R2, the "radius" of the superellipse in the major
        //    and minor axis directions.  A circle has these values equal.
        //
        //    Input, double EXPO, the exponent of the superellipse.
        //    0 = a rectangle;
        //    between 0 and 1, a "rounded" rectangle;
        //    1.0 = an ellipse;
        //    2.0 = a diamond;
        //    > 2.0 a pinched shape.
        //
        //    Input, double PSI, the angle that the major axis of the superellipse
        //    makes with the X axis.  A value of 0.0 means that the major and
        //    minor axes of the superellipse will be the X and Y coordinate axes.
        //
        //    Input, int N, the number of points desired.  N must be at least 1.
        //
        //    Output, double P[2*N], the coordinates of points on the superellipse.
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            double theta = 2.0 * Math.PI * i / n;

            double act = Math.Abs(Math.Cos(theta));
            double sct = typeMethods.r8_sign(Math.Cos(theta));
            double ast = Math.Abs(Math.Sin(theta));
            double sst = typeMethods.r8_sign(Math.Sin(theta));

            p[0 + i * 2] = pc[0] + r1 * Math.Cos(psi) * sct * Math.Pow(act, expo)
                           - r2 * Math.Sin(psi) * sst * Math.Pow(ast, expo);

            p[1 + i * 2] = pc[1] + r1 * Math.Sin(psi) * sct * Math.Pow(act, expo)
                                 + r2 * Math.Cos(psi) * sst * Math.Pow(ast, expo);

        }
    }


}