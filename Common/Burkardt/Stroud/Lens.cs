﻿using System;
using Burkardt.Quadrature;

namespace Burkardt.Stroud;

public static class Lens
{

    public static double lens_half_2d(Func<double, double, double> func, double[] center,
            double r, double theta1, double theta2, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LENS_HALF_2D approximates an integral in a circular half lens in 2D.
        //
        //  Discussion:
        //
        //    A circular half lens is formed by drawing a circular arc,
        //    and joining its endpoints.
        //
        //    This rule for a circular half lens simply views the region as 
        //    a product region, with a coordinate "S" that varies along the
        //    radial direction, and a coordinate "T" that varies in the perpendicular
        //    direction, and whose extent varies as a function of S.
        //
        //    A Gauss-Legendre rule is used to construct a product rule that is
        //    applied to the region.  The accuracy of the Gauss-Legendre rule,
        //    which is valid for a rectangular product rule region, does not
        //    apply straightforwardly to this region, since the limits in the
        //    "T" coordinate are being handled implicitly.
        //
        //    This is simply an application of the QMULT_2D algorithm of Stroud.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, Func < double, double, double > func, the name of the 
        //    user supplied function to be integrated.
        //
        //    Input, double[] center, the center of the circle.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double THETA1, THETA2, the angles of the rays
        //    that begin and end the arc.
        //
        //    Input, int ORDER, the order of the Gauss-Legendre rule
        //    to be used.  Legal values include 1 through 20, 32 or 64.
        //
        //    Output, double LENS_HALF_2D, the approximate value
        //    of the integral of the function over the half lens.
        //
    {
        int i;
        double quad;
        //
        //  Determine the points A (on the secant) and B (on the circumference)
        //  that will form the "S" direction.
        //
        double ax = center[0] + r * 0.5 * (Math.Cos(theta1) + Math.Cos(theta2));
        double ay = center[1] + r * 0.5 * (Math.Sin(theta1) + Math.Sin(theta2));

        double bx = center[0] + r * Math.Cos(0.5 * (theta1 + theta2));
        double by = center[1] + r * Math.Sin(0.5 * (theta1 + theta2));
        //
        //  Find the length of the line between A and B.
        //
        double s_length = Math.Sqrt(Math.Pow(ax - bx, 2) + Math.Pow(ay - by, 2));

        switch (s_length)
        {
            case 0.0:
                quad = 0.0;
                return quad;
        }

        //
        //  Retrieve the Legendre rule of the given order.
        //
        double[] xtab = new double[order];
        double[] weight = new double[order];

        LegendreQuadrature.legendre_set(order, ref xtab, ref weight);
        //
        //  Determine the unit vector in the T direction.
        //
        double tdirx = (ay - by) / s_length;
        double tdiry = (bx - ax) / s_length;

        quad = 0.0;

        for (i = 0; i < order; i++)
        {
            double w1 = 0.5 * s_length * weight[i];
            //
            //  Map the quadrature point to an S coordinate.
            //
            double sx = ((1.0 - xtab[i]) * ax
                         + (1.0 + xtab[i]) * bx)
                        / 2.0;
            double sy = ((1.0 - xtab[i]) * ay
                         + (1.0 + xtab[i]) * by)
                        / 2.0;
            //
            //  Determine the length of the line in the T direction, from the
            //  S axis to the circle circumference.
            //
            double thi = Math.Sqrt((r - 0.25 * (1.0 - xtab[i]) * s_length)
                                   * (1.0 - xtab[i]) * s_length);
            // 
            //  Determine the maximum and minimum T coordinates by going
            //  up and down in the T direction from the S axis.
            //
            double cx = sx + tdirx * thi;
            double cy = sy + tdiry * thi;
            double dx = sx - tdirx * thi;
            double dy = sy - tdiry * thi;
            //
            //  Find the length of the T direction.
            //
            double t_length = Math.Sqrt(Math.Pow(cx - dx, 2) + Math.Pow(cy - dy, 2));

            int j;
            for (j = 0; j < order; j++)
            {
                double w2 = 0.5 * t_length * weight[j];
                //
                //  Map the quadrature point to a T coordinate.
                //
                double tx = ((1.0 - xtab[j]) * cx
                             + (1.0 + xtab[j]) * dx)
                            / 2.0;
                double ty = ((1.0 - xtab[j]) * cy
                             + (1.0 + xtab[j]) * dy)
                            / 2.0;

                quad += w1 * w2 * func(tx, ty);
            }
        }

        return quad;
    }

    public static double lens_half_area_2d(double r, double theta1, double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LENS_HALF_AREA_2D returns the area of a circular half lens in 2D.
        //
        //  Discussion:
        //
        //    A circular half lens is formed by drawing a circular arc, 
        //    and joining its endpoints.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double THETA1, THETA2, the angles of the rays
        //    that begin and end the arc.
        //
        //    Output, double LENS_HALF_AREA_2D, the area of the half lens.
        //
    {
        double sector = Circle.circle_sector_area_2d(r, theta1, theta2);
        double triangle = Circle.circle_triangle_area_2d(r, theta1, theta2);
        double value = sector - triangle;

        return value;
    }

    public static double lens_half_h_area_2d(double r, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LENS_HALF_H_AREA_2D returns the area of a circular half lens in 2D.
        //
        //  Discussion:
        //
        //    A circular half lens is formed by drawing a circular arc, and joining 
        //    its endpoints.
        //
        //    This particular half lens is described by the "height" of the region.  
        //    In other words, the half lens is the region that would be submerged 
        //    if a circle of radius R were standing in water of depth H.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double H, the height of the half lens region.
        //
        //    Output, double LENS_HALF_H_AREA_2D, the area of the half lens.
        //
    {
        double area;

        switch (h)
        {
            case <= 0.0:
                area = 0.0;
                break;
            default:
            {
                if (2.0 * r <= h)
                {
                    area = Math.PI * r * r;
                }
                else
                {
                    double half_width = Math.Sqrt(h * (2.0 * r - h));
                    double angle = 2.0 * Math.Atan2(half_width, r - h);
                    double sector = r * r * angle / 2.0;
                    double triangle = (r - h) * half_width;
                    area = sector - triangle;
                }

                break;
            }
        }

        return area;
    }

    public static double lens_half_w_area_2d(double r, double w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LENS_HALF_W_AREA_2D returns the area of a circular half lens in 2D.
        //
        //  Discussion:
        //
        //    A half lens is formed by drawing a circular arc, and joining its endpoints.
        //    This half lens is described by the "width" of the region.  In other words,
        //    it is the portion of the circle under water if the width
        //    of the water surface is W.  There are two possible values for this
        //    area, A and (PI*R*R-A).  The routine returns the smaller of the 
        //    two values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double W, the width of the half lens region.
        //
        //    Output, double LENS_HALF_W_AREA_2D, the area of the half lens.
        //
    {
        double area;

        switch (w)
        {
            case <= 0.0:
                area = 0.0;
                break;
            default:
            {
                if (2.0 * r <= w)
                {
                    area = 0.5 * Math.PI * r * r;
                }
                else
                {
                    double half_width = 0.5 * w;
                    double h = r - Math.Sqrt(r * r - half_width * half_width);
                    double angle = 2.0 * Math.Atan2(half_width, r - h);
                    double sector = r * r * angle / 2.0;
                    double triangle = (r - h) * half_width;
                    area = sector - triangle;
                }

                break;
            }
        }

        return area;
    }


}