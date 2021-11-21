﻿using System;

namespace Burkardt.AnnulusNS;

public static class Annulus
{
    public static double annulus_area_2d(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANNULUS_AREA_2D computes the area of a circular annulus in 2D.
        //
        //  Discussion:
        //
        //    A circular annulus with center (XC,YC), inner radius R1 and
        //    outer radius R2, is the set of points (X,Y) so that
        //
        //      R1^2 <= (X-XC)^2 + (Y-YC)^2 <= R2^2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Output, double ANNULUS_AREA_2D, the area.
        //
    {
        double area = Math.PI * (r2 + r1) * (r2 - r1);

        return area;
    }

    public static double annulus_sector_area_2d(double r1, double r2, double theta1,
            double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANNULUS_SECTOR_AREA_2D computes the area of an annular sector in 2D.
        //
        //  Discussion:
        //
        //    An annular sector with center PC, inner radius R1 and
        //    outer radius R2, and angles THETA1, THETA2, is the set of points
        //    P so that
        //
        //      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
        //
        //    and
        //
        //      THETA1 <= THETA ( P - PC ) <= THETA2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Input, double THETA1, THETA2, the angles.
        //
        //    Output, double ANNULUS_SECTOR_AREA_2D, the area.
        //
    {
        double area = 0.5 * (theta2 - theta1) * (r2 + r1) * (r2 - r1);

        return area;
    }

    public static double[] annulus_sector_centroid_2d(double[] pc, double r1, double r2,
            double theta1, double theta2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANNULUS_SECTOR_CENTROID_2D computes the centroid of an annular sector in 2D.
        //
        //  Discussion:
        //
        //    An annular sector with center PC, inner radius R1 and
        //    outer radius R2, and angles THETA1, THETA2, is the set of points
        //    P so that
        //
        //      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
        //
        //    and
        //
        //      THETA1 <= THETA ( P - PC ) <= THETA2
        //
        //    Thanks to Ed Segall for pointing out a mistake in the computation
        //    of the angle THETA assciated with the centroid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 December 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Harris, Horst Stocker,
        //    Handbook of Mathematics and Computational Science,
        //    Springer, 1998, QA40.S76
        //
        //  Parameters:
        //
        //    Input, double PC[2], the center.
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Input, double THETA1, THETA2, the angles.
        //
        //    Output, double ANNULUS_SECTOR_CENTROID_2D[2], the centroid.
        //
    {
        double theta = theta2 - theta1;

        double r = 4.0 * Math.Sin(theta / 2.0) / (3.0 * theta)
            * (r1 * r1 + r1 * r2 + r2 * r2) / (r1 + r2);

        double[] centroid = new double[2];

        centroid[0] = pc[0] + r * Math.Cos(theta1 + theta / 2.0);
        centroid[1] = pc[1] + r * Math.Sin(theta1 + theta / 2.0);

        return centroid;
    }
}