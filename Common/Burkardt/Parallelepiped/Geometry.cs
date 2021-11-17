using System;

namespace Burkardt.Parallelepiped;

public class Geometry
{
    public static bool parallelepiped_contains_point_3d(double[] p1, double[] p2, double[] p3,
            double[] p4, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELEPIPED_CONTAINS_POINT_3D determines if a point is inside a parallelepiped in 3D.
        //
        //  Discussion:
        //
        //    A parallelepiped is a "slanted box", that is, opposite
        //    sides are parallel planes.
        //
        //         *------------------*
        //        / .                / .
        //       /   .              /   .
        //      /     .            /     .
        //    P4------------------*       .
        //      .        .         .       .
        //       .        .         .       .
        //        .        .         .       .
        //         .       P2..........-------.
        //          .     /            .     .
        //           .   /              .   .
        //            . /                . .
        //             P1----------------P3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], P4[3], four vertices of the parallelepiped.
        //    It is assumed that P2, P3 and P4 are immediate neighbors of P1.
        //
        //    Input, double P, the point to be checked.
        //
        //    Output, bool PARAPP_CONTAINS_POINT_3D, is true if P is inside the
        //    parallelepiped, or on its boundary, and false otherwise.
        //
    {
        double dot;

        dot = (p2[0] - p1[0]) * (p[0] - p1[0])
              + (p2[1] - p1[1]) * (p[1] - p1[1])
              + (p2[2] - p1[2]) * (p[2] - p1[2]);

        switch (dot)
        {
            case < 0.0:
                return false;
        }

        if (Math.Pow(p2[0] - p1[0], 2)
            + Math.Pow(p2[1] - p1[1], 2)
            + Math.Pow(p2[2] - p1[2], 2) < dot)
        {
            return false;
        }

        dot = (p3[0] - p1[0]) * (p[0] - p1[0])
              + (p3[1] - p1[1]) * (p[1] - p1[1])
              + (p3[2] - p1[2]) * (p[2] - p1[2]);

        switch (dot)
        {
            case < 0.0:
                return false;
        }

        if (Math.Pow(p3[0] - p1[0], 2)
            + Math.Pow(p3[1] - p1[1], 2)
            + Math.Pow(p3[2] - p1[2], 2) < dot)
        {
            return false;
        }

        dot = (p4[0] - p1[0]) * (p[0] - p1[0])
              + (p4[1] - p1[1]) * (p[1] - p1[1])
              + (p4[2] - p1[2]) * (p[2] - p1[2]);

        switch (dot)
        {
            case < 0.0:
                return false;
        }

        if (Math.Pow(p4[0] - p1[0], 2)
            + Math.Pow(p4[1] - p1[1], 2)
            + Math.Pow(p4[2] - p1[2], 2) < dot)
        {
            return false;
        }

        return true;
    }

    public static double parallelepiped_point_dist_3d(double[] p1, double[] p2, double[] p3,
            double[] p4, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELEPIPED_POINT_DIST_3D: distance ( parallelepiped, point ) in 3D.
        //
        //  Discussion:
        //
        //    A parallelepiped is a "slanted box", that is, opposite
        //    sides are parallel planes.
        //
        //    A parallelepiped is a "slanted box", that is, opposite
        //    sides are parallel planes.
        //
        //         *------------------*
        //        / .                / .
        //       /   .              /   .
        //      /     .            /     .
        //    P4------------------*       .
        //      .        .         .       .
        //       .        .         .       .
        //        .        .         .       .
        //         .       P2..........-------.
        //          .     /            .     /
        //           .   /              .   /
        //            . /                . /
        //             P1----------------P3
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
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], P4[3], half of
        //    the corners of the box, from which the other corners can be
        //    deduced.  The corners should be chosen so that the first corner
        //    is directly connected to the other three.  The locations of
        //    corners 5, 6, 7 and 8 will be computed by the parallelogram
        //    relation.
        //
        //    Input, double P[3], the point which is to be checked.
        //
        //    Output, double PARAPP_POINT_DIST_3D, the distance from the point to the box.
        //    The distance is zero if the point lies exactly on the box.
        //
    {
        int DIM_NUM = 3;

        double dis;
        double dist;
        double[] p5 = new double[DIM_NUM];
        double[] p6 = new double[DIM_NUM];
        double[] p7 = new double[DIM_NUM];
        double[] p8 = new double[DIM_NUM];
        //
        //  Fill in the other corners
        //
        p5[0] = p2[0] + p3[0] - p1[0];
        p5[1] = p2[1] + p3[1] - p1[1];
        p5[2] = p2[2] + p3[2] - p1[2];

        p6[0] = p2[0] + p4[0] - p1[0];
        p6[1] = p2[1] + p4[1] - p1[1];
        p6[2] = p2[2] + p4[2] - p1[2];

        p7[0] = p3[0] + p4[0] - p1[0];
        p7[1] = p3[1] + p4[1] - p1[1];
        p7[2] = p3[2] + p4[2] - p1[2];

        p8[0] = p2[0] + p3[0] + p4[0] - 2.0 * p1[0];
        p8[1] = p2[1] + p3[1] + p4[1] - 2.0 * p1[1];
        p8[2] = p2[2] + p3[2] + p4[2] - 2.0 * p1[2];
        //
        //  Compute the distance from the point P to each of the six
        //  paralleogram faces.
        //
        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p1, p2, p3, p);

        dist = dis;

        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p1, p2, p4, p);

        if (dis < dist)
        {
            dist = dis;
        }

        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p1, p3, p4, p);

        if (dis < dist)
        {
            dist = dis;
        }

        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p8, p5, p6, p);

        if (dis < dist)
        {
            dist = dis;
        }

        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p8, p5, p7, p);

        if (dis < dist)
        {
            dist = dis;
        }

        dis = Parallelogram.Geometry.parallelogram_point_dist_3d(p8, p6, p7, p);

        if (dis < dist)
        {
            dist = dis;
        }

        return dist;
    }
}