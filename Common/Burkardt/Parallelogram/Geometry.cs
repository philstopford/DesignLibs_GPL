using System;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.Parallelogram;

public static class Geometry
{
    public static double parallelogram_area_2d(double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_AREA_2D computes the area of a parallelogram in 2D.
        //
        //  Discussion:
        //
        //    A parallelogram is a polygon having four sides, with the property
        //    that each pair of opposite sides is paralell.
        //
        //    Given the first three vertices of the parallelogram,
        //    P1, P2, and P3, the fourth vertex must satisfy
        //
        //      P4 = P1 + ( P3 - P2 )
        //
        //    This routine uses the fact that the norm of the cross product
        //    of two vectors is the area of the parallelogram they form:
        //
        //      Area = ( P3 - P2 ) x ( P1 - P2 ).
        //
        //        P4<-----P3
        //        /       /
        //       /       /
        //      P1----->P2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P[2*4], the parallelogram vertices,
        //    given in counterclockwise order.  The fourth vertex is ignored.
        //
        //    Output, double PARALLELOGRAM_AREA_2D, the (signed) area.
        //
    {
        double area;
        //
        //  Compute the cross product vector, which only has a single
        //  nonzero component.
        //
        area = (p[0 + 1 * 2] - p[0 + 0 * 2]) * (p[1 + 2 * 2] - p[1 + 0 * 2])
               - (p[1 + 1 * 2] - p[1 + 0 * 2]) * (p[0 + 2 * 2] - p[0 + 0 * 2]);

        return area;
    }

    public static double parallelogram_area_3d(double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_AREA_3D computes the area of a parallelogram in 3D.
        //
        //  Discussion:
        //
        //    A parallelogram is a polygon having four sides, with the property
        //    that each pair of opposite sides is paralell.
        //
        //    A parallelogram in 3D must have the property that it is "really"
        //    a 2D object, that is, that the four vertices that define it lie
        //    in some plane.
        //
        //    Given the first three vertices of the parallelogram (in 2D or 3D),
        //    P1, P2, and P3, the fourth vertex must satisfy
        //
        //      P4 = P1 + ( P3 - P2 )
        //
        //    This routine uses the fact that the norm of the cross product
        //    of two vectors is the area of the parallelogram they form:
        //
        //      Area = ( P3 - P2 ) x ( P1 - P2 ).
        //
        //        P4<-----P3
        //        /       /
        //       /       /
        //      P1----->P2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P[3*4], the parallelogram vertices,
        //    given in counterclockwise order.  The fourth vertex is ignored.
        //
        //    Output, double PARALLELOGRAM_AREA_3D, the area
        //
    {
        double area;
        double cross;
        //
        //  Compute the cross product vector.
        //
        area = 0.0;

        cross = (p[1 + 1 * 3] - p[1 + 0 * 3]) * (p[2 + 2 * 3] - p[2 + 0 * 3])
                - (p[2 + 1 * 3] - p[2 + 0 * 3]) * (p[1 + 2 * 3] - p[1 + 0 * 3]);

        area += cross * cross;

        cross = (p[2 + 1 * 3] - p[2 + 0 * 3]) * (p[0 + 2 * 3] - p[0 + 0 * 3])
                - (p[0 + 1 * 3] - p[0 + 0 * 3]) * (p[2 + 2 * 3] - p[2 + 0 * 3]);

        area += cross * cross;

        cross = (p[0 + 1 * 3] - p[0 + 0 * 3]) * (p[1 + 2 * 3] - p[1 + 0 * 3])
                - (p[1 + 1 * 3] - p[1 + 0 * 3]) * (p[0 + 2 * 3] - p[0 + 0 * 3]);

        area += cross * cross;

        area = Math.Sqrt(area);

        return area;
    }

    public static bool parallelogram_contains_point_2d(double[] p1, double[] p2, double[] p3,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_CONTAINS_POINT_2D determines if a point is inside a parallelogram in 2D.
        //
        //  Discussion:
        //
        //         P2..............
        //        /              .
        //       /              .
        //      /              .
        //    P1------------->P3
        //
        //    The algorithm used here essentially computes the barycentric
        //    coordinates of the point P, and accepts it if both coordinates
        //    are between 0 and 1.  ( For a triangle, they must be positive,
        //    and sum to no more than 1.)  The same trick works for a parallelepiped.
        //
        //    05 August 2005: Thanks to Gernot Grabmair for pointing out that a previous
        //    version of this routine was incorrect.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], three vertices of the parallelogram.
        //    P1 should be directly connected to P2 and P3.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool PARALLELOGRAM_CONTAINS_POINT_2D is TRUE if P is inside the
        //    parallelogram, or on its boundary, and FALSE otherwise.
        //
    {
        int DIM_NUM = 2;

        double[] a = new double[DIM_NUM * (DIM_NUM + 1)];
        int info;
        bool value;
        //
        //  Set up the linear system
        //
        //    ( X2-X1  X3-X1 ) C1  = X-X1
        //    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
        //
        //  which is satisfied by the barycentric coordinates.
        //
        a[0 + 0 * DIM_NUM] = p2[0] - p1[0];
        a[1 + 0 * DIM_NUM] = p2[1] - p1[1];

        a[0 + 1 * DIM_NUM] = p3[0] - p1[0];
        a[1 + 1 * DIM_NUM] = p3[1] - p1[1];

        a[0 + 2 * DIM_NUM] = p[0] - p1[0];
        a[1 + 2 * DIM_NUM] = p[1] - p1[1];
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(DIM_NUM, 1, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("PARALLELOGRAM_CONTAINS_POINT_2D - Fatal error!");
            Console.WriteLine("  The linear system is singular.");
            Console.WriteLine("  The input data does not form a proper triangle.");
            return false;
        }

        switch (a[0 + 2 * DIM_NUM])
        {
            case < 0.0:
            case > 1.0:
                value = false;
                break;
            default:
            {
                switch (a[1 + 2 * DIM_NUM])
                {
                    case < 0.0:
                    case > 1.0:
                        value = false;
                        break;
                    default:
                        value = true;
                        break;
                }

                break;
            }
        }

        return value;
    }

    public static bool parallelogram_contains_point_3d(double[] p1, double[] p2, double[] p3,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_CONTAINS_POINT_3D determines if a point is inside a parallelogram in 3D.
        //
        //  Diagram:
        //
        //         P2..............
        //        /              .
        //       /              .
        //      /              .
        //    P1------------->P3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], three vertices of the parallelogram.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, int PARALLELOGRAM_CONTAINS_POINT_3D, is true if P is inside the
        //    parallelogram, or on its boundary, and false otherwise.
        //    A slight amount of leeway is allowed for error, since a three
        //    dimensional point may lie exactly in the plane of the parallelogram,
        //    and yet be computationally slightly outside it.
        //
    {
        int DIM_NUM = 3;
        double TOL = 0.00001;

        double dot;
        double dotb;
        double dott;
        double v;
        double[] p21 = new double[DIM_NUM];
        double[] p31 = new double[DIM_NUM];
        double[] pn12 = new double[DIM_NUM];
        double[] pn23;
        double[] pn31;
        //
        //  Compute V3, the vector normal to V1 = P2-P1 and V2 = P3-P1.
        //
        pn12[0] = (p2[1] - p1[1]) * (p3[2] - p1[2])
                  - (p2[2] - p1[2]) * (p3[1] - p1[1]);

        pn12[1] = (p2[2] - p1[2]) * (p3[0] - p1[0])
                  - (p2[0] - p1[0]) * (p3[2] - p1[2]);

        pn12[2] = (p2[0] - p1[0]) * (p3[1] - p1[1])
                  - (p2[1] - p1[1]) * (p3[0] - p1[0]);
        //
        //  If the component of V = P-P1 in the V3 direction is too large,
        //  then it does not lie in the parallelogram.
        //
        dot = (p[0] - p1[0]) * pn12[0]
              + (p[1] - p1[1]) * pn12[1]
              + (p[2] - p1[2]) * pn12[2];

        v = Math.Sqrt(Math.Pow(p2[0] - p[0], 2)
                      + Math.Pow(p2[1] - p[1], 2)
                      + Math.Pow(p2[2] - p[2], 2));

        if (TOL * (1.0 + v) < Math.Abs(dot))
        {
            return false;
        }

        //
        //  Compute V23, the vector normal to V2 and V3.
        //
        p31[0] = p3[0] - p1[0];
        p31[1] = p3[1] - p1[1];
        p31[2] = p3[2] - p1[2];

        pn23 = typeMethods.r8vec_cross_product_3d(p31, pn12);
        //
        //  Compute ALPHA = ( V dot V23 ) / ( V1 dot V23 )
        //
        dott = (p[0] - p1[0]) * pn23[0]
               + (p[1] - p1[1]) * pn23[1]
               + (p[2] - p1[2]) * pn23[2];

        dotb =
            (p2[0] - p1[0]) * pn23[0] +
            (p2[1] - p1[1]) * pn23[1] +
            (p2[2] - p1[2]) * pn23[2];

        switch (dotb)
        {
            case < 0.0:
                dott = -dott;
                dotb = -dotb;
                break;
        }

        if (dott < 0.0 || dotb < dott)
        {
            return false;
        }

        //
        //  Compute V31, the vector normal to V3 and V1.
        //
        p21[0] = p2[0] - p1[0];
        p21[1] = p2[1] - p1[1];
        p21[2] = p2[2] - p1[2];

        pn31 = typeMethods.r8vec_cross_product_3d(pn12, p21);
        //
        //  Compute BETA = ( V dot V31 ) / ( V2 dot V31 )
        //
        dott = (p[0] - p1[0]) * pn31[0]
               + (p[1] - p1[1]) * pn31[1]
               + (p[2] - p1[2]) * pn31[2];

        dotb =
            (p3[0] - p1[0]) * pn31[0] +
            (p3[1] - p1[1]) * pn31[1] +
            (p3[2] - p1[2]) * pn31[2];

        switch (dotb)
        {
            case < 0.0:
                dott = -dott;
                dotb = -dotb;
                break;
        }

        if (dott < 0.0 || dotb < dott)
        {
            return false;
        }

        //
        //  V = ALPHA * V1 + BETA * V2, where both ALPHA and BETA are between
        //  0 and 1.
        //
        return true;
    }

    public static double parallelogram_point_dist_3d(double[] p1, double[] p2, double[] p3,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_POINT_DIST_3D: distance ( parallelogram, point ) in 3D.
        //
        //  Diagram:
        //
        //         P2..............
        //        /              .
        //       /              .
        //      /              .
        //    P1------------->P3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], P3[3], three vertices of the parallelogram.
        //
        //    Input, double P[3], the point to be checked.
        //
        //    Output, double PARALLELOGRAM_POINT_DIST_3D, the distance from the point
        //    to the parallelogram.  DIST is zero if the point lies exactly on the
        //    parallelogram.
        //
    {
        int DIM_NUM = 3;

        double dis13;
        double dis21;
        double dis34;
        double dis42;
        double dist;
        bool inside;
        double t;
        double temp;
        double[] p4 = new double[DIM_NUM];
        double[] pn = new double[DIM_NUM];
        double[] pp = new double[DIM_NUM];
        //
        //  Compute P, the unit normal to P2-P1 and P3-P1:
        //
        pp[0] = (p2[1] - p1[1]) * (p3[2] - p1[2])
                - (p2[2] - p1[2]) * (p3[1] - p1[1]);

        pp[1] = (p2[2] - p1[2]) * (p3[0] - p1[0])
                - (p2[0] - p1[0]) * (p3[2] - p1[2]);

        pp[2] = (p2[0] - p1[0]) * (p3[1] - p1[1])
                - (p2[1] - p1[1]) * (p3[0] - p1[0]);

        temp = Math.Sqrt(pp[0] * pp[0] + pp[1] * pp[1] + pp[2] * pp[2]);

        switch (temp)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("PARALLELOGRAM_POINT_DIST_3D - Fatal error!");
                Console.WriteLine("  The normal vector is zero.");
                return 1;
        }

        pp[0] /= temp;
        pp[1] /= temp;
        pp[2] /= temp;
        //
        //  Find PN, the nearest point to P in the plane.
        //
        t = pp[0] * (p[0] - p1[0])
            + pp[1] * (p[1] - p1[1])
            + pp[2] * (p[2] - p1[2]);

        pn[0] = p[0] - pp[0] * t;
        pn[1] = p[1] - pp[1] * t;
        pn[2] = p[2] - pp[2] * t;
        //
        //  if PN lies WITHIN the parallelogram, we're done.
        //
        inside = parallelogram_contains_point_3d(p1, p2, p3, p);

        switch (inside)
        {
            case true:
                dist = Math.Sqrt(Math.Pow(pn[0] - p[0], 2)
                                 + Math.Pow(pn[1] - p[1], 2)
                                 + Math.Pow(pn[2] - p[2], 2));
                return dist;
        }

        //
        //  Otherwise, find the distance between P and each of the
        //  four line segments that make up the boundary of the parallelogram.
        //
        p4[0] = p2[0] + p3[0] - p1[0];
        p4[1] = p2[1] + p3[1] - p1[1];
        p4[2] = p2[2] + p3[2] - p1[2];

        dis13 = Segments.segment_point_dist_3d(p1, p3, p);

        dist = dis13;

        dis34 = Segments.segment_point_dist_3d(p3, p4, p);

        if (dis34 < dist)
        {
            dist = dis34;
        }

        dis42 = Segments.segment_point_dist_3d(p4, p2, p);

        if (dis42 < dist)
        {
            dist = dis42;
        }

        dis21 = Segments.segment_point_dist_3d(p2, p1, p);

        if (dis21 < dist)
        {
            dist = dis21;
        }

        return dist;
    }

}