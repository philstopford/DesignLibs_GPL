using Burkardt.Types;

namespace Burkardt.Geometry;

public static class Half
{

    public static bool halfplane_contains_point_2d(double[] pa, double[] pb, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALFPLANE_CONTAINS_POINT_2D reports if a half-plane contains a point in 2d.
        //
        //  Discussion:
        //
        //    The halfplane is assumed to be all the points "to the left" of the
        //    line segment from PA = (XA,YA) to PB = (XB,YB).  Thus, one way to
        //    understand where the point P  is, is to compute the signed
        //    area of the triangle ( PA, PB, P ).
        //
        //    If this area is
        //      positive, the point is strictly inside the halfplane,
        //      zero, the point is on the boundary of the halfplane,
        //      negative, the point is strictly outside the halfplane.
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
        //    Input, double PA[2], PB[2], two points on the line defining the half plane.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool HALFPLANE_CONTAINS_POINT_2D, is TRUE if the halfplane
        //    contains the point, and FALSE otherwise.
        //
    {
        double area_signed;

        area_signed = 0.5 *
                      (pa[0] * (pb[1] - p[1])
                       + pb[0] * (p[1] - pa[1])
                       + p[0] * (pa[1] - pb[1]));

        return 0.0 <= area_signed;
    }

    public static int halfspace_imp_triangle_int_3d(double a, double b, double c, double d,
            double[] t, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALFSPACE_IMP_TRIANGLE_INT_3D: intersection ( implicit halfspace, triangle ) in 3D.
        //
        //  Discussion:
        //
        //    The implicit form of a half-space in 3D may be described as the set
        //    of points on or "above" an implicit plane:
        //
        //      0 <= A * X + B * Y + C * Z + D
        //
        //    The triangle is specified by listing its three vertices.
        //
        //    The intersection may be described by the number of vertices of the
        //    triangle that are included in the halfspace, and by the location of
        //    points between vertices that separate a side of the triangle into
        //    an included part and an unincluded part.
        //
        //    0 vertices, 0 separators    (no intersection)
        //    1 vertex,   0 separators    (point intersection)
        //    2 vertices, 0 separators    (line intersection)
        //    3 vertices, 0 separators    (triangle intersection)
        //
        //    1 vertex,   2 separators,     (intersection is a triangle)
        //    2 vertices, 2 separators,   (intersection is a quadrilateral).
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
        //    Input, double A, B, C, D, the parameters that define the implicit plane,
        //    which in turn define the implicit halfspace.
        //
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double P[3*4], the coordinates of the
        //    intersection points.  The points will lie in sequence on the triangle.
        //    Some points will be vertices, and some may be separators.
        //
        //    Output, int HALFSPACE_IMP_TRIANGLE_INT_3D, the number of intersection
        //    points returned, which will always be between 0 and 4.
        //
    {
        double dist1;
        double dist2;
        double dist3;
        int int_num;
        //
        //  Compute the signed distances between the vertices and the plane.
        //
        dist1 = a * t[0 + 0 * 3] + b * t[1 + 0 * 3] + c * t[2 + 0 * 3] + d;
        dist2 = a * t[0 + 1 * 3] + b * t[1 + 1 * 3] + c * t[2 + 1 * 3] + d;
        dist3 = a * t[0 + 2 * 3] + b * t[1 + 2 * 3] + c * t[2 + 2 * 3] + d;
        //
        //  Now we can find the intersections.
        //
        int_num = halfspace_triangle_int_3d(dist1, dist2, dist3, t, ref p);

        return int_num;
    }

    public static int halfspace_norm_triangle_int_3d(double[] pp, double[] pn, double[] t,
            ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALFSPACE_NORM_TRIANGLE_INT_3D: intersection ( normal halfspace, triangle ) in 3D.
        //
        //  Discussion:
        //
        //    The normal form of a halfspace in 3D may be described as the set
        //    of points P on or "above" a plane described in normal form:
        //
        //      PP is a point on the plane,
        //      PN is the unit normal vector, pointing "out" of the halfspace
        //
        //    The triangle is specified by listing its three vertices.
        //
        //    The intersection may be described by the number of vertices of the
        //    triangle that are included in the halfspace, and by the location of
        //    points between vertices that separate a side of the triangle into
        //    an included part and an unincluded part.
        //
        //    0 vertices, 0 separators    (no intersection)
        //    1 vertex, 0 separators  (point intersection)
        //    2 vertices, 0 separators    (line intersection)
        //    3 vertices, 0 separators    (triangle intersection)
        //
        //    1 vertex, 2 separators,     (intersection is a triangle)
        //    2 vertices, 2 separators,   (intersection is a quadrilateral).
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
        //    Input, double PP[3], a point on the bounding plane that defines
        //    the halfspace.
        //
        //    Input, double PN[3], the components of the normal vector to the
        //    bounding plane that defines the halfspace.  By convention, the
        //    normal vector points "outwards" from the halfspace.
        //
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double P[3*4], the intersection points.  The points will lie
        //    in sequence on the triangle.  Some points will be vertices, and some
        //    may be separators.
        //
        //    Output, int HALFSPACE_NORM_TRIANGLE_INT_3D, the number of intersection
        //    points returned, which will always be between 0 and 4.
        //
    {
        double dist1;
        double dist2;
        double dist3;
        int int_num;
        //
        //  Compute the signed distances between the vertices and the plane.
        //
        dist1 = typeMethods.r8vec_dot_product(3, pn, t, a2Index: +0 * 3);
        dist2 = typeMethods.r8vec_dot_product(3, pn, t, a2Index: +1 * 3);
        dist3 = typeMethods.r8vec_dot_product(3, pn, t, a2Index: +2 * 3);
        //
        //  Now we can find the intersections.
        //
        int_num = halfspace_triangle_int_3d(dist1, dist2, dist3, t, ref p);

        return int_num;
    }

    public static int halfspace_triangle_int_3d(double dist1, double dist2, double dist3,
            double[] t, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALFSPACE_TRIANGLE_INT_3D: intersection ( halfspace, triangle ) in 3D.
        //
        //  Discussion:
        //
        //    The triangle is specified by listing its three vertices.
        //
        //    The halfspace is not described in the input data.  Rather, the
        //    distances from the triangle vertices to the halfspace are given.
        //
        //    The intersection may be described by the number of vertices of the
        //    triangle that are included in the halfspace, and by the location of
        //    points between vertices that separate a side of the triangle into
        //    an included part and an unincluded part.
        //
        //    0 vertices, 0 separators    (no intersection)
        //    1 vertex, 0 separators  (point intersection)
        //    2 vertices, 0 separators    (line intersection)
        //    3 vertices, 0 separators    (triangle intersection)
        //
        //    1 vertex, 2 separators,     (intersection is a triangle)
        //    2 vertices, 2 separators,   (intersection is a quadrilateral).
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
        //    Input, double DIST1, DIST2, DIST3, the distances from each of the
        //    three vertices of the triangle to the halfspace.  The distance is
        //    zero if a vertex lies within the halfspace, or on the plane that
        //    defines the boundary of the halfspace.  Otherwise, it is the
        //    distance from that vertex to the bounding plane.
        //
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double P[3*4], the coordinates of the
        //    intersection points.  The points will lie in sequence on the triangle.
        //    Some points will be vertices, and some may be separators.
        //
        //    Output, int HALFSPACE_TRIANGLE_INT_3D, the number of intersection points
        //    returned, which will always be between 0 and 4.
        //
    {
        int int_num;
        //
        //  Walk around the triangle, looking for vertices that are included,
        //  and points of separation.
        //
        int_num = 0;

        switch (dist1)
        {
            case <= 0.0:
                p[0 + int_num * 3] = t[0 + 0 * 3];
                p[1 + int_num * 3] = t[1 + 0 * 3];
                p[2 + int_num * 3] = t[2 + 0 * 3];
                int_num += 1;
                break;
        }

        switch (dist1 * dist2)
        {
            case < 0.0:
                p[0 + int_num * 3] = (dist1 * t[0 + 1 * 3] - dist2 * t[0 + 0 * 3]) / (dist1 - dist2);
                p[1 + int_num * 3] = (dist1 * t[1 + 1 * 3] - dist2 * t[1 + 0 * 3]) / (dist1 - dist2);
                p[2 + int_num * 3] = (dist1 * t[2 + 1 * 3] - dist2 * t[2 + 0 * 3]) / (dist1 - dist2);
                int_num += 1;
                break;
        }

        switch (dist2)
        {
            case <= 0.0:
                p[0 + int_num * 3] = t[0 + 1 * 3];
                p[1 + int_num * 3] = t[1 + 1 * 3];
                p[2 + int_num * 3] = t[2 + 1 * 3];
                int_num += 1;
                break;
        }

        switch (dist2 * dist3)
        {
            case < 0.0:
                p[0 + int_num * 3] = (dist2 * t[0 + 2 * 3] - dist3 * t[0 + 1 * 3]) / (dist2 - dist3);
                p[1 + int_num * 3] = (dist2 * t[1 + 2 * 3] - dist3 * t[1 + 1 * 3]) / (dist2 - dist3);
                p[2 + int_num * 3] = (dist2 * t[2 + 2 * 3] - dist3 * t[2 + 0 * 3]) / (dist2 - dist3);
                int_num += 1;
                break;
        }

        switch (dist3)
        {
            case <= 0.0:
                p[0 + int_num * 3] = t[0 + 2 * 3];
                p[1 + int_num * 3] = t[1 + 2 * 3];
                p[2 + int_num * 3] = t[2 + 2 * 3];
                int_num += 1;
                break;
        }

        switch (dist3 * dist1)
        {
            case < 0.0:
                p[0 + int_num * 3] = (dist3 * t[0 + 0 * 3] - dist1 * t[0 + 2 * 3]) / (dist3 - dist1);
                p[1 + int_num * 3] = (dist3 * t[1 + 0 * 3] - dist1 * t[1 + 2 * 3]) / (dist3 - dist1);
                p[2 + int_num * 3] = (dist3 * t[2 + 0 * 3] - dist1 * t[2 + 2 * 3]) / (dist3 - dist1);
                int_num += 1;
                break;
        }

        return int_num;
    }

}