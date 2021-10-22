using System;
using Burkardt.Geometry;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Quadrilateral
{
    public static class Geometry
    {
        public static double quad_area_2d(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_AREA_2D computes the area of a quadrilateral in 2D.
            //
            //  Discussion:
            //
            //    This algorithm should be able to handle nonconvex quadrilaterals.
            //
            //    3----2
            //    |   /|
            //    |  / |    We subdivide the quadrilateral into triangles (0,1,2)
            //    | /  |    and (2,3,0), computer their areas, and add.
            //    |/   |
            //    0----1
            //
            //    Thanks to Eduardo Olmedo of Universidad Politecnica de Madrid for
            //    pointing out an error in a previous version of this routine!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 December 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[2*4], the vertices of the quadrilateral,
            //    in counter clockwise order.
            //
            //    Output, double QUAD_AREA_2D, the area of the quadrilateral.
            //
        {
            int DIM_NUM = 2;

            double area;
            double[] t = new double[DIM_NUM * 3];

            area = 0.0;

            t[0 + 0 * 2] = q[0 + 0 * 2];
            t[1 + 0 * 2] = q[1 + 0 * 2];
            t[0 + 1 * 2] = q[0 + 1 * 2];
            t[1 + 1 * 2] = q[1 + 1 * 2];
            t[0 + 2 * 2] = q[0 + 2 * 2];
            t[1 + 2 * 2] = q[1 + 2 * 2];

            area = area + typeMethods.triangle_area_2d(t);

            t[0 + 0 * 2] = q[0 + 2 * 2];
            t[1 + 0 * 2] = q[1 + 2 * 2];
            t[0 + 1 * 2] = q[0 + 3 * 2];
            t[1 + 1 * 2] = q[1 + 3 * 2];
            t[0 + 2 * 2] = q[0 + 0 * 2];
            t[1 + 2 * 2] = q[1 + 0 * 2];

            area = area + typeMethods.triangle_area_2d(t);

            return area;
        }

        public static double quad_area2_2d(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_AREA2_2D computes the area of a quadrilateral in 2D.
            //
            //  Discussion:
            //
            //    A quadrilateral is a polygon defined by 4 vertices.
            //
            //    This algorithm computes the area of the related
            //    Varignon parallelogram first.
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
            //    Input, double Q[2*4], the vertices, specified in
            //    counter clockwise order.
            //
            //    Output, double QUAD_AREA2_2D, the area of the quadrilateral.
            //
        {
            double area;
            int i;
            int j;
            double[] p;

            p = new double[2 * 4];
            //
            //  Define a parallelogram by averaging consecutive vertices.
            //
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    p[i + j * 2] = (q[i + j * 2] + q[i + (j + 1) * 2]) / 2.0;
                }
            }

            for (i = 0; i < 2; i++)
            {
                p[i + 3 * 2] = (q[i + 3 * 2] + q[i + 0 * 2]) / 2.0;
            }

            //
            //  Compute the area.
            //
            area = Burkardt.Parallelogram.Geometry.parallelogram_area_2d(p);
            //
            //  The quadrilateral's area is twice that of the parallelogram.
            //
            area = 2.0 * area;

            return area;
        }

        public static double quad_area_3d(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_AREA_3D computes the area of a quadrilateral in 3D.
            //
            //  Discussion:
            //
            //    A quadrilateral is a polygon defined by 4 vertices.
            //
            //    It is assumed that the four vertices of the quadrilateral
            //    are coplanar.
            //
            //    This algorithm computes the area of the related
            //    Varignon parallelogram first.
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
            //    Input, double Q[3*4], the vertices, specified in
            //    counter clockwise order.
            //
            //    Output, double QUAD_AREA_3D, the area of the quadrilateral.
            //
        {
            double area;
            int i;
            int j;
            double[] p;

            p = new double[3 * 4];
            //
            //  Define a parallelogram by averaging consecutive vertices.
            //
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    p[i + j * 3] = (q[i + j * 3] + q[i + (j + 1) * 3]) / 2.0;
                }
            }

            for (i = 0; i < 3; i++)
            {
                p[i + 3 * 3] = (q[i + 3 * 3] + q[i + 0 * 3]) / 2.0;
            }

            //
            //  Compute the area.
            //
            area = Burkardt.Parallelogram.Geometry.parallelogram_area_3d(p);
            //
            //  The quadrilateral's area is twice that of the parallelogram.
            //
            area = 2.0 * area;

            return area;
        }

        public static bool quad_contains_point_2d(double[] q, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_CONTAINS_POINT_2D finds if a point is inside a convex quadrilateral in 2D.
            //
            //  Discussion:
            //
            //    This method will only handle convex quadrilaterals.
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
            //    Input, double Q[2*4], the vertices of the quadrilateral, in counter clockwise order.
            //
            //    Input, double P[2], the point to be checked.
            //
            //    Output, bool QUAD_CONTAINS_POINT, is TRUE if the point is inside
            //    the quadrilateral or on its boundary, and FALSE otherwise.
            //
        {
            int DIM_NUM = 2;

            if (Angle.anglei_rad_2d(q, q, q, +0 * 2, +1 * 2, +2 * 2) <
                Angle.anglei_rad_2d(q, q, p, +0 * 2, +1 * 2, 0))
            {
                return false;
            }

            if (Angle.anglei_rad_2d(q, q, q, +1 * 2, +2 * 2, +3 * 2) <
                Angle.anglei_rad_2d(q, q, p, +1 * 2, +2 * 2, 0))
            {
                return false;
            }

            if (Angle.anglei_rad_2d(q, q, q, +2 * 2, +3 * 2, +0 * 2) <
                Angle.anglei_rad_2d(q, q, p, +2 * 2, +3 * 2, 0))
            {
                return false;
            }

            if (Angle.anglei_rad_2d(q, q, q, +3 * 2, +0 * 2, +1 * 2) <
                Angle.anglei_rad_2d(q, q, p, +3 * 2, +0 * 2, 0))
            {
                return false;
            }

            return true;
        }

        public static void quad_convex_random(ref int seed, ref double[] xy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
            //
            //  Description:
            //
            //    The quadrilateral is constrained in that the vertices must all lie
            //    with the unit square.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number
            //    generator.
            //
            //    Output, double XY[2*NODE_NUM], the coordinates of the
            //    nodes of the quadrilateral, given in counterclockwise order.
            //
        {
            int[] hull = new int[4];
            int hull_num = 0;
            int i;
            int j;
            double[] xy_random;

            for (;;)
            {
                //
                //  Generate 4 random points.
                //
                xy_random = UniformRNG.r8mat_uniform_01_new(2, 4, ref seed);
                //
                //  Determine the convex hull.
                //
                Burkardt.PointsNS.Geometry.points_hull_2d(4, xy_random, ref hull_num, ref hull);
                //
                //  If HULL_NUM < 4, then our convex hull is a triangle.
                //  Try again.
                //
                if (hull_num == 4)
                {
                    break;
                }
            }

            //
            //  Make an ordered copy of the random points.
            //
            for (j = 0; j < 4; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    xy[i + j * 2] = xy_random[i + (hull[j] - 1) * 2];
                }
            }
        }

        public static double quad_point_dist_2d(double[] q, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_POINT_DIST_2D finds the distance from a point to a quadrilateral in 2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[2*4], the vertices of the quadrilateral.
            //
            //    Input, double P[2], the point to be checked.
            //
            //    Output, double QUAD_POINT_DIST_2D, the distance from the point to the quadrilateral.
            //    DIST is zero if the point lies exactly on the quadrilateral.
            //
        {
            int DIM_NUM = 2;

            double dist;
            double dist2;
            int j;
            int jp1;
            int side_num = 4;
            //
            //  Find the distance to each of the line segments.
            //
            dist = typeMethods.r8_huge();

            for (j = 0; j < side_num; j++)
            {
                jp1 = typeMethods.i4_wrap(j + 1, 0, side_num - 1);

                dist2 = Segments.segment_point_dist_2d(q, q, p, + j * 2, + jp1 * 2 );

                if (dist2 < dist)
                {
                    dist = dist2;
                }
            }

            return dist;
        }

        public static double quad_point_dist_signed_2d(double[] q, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_POINT_DIST_SIGNED_2D: signed distanct ( quadrilateral, point ) in 2D.
            //
            //  Discussion:
            //
            //    The quadrilateral must be convex.  DIST_SIGNED is actually the maximum
            //    of the signed distances from the point to each of the four lines that
            //    make up the quadrilateral.
            //
            //    Essentially, if the point is outside the convex quadrilateral,
            //    only one of the signed distances can be positive, or two can
            //    be positive and equal.
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
            //    Input, double Q[2*4], the vertices of the quadrilateral.
            //
            //    Input, double P[2], the point which is to be checked.
            //
            //    Output, double QUAD_POINT_DIST_SIGNED_2D, the signed distance from
            //    the point to the convex quadrilateral.  If the distance is
            //    0.0, the point is on the boundary;
            //    negative, the point is in the interior;
            //    positive, the point is in the exterior.
            //
        {
            int DIM_NUM = 2;

            double dis1;
            double dis2;
            double dist_signed;
            double[] pm = new double[DIM_NUM];
            //
            //  Compare the signed distance from each line segment to the point,
            //  with the signed distance to the midpoint of the opposite line.
            //
            //  The signed distances should all be negative if the point is inside.
            //
            //  Side 12
            //
            dis1 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, p, +0 * 2, +1 * 2);

            pm[0] = 0.5 * (q[0 + 2 * 2] + q[0 + 3 * 2]);
            pm[1] = 0.5 * (q[1 + 2 * 2] + q[1 + 3 * 2]);

            dis2 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, pm, +0 * 2, +1 * 2);

            if (0.0 < dis2)
            {
                dis1 = -dis1;
            }

            dist_signed = dis1;
            //
            //  Side 23
            //
            dis1 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, p, +1 * 2, +2 * 2);

            pm[0] = 0.5 * (q[0 + 3 * 2] + q[0 + 0 * 2]);
            pm[1] = 0.5 * (q[1 + 3 * 2] + q[1 + 0 * 2]);

            dis2 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, pm, +1 * 2, +2 * 2);

            if (0.0 < dis2)
            {
                dis1 = -dis1;
            }

            dist_signed = Math.Max(dist_signed, dis1);
            //
            //  Side 34
            //
            dis1 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, p, +2 * 2, +3 * 2);

            pm[0] = 0.5 * (q[0 + 0 * 2] + q[0 + 1 * 2]);
            pm[1] = 0.5 * (q[1 + 0 * 2] + q[1 + 1 * 2]);

            dis2 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, pm, +2 * 2, +3 * 2);

            if (0.0 < dis2)
            {
                dis1 = -dis1;
            }

            dist_signed = Math.Max(dist_signed, dis1);
            //
            //  Side 41
            //
            dis1 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, p, +3 * 2, +0 * 2);

            pm[0] = 0.5 * (q[0 + 1 * 2] + q[0 + 2 * 2]);
            pm[1] = 0.5 * (q[1 + 1 * 2] + q[1 + 2 * 2]);

            dis2 = Burkardt.LineNS.Geometry.line_exp_point_dist_signed_2d(q, q, pm, +3 * 2, +0 * 2);

            if (0.0 < dis2)
            {
                dis1 = -dis1;
            }

            dist_signed = Math.Max(dist_signed, dis1);

            return dist_signed;
        }

        public static void quad_point_near_2d(double[] q, double[] p, ref double[] pn,
                ref double dist)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    QUAD_POINT_NEAR_2D computes the nearest point on a quadrilateral in 2D.
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
            //    Input, double Q[2*4], the quadrilateral vertices.
            //
            //    Input, double P[2], the point whose nearest quadrilateral point
            //    is to be determined.
            //
            //    Output, double PN[2], the nearest point to P.
            //
            //    Output, double *DIST, the distance from the point to the
            //    quadrilateral.
            //
        {
            int DIM_NUM = 2;

            double dist2 = 0;
            int j;
            int jp1;
            double[] pn2 = new double[DIM_NUM];
            int side_num = 4;
            double tval = 0;

            dist = typeMethods.r8_huge();
            typeMethods.r8vec_zero(DIM_NUM, ref pn);

            for (j = 0; j < side_num; j++)
            {
                jp1 = typeMethods.i4_wrap(j + 1, 0, side_num - 1);

                Segments.segment_point_near_2d(q, q, p, ref pn2, ref dist2, ref tval,  + j * 2,  + jp1 * 2);

                if (dist2 < dist)
                {
                    dist = dist2;
                    typeMethods.r8vec_copy(DIM_NUM, pn2, ref pn);
                }
            }
        }

    }
}