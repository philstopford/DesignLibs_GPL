using System;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.Cube;

public static class Geometry
{
    public static bool box_contains_point_2d(double[] p1, double[] p2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_CONTAINS_POINT_2D reports if a point is contained in a box in 2D.
        //
        //  Discussion:
        //
        //    A box in 2D is a rectangle with sides aligned on coordinate
        //    axes.  It can be described by its low and high corners, P1 and P2
        //    as the set of points P satisfying:
        //
        //      P1(1:2) <= P(1:2) <= P2(1:2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the minimum and maximum X and Y
        //    values, which define the box.
        //
        //    Input, double P[2], the coordinates of the point.
        //
        //    Output, bool BOX_CONTAINS_POINT_2D, is TRUE if the box contains
        //    the point.
        //
    {
        return p1[0] <= p[0] && p[0] <= p2[0] && p1[1] <= p[1] && p[1] <= p2[1];
    }

    public static bool box_contains_point_nd(int dim_num, double[] p1, double[] p2, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_CONTAINS_POINT_ND reports if a point is contained in a box in ND.
        //
        //  Discussion:
        //
        //    A box in ND is a rectangle with sides aligned on coordinate
        //    axes.  It can be described by its low and high corners, P1 and P2
        //    as the set of points P satisfying:
        //
        //      P1(1:DIM_NUM) <= P(1:DIM_NUM) <= P2(1:DIM_NUM).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double P1[DIM_NUM], P2[DIM_NUM], the minimum and maximum X and Y
        //    values, which define the box.
        //
        //    Input, double P[DIM_NUM], the coordinates of the point.
        //
        //    Output, bool BOX_CONTAINS_POINT_ND, is TRUE if the box contains
        //    the point.
        //
    {
        int i;

        for (i = 0; i < dim_num; i++)
        {
            if (p[i] < p1[i] || p2[i] < p[i])
            {
                return false;
            }
        }

        return true;
    }

    public static void box_ray_int_2d(double[] p1, double[] p2, double[] pa,
            double[] pb, double[] pint)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_RAY_INT_2D: intersection ( box, ray ) in 2D.
        //
        //  Discussion:
        //
        //    A box in 2D is a rectangle with sides aligned on coordinate
        //    axes.  It can be described by its low and high corners, P1 and P2
        //    as the set of points P satisfying:
        //
        //      P1(1:2) <= P(1:2) <= P2(1:2).
        //
        //    The origin of the ray is assumed to be inside the box.  This
        //    guarantees that the ray will intersect the box in exactly one point.
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
        //    Input, double P1[2], the lower left corner of the box.
        //
        //    Input, double P2[2], the upper right corner of the box.
        //
        //    Input, double PA[2], the origin of the ray, which should be
        //    inside the box.
        //
        //    Input, double PB[2], a second point on the ray.
        //
        //    Output, double PINT[2], the point on the box intersected by the ray.
        //
    {
        const int DIM_NUM = 2;

        int ival = 0;
        double[] pc = new double[DIM_NUM];
        double[] pd = new double[DIM_NUM];
        int side;

        for (side = 1; side <= 4; side++)
        {
            switch (side)
            {
                case 1:
                    pc[0] = p1[0];
                    pc[1] = p1[1];
                    pd[0] = p2[0];
                    pd[1] = p1[1];
                    break;
                case 2:
                    pc[0] = p2[0];
                    pc[1] = p1[1];
                    pd[0] = p2[0];
                    pd[1] = p2[1];
                    break;
                case 3:
                    pc[0] = p2[0];
                    pc[1] = p2[1];
                    pd[0] = p1[0];
                    pd[1] = p2[1];
                    break;
                case 4:
                    pc[0] = p1[0];
                    pc[1] = p2[1];
                    pd[0] = p1[0];
                    pd[1] = p1[1];
                    break;
            }

            bool inside = Angle.angle_contains_ray_2d(pc, pa, pd, pb);

            if (inside)
            {
                break;
            }

            switch (side)
            {
                case 4:
                    Console.WriteLine("");
                    Console.WriteLine("BOX_RAY_INT_2D - Fatal error!");
                    Console.WriteLine("  No intersection could be found.");
                    return;
            }

        }

        LineNS.Geometry.lines_exp_int_2d(pa, pb, pc, pd, ref ival, ref pint);
    }

    public static int box_segment_clip_2d(double[] p1, double[] p2, double[] pa,
            double[] pb)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_SEGMENT_CLIP_2D uses a box to clip a line segment in 2D.
        //
        //  Discussion:
        //
        //    A box in 2D is a rectangle with sides aligned on coordinate
        //    axes.  It can be described by its low and high corners, P1 and P2
        //    as the set of points P satisfying:
        //
        //      P1(1:2) <= P(1:2) <= P2(1:2).
        //
        //    A line segment is the finite portion of a line that lies between
        //    two points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], the minimum and maximum X and Y
        //    values, which define the box.
        //
        //    Input/output, double PA[2], PB[2]; on input, the endpoints
        //    of a line segment.  On output, the endpoints of the portion of the
        //    line segment that lies inside the box.  However, if no part of the
        //    initial line segment lies inside the box, the output value is the
        //    same as the input value.
        //
        //    Output, int BOX_SEGMENT_CLIP_LINE_2D:
        //    -1, no part of the line segment is within the box.
        //     0, no clipping was necessary.
        //     1, P1 was clipped.
        //     2, P2 was clipped.
        //     3, P1 and P2 were clipped.
        //
    {
        int ival;
        double x;
        double y;

        bool clip_a = false;
        bool clip_b = false;
        //
        //  Require that XMIN <= X.
        //
        if (pa[0] < p1[0] && pb[0] < p1[0])
        {
            ival = -1;
            return ival;
        }

        if (pa[0] < p1[0] && p1[0] <= pb[0])
        {
            x = p1[0];
            y = pa[1] + (pb[1] - pa[1]) * (x - pa[0]) / (pb[0] - pa[0]);
            pa[0] = x;
            pa[1] = y;
            clip_a = true;
        }
        else if (p1[0] <= pa[0] && pb[0] < p1[0])
        {
            x = p1[0];
            y = pa[1] + (pb[1] - pa[1]) * (x - pa[0]) / (pb[0] - pa[0]);
            pb[0] = x;
            pb[1] = y;
            clip_b = true;
        }

        //
        //  Require that X <= XMAX.
        //
        if (p2[0] < pa[0] && p2[0] < pb[0])
        {
            ival = -1;
            return ival;
        }

        if (p2[0] < pa[0] && pb[0] <= p2[0])
        {
            x = p2[0];
            y = pa[1] + (pb[1] - pa[1]) * (x - pa[0]) / (pb[0] - pa[0]);
            pa[0] = x;
            pa[1] = y;
            clip_a = true;
        }
        else if (pa[0] <= p2[0] && p2[0] < pb[0])
        {
            x = p2[0];
            y = pa[1] + (pb[1] - pa[1]) * (x - pa[0]) / (pb[0] - pa[0]);
            pb[0] = x;
            pb[1] = y;
            clip_b = true;
        }

        //
        //  Require that YMIN <= Y.
        //
        if (pa[1] < p1[1] && pb[1] < p1[1])
        {
            ival = -1;
            return ival;
        }

        if (pa[1] < p1[1] && p1[1] <= pb[1])
        {
            y = p1[1];
            x = pa[0] + (pb[0] - pa[0]) * (y - pa[1]) / (pb[1] - pa[1]);

            pa[0] = x;
            pa[1] = y;
            clip_a = true;
        }
        else if (p1[1] <= pa[1] && pb[1] < p1[1])
        {
            y = p1[1];
            x = pa[0] + (pb[0] - pa[0]) * (y - pa[1]) / (pb[1] - pa[1]);
            pb[0] = x;
            pb[1] = y;
            clip_b = true;
        }

        //
        //  Require that Y <= YMAX.
        //
        if (p2[1] < pa[1] && p2[1] < pb[1])
        {
            ival = -1;
            return ival;
        }

        if (p2[1] < pa[1] && pb[1] <= p2[1])
        {
            y = p2[1];
            x = pa[0] + (pb[0] - pa[0]) * (y - pa[1]) / (pb[1] - pa[1]);
            pa[0] = x;
            pa[1] = y;
            clip_a = true;
        }
        else if (pa[1] <= p2[1] && p2[1] < pb[1])
        {
            y = p2[1];
            x = pa[0] + (pb[0] - pa[0]) * (y - pa[1]) / (pb[1] - pa[1]);
            pb[0] = x;
            pb[1] = y;
            clip_b = true;
        }

        ival = 0;

        switch (clip_a)
        {
            case true:
                ival += 1;
                break;
        }

        switch (clip_b)
        {
            case true:
                ival += 2;
                break;
        }

        return ival;
    }

    public static bool box01_contains_point_2d(double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX01_CONTAINS_POINT_2D reports if a point is contained in the unit box in 2D.
        //
        //  Discussion:
        //
        //    A unit box in 2D is a rectangle with sides aligned on coordinate
        //    axes.  It can be described as the set of points P satisfying:
        //
        //      0 <= P(1:2) <= 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P[2], the coordinates of the point.
        //
        //    Output, bool BOX01_CONTAINS_POINT_2D, is TRUE if the box contains
        //    the point.
        //
    {
        return p[0] switch
        {
            >= 0.0 and <= 1.0 when 0.0 <= p[1] && p[1] <= 1.0 => true,
            _ => false
        };
    }

    public static bool box01_contains_point_nd(int dim_num, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX01_CONTAINS_POINT_ND reports if a point is contained in the unit box in ND.
        //
        //  Discussion:
        //
        //    A unit box is assumed to be a rectangle with sides aligned on coordinate
        //    axes.  It can be described as the set of points P satisfying:
        //
        //      0.0 <= P(1:DIM_NUM) <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double P[DIM_NUM], the coordinates of the point.
        //
        //    Output, bool BOX01_CONTAINS_POINT_ND, is TRUE if the box contains
        //    the point.
        //
    {
        int i;

        for (i = 0; i < dim_num; i++)
        {
            switch (p[i])
            {
                case < 0.0:
                case > 1.0:
                    return false;
            }
        }

        return true;
    }

    public static void cube_shape_3d(int point_num, int face_num, int face_order_max,
            ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_SHAPE_3D describes a cube in 3D.
        //
        //  Discussion:
        //
        //    The vertices lie on the unit sphere.
        //
        //    The dual of the octahedron is the cube.
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
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int FACE_NUM, the number of faces.
        //
        //    Input, int FACE_ORDER_MAX, the maximum number of vertices
        //    per face.
        //
        //    Output, double POINT_COORD[3*POINT_NUM], the point coordinates.
        //
        //    Output, int FACE_ORDER[FACE_NUM], the number of vertices per face.
        //
        //    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.
        //
    {
        int DIM_NUM = 3;

        double a = Math.Sqrt(1.0 / 3.0);

        int[] face_order_save =
        {
            4, 4, 4, 4, 4, 4
        };
        int[] face_point_save =
        {
            1, 4, 3, 2,
            1, 2, 6, 5,
            2, 3, 7, 6,
            3, 4, 8, 7,
            1, 5, 8, 4,
            5, 6, 7, 8
        };
        double[] point_coord_save =
        {
            -a, -a, -a,
            a, -a, -a,
            a, a, -a,
            -a, a, -a,
            -a, -a, a,
            a, -a, a,
            a, a, a,
            -a, a, a
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);
    }

    public static void cube_size_3d(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE_SIZE_3D gives "sizes" for a cube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int *POINT_NUM, the number of points.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
        //    Output, int *FACE_NUM, the number of faces.
        //
        //    Output, int *FACE_ORDER_MAX, the maximum order of any face.
        //
    {
        point_num = 8;
        edge_num = 12;
        face_num = 6;
        face_order_max = 4;
    }

    public static double cube01_volume()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBE01_VOLUME returns the volume of the unit cube in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CUBE01_VOLUME, the volume.
        //
    {
        const double volume = 1.0;

        return volume;
    }

}