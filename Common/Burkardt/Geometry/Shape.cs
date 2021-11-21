using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.Geometry;

public static class Shape
{
    public static double shape_point_dist_2d(double[] pc, double[] p1, int side_num,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHAPE_POINT_DIST_2D: distance ( regular shape, point ) in 2D.
        //
        //  Discussion:
        //
        //    The "regular shape" is assumed to be an equilateral and equiangular
        //    polygon, such as the standard square, pentagon, hexagon, and so on.
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
        //    Input, double PC[2], the center of the shape.
        //
        //    Input, double P1[2], the first vertex of the shape.
        //
        //    Input, int SIDE_NUM, the number of sides.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double SHAPE_POINT_DIST_2D, the distance from the point
        //    to the shape.
        //
    {
        const int DIM_NUM = 2;

        double dist;
        double[] pa = new double[DIM_NUM];
        double[] pb = new double[DIM_NUM];
        //
        //  Determine the angle subtended by a single side.
        //
        double sector_angle = 360.0 / side_num;
        //
        //  How long is the half-diagonal?
        //
        double radius = Math.Sqrt(Math.Pow(p1[0] - pc[0], 2) + Math.Pow(p1[1] - pc[1], 2));
        switch (radius)
        {
            //
            //  If the radius is zero, then the shape is a point and the computation is easy.
            //
            case 0.0:
                dist = Math.Sqrt(Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2));
                return dist;
        }

        switch (Math.Sqrt(Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2)))
        {
            //
            //  If the test point is at the center, then the computation is easy.
            //  The angle subtended by any side is ( 2 * PI / SIDE_NUM ) and the
            //  nearest distance is the midpoint of any such side.
            //
            case 0.0:
                dist = radius * Math.Cos(Math.PI / side_num);
                return dist;
        }

        //
        //  Determine the angle between the ray to the first corner,
        //  and the ray to the test point.
        //
        double angle = Angle.angle_deg_2d(p1, pc, p);
        //
        //  Determine the sector of the point.
        //
        int sector_index = (int) (angle / sector_angle) + 1;
        //
        //  Generate the two corner points that terminate the SECTOR-th side.
        //
        double angle2 = (sector_index - 1) * sector_angle;
        angle2 = Helpers.degrees_to_radians(angle2);

        Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pa);

        angle2 = sector_index * sector_angle;
        angle2 = Helpers.degrees_to_radians(angle2);

        Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pb);
        //
        //  Determine the distance from the test point to the line segment that
        //  is the SECTOR-th side.
        //
        dist = Segments.segment_point_dist_2d(pa, pb, p);

        return dist;
    }

    public static void shape_point_near_2d(double[] pc, double[] p1, int side_num,
            double[] p, ref double[] pn, ref double dist)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHAPE_POINT_NEAR_2D: nearest point ( regular shape, point ) in 2D.
        //
        //  Discussion:
        //
        //    The "regular shape" is assumed to be an equilateral and equiangular
        //    polygon, such as the standard square, pentagon, hexagon, and so on.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double PC[2], the center of the shape.
        //
        //    Input, double P1[2], the first vertex of the shape.
        //
        //    Input, int SIDE_NUM, the number of sides.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double PN[2], the point on the shape that is nearest
        //    to the given point.
        //
        //    Output, double *DIST, the distance between the points.
        //
    {
        const int DIM_NUM = 2;

        double angle;
        double t = 0;
        double[] pa = new double[DIM_NUM];
        double[] pb = new double[DIM_NUM];
        double[] pd = new double[DIM_NUM];
        //
        //  Determine the angle subtended by a single side.
        //
        double sector_angle = 360.0 / side_num;
        //
        //  How long is the half-diagonal?
        //
        double radius = Math.Sqrt(Math.Pow(p1[0] - pc[0], 2) + Math.Pow(p1[1] - pc[1], 2));
        switch (radius)
        {
            //
            //  If the radius is zero, then the shape is a point and the computation is easy.
            //
            case 0.0:
                typeMethods.r8vec_copy(DIM_NUM, pc, ref pn);
                dist = Math.Sqrt(Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2));
                return;
        }

        switch (Math.Sqrt(Math.Pow(p[0] - pc[0], 2) + Math.Pow(p[1] - pc[1], 2)))
        {
            //
            //  If the test point is at the center, then the computation is easy.
            //  The angle subtended by any side is ( 2 * PI / SIDE_NUM ) and the
            //  nearest distance is the midpoint of any such side.
            //
            case 0.0:
                angle = Math.PI / side_num;
                pd[0] = (p1[0] - pc[0]) * Math.Cos(angle)
                        + (p1[1] - pc[1]) * Math.Sin(angle);
                pd[1] = -(p1[0] - pc[0]) * Math.Sin(angle)
                        + (p1[1] - pc[1]) * Math.Cos(angle);
                pn[0] = pc[0] + pd[0] * Math.Cos(angle);
                pn[1] = pc[1] + pd[1] * Math.Cos(angle);
                dist = radius * Math.Cos(angle);
                return;
        }

        //
        //  Determine the angle between the ray to the first corner,
        //  and the ray to the test point.
        //
        angle = Angle.angle_deg_2d(p1, pc, p);
        //
        //  Determine the sector of the point.
        //
        int sector_index = (int) (angle / sector_angle) + 1;
        //
        //  Generate the two corner points that terminate the SECTOR-th side.
        //
        double angle2 = (sector_index - 1) * sector_angle;
        angle2 = Helpers.degrees_to_radians(angle2);

        Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pa);

        angle2 = sector_index * sector_angle;
        angle2 = Helpers.degrees_to_radians(angle2);

        Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pb);
        //
        //  Determine the point on the SECTOR-th side of the shape which is
        //  nearest.
        //
        Segments.segment_point_near_2d(pa, pb, p, ref pn, ref dist, ref t);

    }

    public static void shape_print_3d(int point_num, int face_num, int face_order_max,
            double[] point_coord, int[] face_order, int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHAPE_PRINT_3D prints information about a polyhedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 January 2007
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
        //    Input, int FACE_ORDER_MAX, the number of vertices per face.
        //
        //    Input, double POINT_COORD[DIM_NUM*POINT_NUM], the point coordinates.
        //
        //    Input, int FACE_ORDER[FACE_NUM], the number of vertices per face.
        //
        //    Input, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.
        //
    {
        int DIM_NUM = 3;

        int i;
        int j;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("SHAPE_PRINT_3D");
        Console.WriteLine("  Information about a polytope.");
        Console.WriteLine("");
        Console.WriteLine("  The number of vertices is " + point_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Vertices:");
        Console.WriteLine("");
        Console.WriteLine("     Index          X               Y               Z");
        Console.WriteLine("");
        for (j = 0; j < point_num; j++)
        {
            cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += point_coord[i + j * DIM_NUM].ToString("0.########").PadLeft(16);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  The number of faces is " + face_num + "");
        Console.WriteLine("  The maximum order of any face is " + face_order_max + "");
        Console.WriteLine("");
        Console.WriteLine("     Index     Order         Indices of Nodes in Face");
        cout = "";
        for (j = 1; j <= face_order_max; j++)
        {
            cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(8);
        }

        Console.WriteLine(cout);

        cout = "                      ";
        Console.WriteLine(cout);

        for (j = 0; j < face_num; j++)
        {
            cout = "  " + (j + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                        + "  " + face_order[j].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                        + "  ";
            for (i = 0; i < face_order[j]; i++)
            {
                cout += face_point[i + j * face_order_max].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);
        }
    }

    public static void shape_ray_int_2d(double[] pc, double[] p1, int side_num, double[] pa,
            double[] pb, ref double[] pint)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHAPE_RAY_INT_2D: intersection ( regular shape, ray ) in 2D.
        //
        //  Discussion:
        //
        //    The "regular shape" is assumed to be an equilateral and equiangular
        //    polygon, such as the standard square, pentagon, hexagon, and so on.
        //
        //    The origin of the ray is assumed to be inside the shape.  This
        //    guarantees that the ray will intersect the shape in exactly one point.
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
        //    Input, double PC[2], the center of the shape.
        //
        //    Input, double P1[2], the first vertex of the shape.
        //
        //    Input, int SIDE_NUM, the number of sides.
        //
        //    Input, double PA[2], the origin of the ray.
        //
        //    Input, double PB[2], a second point on the ray.
        //
        //    Output, double PI[2], the point on the shape intersected by the ray.
        //
    {
        const int DIM_NUM = 2;

        int ival = 0;
        double[] pv1 = new double[DIM_NUM];
        double[] pv2 = new double[DIM_NUM];
        int sector_index;
        //
        //  Warning!
        //  No check is made to ensure that the ray origin is inside the shape.
        //  These calculations are not valid if that is not true!
        //
        //  Determine the angle subtended by a single side.
        //
        double sector_angle = 360.0 / side_num;
        //
        //  How long is the half-diagonal?
        //
        double radius = Math.Sqrt(Math.Pow(p1[0] - pc[0], 2) + Math.Pow(p1[1] - pc[1], 2));
        switch (radius)
        {
            //
            //  If the radius is zero, refuse to continue.
            //
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("SHAPE_RAY_INT_2D - Fatal error!");
                Console.WriteLine("  The shape has radius zero.");
                return;
        }

        //
        //  Determine which sector side intersects the ray.
        //
        pv2[0] = 0.0;
        pv2[1] = 0.0;

        for (sector_index = 1; sector_index <= side_num; sector_index++)
        {
            double angle2;
            switch (sector_index)
            {
                //
                //  Determine the two vertices that define this sector.
                //
                case 1:
                    angle2 = ((double) sector_index - 1) * sector_angle;
                    angle2 = Helpers.degrees_to_radians(angle2);

                    Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pv1);
                    break;
                default:
                    typeMethods.r8vec_copy(DIM_NUM, pv2, ref pv1);
                    break;
            }

            angle2 = sector_index * sector_angle;
            angle2 = Helpers.degrees_to_radians(angle2);

            Vector.Geometry.vector_rotate_base_2d(p1, pc, angle2, ref pv2);
            //
            //  Draw the angle from one vertex to the ray origin to the next vertex,
            //  and see if that angle contains the ray.  If so, then the ray
            //  must intersect the shape side of that sector.
            //
            bool inside = Angle.angle_contains_ray_2d(pv1, pa, pv2, pb);

            switch (inside)
            {
                case true:
                    //
                    //  Determine the intersection of the lines defined by the ray and the
                    //  sector side.  (We're already convinced that the ray and sector line
                    //  segment intersect, so we can use the simpler code that treats them
                    //  as full lines).
                    //
                    LineNS.Geometry.lines_exp_int_2d(pa, pb, pv1, pv2, ref ival, ref pint);

                    return;
            }
        }

        //
        //  If the calculation fell through the loop, then something's wrong.
        //
        Console.WriteLine("");
        Console.WriteLine("SHAPE_RAY_INT_2D - Fatal error!");
        Console.WriteLine("  Cannot find intersection of ray and shape.");
    }

}