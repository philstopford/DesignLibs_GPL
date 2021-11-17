using System;
using Burkardt.Types;

namespace Burkardt.Dodecahedron;

public static class Geometry
{
    public static void dodec_shape_3d(int point_num, int face_num, int face_order_max,
            ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DODEC_SHAPE_3D describes a dodecahedron in 3D.
        //
        //  Discussion:
        //
        //    The vertices lie on the unit sphere.
        //
        //    The dual of a dodecahedron is the icosahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2003
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

        double phi = 0.5 * (Math.Sqrt(5.0) + 1.0);

        double a = 1.0 / Math.Sqrt(3.0);
        double b = phi / Math.Sqrt(3.0);
        double c = (phi - 1.0) / Math.Sqrt(3.0);
        double z = 0.0;

        int[] face_order_save =
        {
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5
        };
        int[] face_point_save =
        {
            2, 9, 1, 13, 14,
            5, 10, 6, 16, 15,
            3, 11, 4, 14, 13,
            8, 12, 7, 15, 16,
            3, 13, 1, 17, 18,
            2, 14, 4, 20, 19,
            5, 15, 7, 18, 17,
            8, 16, 6, 19, 20,
            5, 17, 1, 9, 10,
            3, 18, 7, 12, 11,
            2, 19, 6, 10, 9,
            8, 20, 4, 11, 12
        };
        double[] point_coord_save =
        {
            a, a, a,
            a, a, -a,
            a, -a, a,
            a, -a, -a,
            -a, a, a,
            -a, a, -a,
            -a, -a, a,
            -a, -a, -a,
            c, b, z,
            -c, b, z,
            c, -b, z,
            -c, -b, z,
            b, z, c,
            b, z, -c,
            -b, z, c,
            -b, z, -c,
            z, c, b,
            z, -c, b,
            z, c, -b,
            z, -c, -b
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);

    }

    public static void dodec_size_3d(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DODEC_SIZE_3D gives "sizes" for a dodecahedron in 3D.
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
        point_num = 20;
        edge_num = 30;
        face_num = 12;
        face_order_max = 5;

    }
}