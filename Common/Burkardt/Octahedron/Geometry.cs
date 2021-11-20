using Burkardt.Types;

namespace Burkardt.Octahedron;

public static class Geometry
{
    public static void octahedron_shape_3d(int point_num, int face_num, int face_order_max,
            ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OCTAHEDRON_SHAPE_3D describes an octahedron in 3D.
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
        const int DIM_NUM = 3;

        int[] face_order_save =
        {
            3, 3, 3, 3, 3, 3, 3, 3
        };
        int[] face_point_save =
        {
            1, 3, 2,
            1, 4, 3,
            1, 5, 4,
            1, 2, 5,
            2, 3, 6,
            3, 4, 6,
            4, 5, 6,
            5, 2, 6
        };
        double[] point_coord_save =
        {
            0.0, 0.0, -1.0,
            0.0, -1.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            -1.0, 0.0, 0.0,
            0.0, 0.0, 1.0
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);
    }

    public static void octahedron_size_3d(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OCTAHEDRON_SIZE_3D returns size information for an octahedron in 3D.
        //
        //  Discussion:
        //
        //    This routine can be called before calling OCTAHEDRON_SHAPE_3D,
        //    so that space can be allocated for the arrays.
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
        //    Output, int *FACE_ORDER_MAX, the maximum number of vertices
        //    per face.
        //
    {
        point_num = 6;
        edge_num = 12;
        face_num = 8;
        face_order_max = 3;

    }

    public static void truncated_octahedron_shape_3d(int point_num, int face_num,
            int face_order_max, ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_OCTAHEDRON_SHAPE_3D describes a truncated octahedron in 3D.
        //
        //  Discussion:
        //
        //    The shape is a truncated octahedron.  There are 8 hexagons and 6
        //    squares.
        //
        //    The truncated octahedron is an interesting shape because it
        //    is "space filling".  In other words, all of 3D space can be
        //    filled by a regular lattice of these shapes.
        //
        //    Call TRUNCATED_OCTAHEDRON_SIZE_3D to get the values of POINT_NUM,
        //    FACE_NUM, and FACE_ORDER_MAX, so you can allocate space for the arrays.
        //
        //    For each face, the face list must be of length FACE_ORDER_MAX.
        //    In cases where a face is of lower than maximum order (the
        //    squares, in this case), the extra entries are listed as "-1".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of points (24).
        //
        //    Input, int FACE_NUM, the number of faces (14).
        //
        //    Input, int FACE_ORDER_MAX, the maximum order of any face (6).
        //
        //    Output, double POINT_COORD[3*POINT_NUM], the vertices.
        //
        //    Output, int FACE_ORDER[FACE_NUM], the number of vertices per face.
        //
        //    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.
        //
    {
        const int DIM_NUM = 3;

        int[] face_order_save =
        {
            4, 4, 4, 4, 4, 4, 6, 6, 6, 6,
            6, 6, 6, 6
        };
        int[] face_point_save =
        {
            17, 11, 9, 15, -1, -1,
            14, 8, 10, 16, -1, -1,
            22, 24, 21, 18, -1, -1,
            12, 5, 2, 6, -1, -1,
            13, 19, 23, 20, -1, -1,
            4, 1, 3, 7, -1, -1,
            19, 13, 7, 3, 8, 14,
            15, 9, 4, 7, 13, 20,
            16, 10, 5, 12, 18, 21,
            22, 18, 12, 6, 11, 17,
            20, 23, 24, 22, 17, 15,
            14, 16, 21, 24, 23, 19,
            9, 11, 6, 2, 1, 4,
            3, 1, 2, 5, 10, 8
        };
        double[] point_coord_save =
        {
            -1.5, -0.5, 0.0,
            -1.5, 0.5, 0.0,
            -1.0, -1.0, -0.70710677,
            -1.0, -1.0, 0.70710677,
            -1.0, 1.0, -0.70710677,
            -1.0, 1.0, 0.70710677,
            -0.5, -1.5, 0.0,
            -0.5, -0.5, -1.4142135,
            -0.5, -0.5, 1.4142135,
            -0.5, 0.5, -1.4142135,
            -0.5, 0.5, 1.4142135,
            -0.5, 1.5, 0.0,
            0.5, -1.5, 0.0,
            0.5, -0.5, -1.4142135,
            0.5, -0.5, 1.4142135,
            0.5, 0.5, -1.4142135,
            0.5, 0.5, 1.4142135,
            0.5, 1.5, 0.0,
            1.0, -1.0, -0.70710677,
            1.0, -1.0, 0.70710677,
            1.0, 1.0, -0.70710677,
            1.0, 1.0, 0.70710677,
            1.5, -0.5, 0.0,
            1.5, 0.5, 0.0
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);

    }

    public static void truncated_octahedron_size_3d(ref int point_num, ref int edge_num,
            ref int face_num, ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRUNCATED_OCTAHEDRON_SIZE_3D gives "sizes" for a truncated octahedron in 3D.
        //
        //  Discussion:
        //
        //    The truncated octahedron is "space-filling".
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
        point_num = 24;
        edge_num = 36;
        face_num = 14;
        face_order_max = 6;
    }


}