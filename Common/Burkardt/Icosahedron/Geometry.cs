using System;
using Burkardt.Types;

namespace Burkardt.Icosahedron;

public static class Geometry
{
    public static void icos_shape(int point_num, int edge_num, int face_num,
            int face_order_max, ref double[] point_coord, ref int[] edge_point, ref int[] face_order,
            ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ICOS_SHAPE describes a icosahedron.
        //
        //  Discussion:
        //
        //    The input data required for this routine can be retrieved from ICOS_NUM.
        //
        //    The vertices lie on the unit sphere.
        //
        //    The dual of an icosahedron is the dodecahedron.
        //
        //    The data has been rearranged from a previous assignment.  
        //    The STRIPACK program refuses to triangulate data if the first
        //    three nodes are "collinear" on the sphere.
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
        //    Input, int POINT_NUM, the number of points (12).
        //
        //    Input, int EDGE_NUM, the number of edges (30).
        //
        //    Input, int FACE_NUM, the number of faces (20).
        //
        //    Input, int FACE_ORDER_MAX, the maximum number of vertices 
        //    per face (3).
        //
        //    Output, double POINT_COORD[3*POINT_NUM], the point coordinates.
        //
        //    Output, int EDGE_POINT[2*EDGE_NUM], the points that make up each 
        //    edge, listed in ascending order of their indexes.
        //
        //    Output, int FACE_ORDER[FACE_NUM], the number of vertices per face.
        //
        //    Output, int FACE_POINT[FACE_ORDER_MAX*FACE_NUM]; FACE_POINT(I,J)
        //    contains the index of the I-th point in the J-th face.  The
        //    points are listed in the counter clockwise direction defined
        //    by the outward normal at the face.  The nodes of each face are 
        //    ordered so that the lowest index occurs first.  The faces are 
        //    then sorted by nodes.
        //
    {
        const int DIM_NUM = 3;
        const int EDGE_ORDER = 2;

        double phi = 0.5 * (Math.Sqrt(5.0) + 1.0);

        double a = phi / Math.Sqrt(1.0 + phi * phi);
        double b = 1.0 / Math.Sqrt(1.0 + phi * phi);
        const double z = 0.0;

        int[] edge_point_save =
        {
            0, 1,
            0, 2,
            0, 3,
            0, 4,
            0, 5,
            1, 2,
            1, 3,
            1, 6,
            1, 7,
            2, 4,
            2, 6,
            2, 8,
            3, 5,
            3, 7,
            3, 9,
            4, 5,
            4, 8,
            4, 10,
            5, 9,
            5, 10,
            6, 7,
            6, 8,
            6, 11,
            7, 9,
            7, 11,
            8, 10,
            8, 11,
            9, 10,
            9, 11,
            10, 11
        };
        int[] face_order_save =
        {
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3
        };
        int[] face_point_save =
        {
            0, 1, 3,
            0, 2, 1,
            0, 3, 5,
            0, 4, 2,
            0, 5, 4,
            1, 2, 6,
            1, 6, 7,
            1, 7, 3,
            2, 4, 8,
            2, 8, 6,
            3, 7, 9,
            3, 9, 5,
            4, 5, 10,
            4, 10, 8,
            5, 9, 10,
            6, 8, 11,
            6, 11, 7,
            7, 11, 9,
            8, 10, 11,
            9, 11, 10
        };
        double[] point_coord_save =
        {
            a, b, z,
            a, -b, z,
            b, z, a,
            b, z, -a,
            z, a, b,
            z, a, -b,
            z, -a, b,
            z, -a, -b,
            -b, z, a,
            -b, z, -a,
            -a, b, z,
            -a, -b, z
        };

        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);
        typeMethods.i4vec_copy(EDGE_ORDER * edge_num, edge_point_save, ref edge_point);
        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);

    }

    public static void icos_num(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ICOS_NUM gives "sizes" for an icosahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2007
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
        point_num = 12;
        edge_num = 30;
        face_num = 20;
        face_order_max = 3;
    }

    public static void icos_size(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ICOS_SIZE gives "sizes" for an icosahedron in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2007
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
        point_num = 12;
        edge_num = 30;
        face_num = 20;
        face_order_max = 3;
    }

    public static void soccer_shape_3d(int point_num, int face_num, int face_order_max,
            ref double[] point_coord, ref int[] face_order, ref int[] face_point)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOCCER_SHAPE_3D describes a truncated icosahedron in 3D.
        //
        //  Discussion:
        //
        //    The shape is a truncated icosahedron, which is the design used
        //    on a soccer ball.  There are 12 pentagons and 20 hexagons.
        //
        //    Call SOCCER_SIZE_3D to get the values of POINT_NUM, FACE_NUM, and
        //    FACE_ORDER_MAX, so you can allocate space for the arrays.
        //
        //    For each face, the face list must be of length FACE_ORDER_MAX.
        //    In cases where a face is of lower than maximum order (the
        //    12 pentagons, in this case), the extra entries are listed as
        //    "-1".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    http://mathworld.wolfram.com/TruncatedIcosahedron.html
        //
        //  Parameters:
        //
        //    Input, int POINT_NUM, the number of points (60).
        //
        //    Input, int FACE_NUM, the number of faces (32).
        //
        //    Input, int FACE_ORDER_MAX, the maximum order of any face (6).
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
            6, 6, 5, 6, 5, 6, 5, 6, 6, 6,
            5, 6, 5, 6, 5, 6, 6, 6, 5, 6,
            5, 5, 6, 6, 6, 5, 6, 5, 6, 6,
            5, 6
        };
        int[] face_point_save =
        {
            30, 43, 47, 41, 29, 23,
            30, 23, 12, 9, 15, 27,
            30, 27, 35, 45, 43, -1,
            43, 45, 53, 59, 56, 47,
            23, 29, 21, 11, 12, -1,
            27, 15, 13, 22, 33, 35,
            47, 56, 54, 44, 41, -1,
            45, 35, 33, 42, 51, 53,
            12, 11, 4, 1, 3, 9,
            29, 41, 44, 37, 25, 21,
            15, 9, 3, 6, 13, -1,
            56, 59, 60, 58, 55, 54,
            53, 51, 57, 60, 59, -1,
            11, 21, 25, 19, 10, 4,
            33, 22, 24, 36, 42, -1,
            13, 6, 7, 17, 24, 22,
            54, 55, 48, 39, 37, 44,
            51, 42, 36, 40, 50, 57,
            4, 10, 8, 2, 1, -1,
            3, 1, 2, 5, 7, 6,
            25, 37, 39, 28, 19, -1,
            55, 58, 52, 46, 48, -1,
            60, 57, 50, 49, 52, 58,
            10, 19, 28, 26, 16, 8,
            36, 24, 17, 20, 32, 40,
            7, 5, 14, 20, 17, -1,
            48, 46, 34, 26, 28, 39,
            50, 40, 32, 38, 49, -1,
            8, 16, 18, 14, 5, 2,
            46, 52, 49, 38, 31, 34,
            16, 26, 34, 31, 18, -1,
            32, 20, 14, 18, 31, 38
        };
        double[] point_coord_save =
        {
            -1.00714, 0.153552, 0.067258,
            -0.960284, 0.0848813, -0.33629,
            -0.95172, -0.153552, 0.33629,
            -0.860021, 0.529326, 0.150394,
            -0.858, -0.290893, -0.470806,
            -0.849436, -0.529326, 0.201774,
            -0.802576, -0.597996, -0.201774,
            -0.7842, 0.418215, -0.502561,
            -0.749174, -0.0848813, 0.688458,
            -0.722234, 0.692896, -0.201774,
            -0.657475, 0.597996, 0.502561,
            -0.602051, 0.290893, 0.771593,
            -0.583675, -0.692896, 0.470806,
            -0.579632, -0.333333, -0.771593,
            -0.52171, -0.418215, 0.771593,
            -0.505832, 0.375774, -0.803348,
            -0.489955, -0.830237, -0.33629,
            -0.403548, 0, -0.937864,
            -0.381901, 0.925138, -0.201774,
            -0.352168, -0.666667, -0.688458,
            -0.317142, 0.830237, 0.502561,
            -0.271054, -0.925138, 0.33629,
            -0.227464, 0.333333, 0.937864,
            -0.224193, -0.993808, -0.067258,
            -0.179355, 0.993808, 0.150394,
            -0.165499, 0.608015, -0.803348,
            -0.147123, -0.375774, 0.937864,
            -0.103533, 0.882697, -0.502561,
            -0.0513806, 0.666667, 0.771593,
            0.0000000, 0, 1.021,
            0.0000000, 0, -1.021,
            0.0513806, -0.666667, -0.771593,
            0.103533, -0.882697, 0.502561,
            0.147123, 0.375774, -0.937864,
            0.165499, -0.608015, 0.803348,
            0.179355, -0.993808, -0.150394,
            0.224193, 0.993808, 0.067258,
            0.227464, -0.333333, -0.937864,
            0.271054, 0.925138, -0.33629,
            0.317142, -0.830237, -0.502561,
            0.352168, 0.666667, 0.688458,
            0.381901, -0.925138, 0.201774,
            0.403548, 0, 0.937864,
            0.489955, 0.830237, 0.33629,
            0.505832, -0.375774, 0.803348,
            0.521710, 0.418215, -0.771593,
            0.579632, 0.333333, 0.771593,
            0.583675, 0.692896, -0.470806,
            0.602051, -0.290893, -0.771593,
            0.657475, -0.597996, -0.502561,
            0.722234, -0.692896, 0.201774,
            0.749174, 0.0848813, -0.688458,
            0.784200, -0.418215, 0.502561,
            0.802576, 0.597996, 0.201774,
            0.849436, 0.529326, -0.201774,
            0.858000, 0.290893, 0.470806,
            0.860021, -0.529326, -0.150394,
            0.951720, 0.153552, -0.33629,
            0.960284, -0.0848813, 0.33629,
            1.007140, -0.153552, -0.067258
        };

        typeMethods.i4vec_copy(face_num, face_order_save, ref face_order);
        typeMethods.i4vec_copy(face_order_max * face_num, face_point_save, ref face_point);
        typeMethods.r8vec_copy(DIM_NUM * point_num, point_coord_save, ref point_coord);

    }

    public static void soccer_size_3d(ref int point_num, ref int edge_num, ref int face_num,
            ref int face_order_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOCCER_SIZE_3D gives "sizes" for a truncated icosahedron in 3D.
        //
        //  Discussion:
        //
        //    The shape is a truncated icosahedron, which is the design used
        //    on a soccer ball.  There are 12 pentagons and 20 hexagons.
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
        //  Reference:
        //
        //    http://polyhedra.wolfram.com/uniform/u25.html
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
        point_num = 60;
        edge_num = 90;
        face_num = 32;
        face_order_max = 6;
    }

}