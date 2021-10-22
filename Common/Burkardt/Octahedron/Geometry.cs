using Burkardt.Types;

namespace Burkardt.Octahedron
{
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
            int DIM_NUM = 3;

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

    }
}