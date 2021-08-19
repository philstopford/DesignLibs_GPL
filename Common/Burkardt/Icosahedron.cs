using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Icosahedron
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
            int DIM_NUM = 3;
            int EDGE_NUM = 30;
            int EDGE_ORDER = 2;
            int FACE_NUM = 20;
            int POINT_NUM = 12;

            double phi = 0.5 * (Math.Sqrt(5.0) + 1.0);

            double a = phi / Math.Sqrt(1.0 + phi * phi);
            double b = 1.0 / Math.Sqrt(1.0 + phi * phi);
            double z = 0.0;

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
        
        public static void icos_size ( ref int point_num, ref int edge_num, ref int face_num, 
                ref int face_order_max )

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
    }
}