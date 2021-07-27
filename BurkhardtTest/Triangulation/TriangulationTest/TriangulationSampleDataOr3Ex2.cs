namespace TriangulationTest
{
    public static partial class TriangulationSampleData
    {
        public static void triangulation_order3_example2(int node_num, int triangle_num,
                ref double[] node_xy, ref int[] triangle_node, ref int[] triangle_neighbor)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_EXAMPLE2 sets up a sample triangulation.
            //
            //  Discussion:
            //
            //    This triangulation is actually a Delaunay triangulation.
            //
            //    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
            //    determined by calling TRIANGULATION_ORDER3_EXAMPLE2_SIZE first.
            //
            //  Diagram:
            //
            //   21-22-23-24-25
            //    |. |. |. |. |
            //    | .| .| .| .|
            //   16-17-18-19-20
            //    |. |. |. |. |
            //    | .| .| .| .|
            //   11-12-13-14-15
            //    |. |. |. |. |
            //    | .| .| .| .|
            //    6--7--8--9-10
            //    |. |. |. |. |
            //    | .| .| .| .|
            //    1--2--3--4--5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
            //
            //    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
            //    triangles.
            //
            //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
            //    on each side.  Negative values indicate edges that lie on the exterior.
            //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 25;
            int TRIANGLE_NUM = 32;
            int TRIANGLE_ORDER = 3;

            int i;
            int[] triangle_neighbor_save =
            {
                -1, 2, -1,
                9, 1, 3,
                -1, 4, 2,
                11, 3, 5,
                -1, 6, 4,
                13, 5, 7,
                -1, 8, 6,
                15, 7, -1,
                2, 10, -1,
                17, 9, 11,
                4, 12, 10,
                19, 11, 13,
                6, 14, 12,
                21, 13, 15,
                8, 16, 14,
                23, 15, -1,
                10, 18, -1,
                25, 17, 19,
                12, 20, 18,
                27, 19, 21,
                14, 22, 20,
                29, 21, 23,
                16, 24, 22,
                31, 23, -1,
                18, 26, -1,
                -1, 25, 27,
                20, 28, 26,
                -1, 27, 29,
                22, 30, 28,
                -1, 29, 31,
                24, 32, 30,
                -1, 31, -1
            };
            int[] triangle_node_save =
            {
                1, 2, 6,
                7, 6, 2,
                2, 3, 7,
                8, 7, 3,
                3, 4, 8,
                9, 8, 4,
                4, 5, 9,
                10, 9, 5,
                6, 7, 11,
                12, 11, 7,
                7, 8, 12,
                13, 12, 8,
                8, 9, 13,
                14, 13, 9,
                9, 10, 14,
                15, 14, 10,
                11, 12, 16,
                17, 16, 12,
                12, 13, 17,
                18, 17, 13,
                13, 14, 18,
                19, 18, 14,
                14, 15, 19,
                20, 19, 15,
                16, 17, 21,
                22, 21, 17,
                17, 18, 22,
                23, 22, 18,
                18, 19, 23,
                24, 23, 19,
                19, 20, 24,
                25, 24, 20
            };
            double[] node_xy_save =
            {
                0.0, 0.0,
                1.0, 0.0,
                2.0, 0.0,
                3.0, 0.0,
                4.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                2.0, 1.0,
                3.0, 1.0,
                4.0, 1.0,
                0.0, 2.0,
                1.0, 2.0,
                2.0, 2.0,
                3.0, 2.0,
                4.0, 2.0,
                0.0, 3.0,
                1.0, 3.0,
                2.0, 3.0,
                3.0, 3.0,
                4.0, 3.0,
                0.0, 4.0,
                1.0, 4.0,
                2.0, 4.0,
                3.0, 4.0,
                4.0, 4.0
            };

            for (i = 0; i < 3 * TRIANGLE_NUM; i++)
            {
                triangle_neighbor[i] = triangle_neighbor_save[i];
            }

            for (i = 0; i < TRIANGLE_ORDER * TRIANGLE_NUM; i++)
            {
                triangle_node[i] = triangle_node_save[i];
            }

            for (i = 0; i < DIM_NUM * NODE_NUM; i++)
            {
                node_xy[i] = node_xy_save[i];
            }
        }

        public static void triangulation_order3_example2_size(ref int node_num, ref int triangle_num,
                ref int hole_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_EXAMPLE2_SIZE sets sizes for a sample triangulation.
            //
            //  Diagram:
            //
            //   21-22-23-24-25
            //    |. |. |. |. |
            //    | .| .| .| .|
            //   16-17-18-19-20
            //    |. |. |. |. |
            //    | .| .| .| .|
            //   11-12-13-14-15
            //    |. |. |. |. |
            //    | .| .| .| .|
            //    6--7--8--9-10
            //    |. |. |. |. |
            //    | .| .| .| .|
            //    1--2--3--4--5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *NODE_NUM, the number of nodes.
            //
            //    Output, int *TRIANGLE_NUM, the number of triangles.
            //
            //    Output, int *HOLE_NUM, the number of holes.
            //
        {
            node_num = 25;
            triangle_num = 32;
            hole_num = 0;
        }
    }
}