namespace TriangulationTest
{
    public static partial class TriangulationSampleData
    {
        public static void triangulation_order3_example1(int node_num, int triangle_num,
            ref double[] node_xy, ref int[] triangle_node, ref int[] triangle_neighbor )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_EXAMPLE1 sets up a sample triangulation.
        //
        //  Discussion:
        //
        //    This triangulation is actually a Delaunay triangulation.
        //
        //    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
        //    determined by calling TRIANGULATION_ORDER3_EXAMPLE1_SIZE first.
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
        //    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
        //    the triangles.
        //
        //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
        //    on each side.  Negative values indicate edges that lie on the exterior.
        //
        {
            int DIM_NUM = 2;
            int NODE_NUM = 13;
            int TRIANGLE_NUM = 16;
            int TRIANGLE_ORDER = 3;

            int i;
            int[] triangle_neighbor_save = {
                -4, -13, 2,
                1, 4, 3,
                2, 5, 7,
                2, -43, 8,
                3, 8, 6,
                5, 9, 7,
                3, 6, -3,
                5, 4, 10,
                6, 10, 12,
                9, 8, 11,
                12, 10, 14,
                9, 11, 13,
                -23, 12, 16,
                11, -47, 15,
                16, 14, -50,
                13, 15, -39
            }
            ;
            int[] triangle_node_save = {
                3, 4, 1,
                3, 1, 2,
                3, 2, 8,
                2, 1, 5,
                8, 2, 13,
                8, 13, 9,
                3, 8, 9,
                13, 2, 5,
                9, 13, 7,
                7, 13, 5,
                6, 7, 5,
                9, 7, 6,
                10, 9, 6,
                6, 5, 12,
                11, 6, 12,
                10, 6, 11
            }
            ;
            double[] node_xy_save = {
                0.0, 0.0,
                2.0, 2.0,
                -1.0, 3.0,
                -2.0, 2.0,
                8.0, 2.0,
                9.0, 5.0,
                7.0, 4.0,
                5.0, 6.0,
                6.0, 7.0,
                8.0, 8.0,
                11.0, 7.0,
                10.0, 4.0,
                6.0, 4.0
            }
            ;

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

        public static void triangulation_order3_example1_size(ref int node_num, ref int triangle_num,
                ref int hole_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_EXAMPLE1_SIZE sets sizes for a sample triangulation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2005
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
            node_num = 13;
            triangle_num = 16;
            hole_num = 0;
        }
    }
}