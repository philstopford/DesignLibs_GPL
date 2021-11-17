namespace TriangulationTest;

public static partial class TriangulationSampleData
{
    public static void triangulation_order6_example1(int node_num, int triangle_num,
            ref double[] node_xy, ref int[] triangle_node, ref int[] triangle_neighbor )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_EXAMPLE1 sets up a sample triangulation.
        //
        //  Discussion:
        //
        //    This triangulation is actually a Delaunay triangulation.
        //
        //    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
        //    determined by calling TRIANGULATION_ORDER6_EXAMPLE1_SIZE first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2005
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
        //    Output, int TRIANGLE_ORDER[6*TRIANGLE_NUM], the nodes that make up
        //    the triangles.
        //
        //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
        //    on each side.  Negative values indicate edges that lie on the exterior.
        //
    {
        int DIM_NUM = 2;
        int NODE_NUM = 48;
        int TRIANGLE_NUM = 16;
        int TRIANGLE_ORDER = 6;

        int i;
        int j;
        double[] node_xy_save = {
                0.0, 0.0,
                1.0, 0.0,
                2.0, 0.0,
                3.0, 0.0,
                4.0, 0.0,
                5.0, 0.0,
                6.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                2.0, 1.0,
                3.0, 1.0,
                4.0, 1.0,
                5.0, 1.0,
                6.0, 1.0,
                0.0, 2.0,
                1.0, 2.0,
                2.0, 2.0,
                3.0, 2.0,
                4.0, 2.0,
                5.0, 2.0,
                6.0, 2.0,
                0.0, 3.0,
                1.0, 3.0,
                2.0, 3.0,
                3.0, 3.0,
                5.0, 3.0,
                6.0, 3.0,
                0.0, 4.0,
                1.0, 4.0,
                2.0, 4.0,
                3.0, 4.0,
                4.0, 4.0,
                5.0, 4.0,
                6.0, 4.0,
                0.0, 5.0,
                1.0, 5.0,
                2.0, 5.0,
                3.0, 5.0,
                4.0, 5.0,
                5.0, 5.0,
                6.0, 5.0,
                0.0, 6.0,
                1.0, 6.0,
                2.0, 6.0,
                3.0, 6.0,
                4.0, 6.0,
                5.0, 6.0,
                6.0, 6.0
            }
            ;
        int[] triangle_node_save = {
                1, 3, 15, 2, 9, 8,
                17, 15, 3, 16, 9, 10,
                5, 19, 3, 12, 11, 4,
                17, 3, 19, 10, 11, 18,
                7, 21, 5, 14, 13, 6,
                19, 5, 21, 12, 13, 20,
                17, 30, 15, 24, 23, 16,
                28, 15, 30, 22, 23, 29,
                30, 17, 32, 24, 25, 31,
                21, 34, 19, 27, 26, 20,
                30, 44, 28, 37, 36, 29,
                42, 28, 44, 35, 36, 43,
                32, 46, 30, 39, 38, 31,
                44, 30, 46, 37, 38, 45,
                32, 34, 46, 33, 40, 39,
                48, 46, 34, 47, 40, 41
            }
            ;
        int[] triangle_neighbor_save = {
                -3, 2, -5,
                7, 1, 4,
                6, 4, -11,
                2, 3, -14,
                -15, 6, -17,
                3, 5, 10,
                9, 8, 2,
                -24, 7, 11,
                7, -28, 13,
                27, -31, 6,
                8, 14, 12,
                -36, 11, -38,
                15, 14, 9,
                11, 13, -44,
                -45, 16, 13,
                -48, 15, -50
            }
            ;

        for (j = 0; j < NODE_NUM; j++)
        {
            for (i = 0; i < DIM_NUM; i++)
            {
                node_xy[i + j * DIM_NUM] = node_xy_save[i + j * DIM_NUM];
            }
        }

        for (j = 0; j < TRIANGLE_NUM; j++)
        {
            for (i = 0; i < TRIANGLE_ORDER; i++)
            {
                triangle_node[i + j * TRIANGLE_ORDER] =
                    triangle_node_save[i + j * TRIANGLE_ORDER];
            }
        }

        for (j = 0; j < TRIANGLE_NUM; j++)
        {
            for (i = 0; i < 3; i++)
            {
                triangle_neighbor[i + j * 3] = triangle_neighbor_save[i + j * 3];
            }
        }
    }

    public static void triangulation_order6_example1_size(ref int node_num, ref int triangle_num,
            ref int hole_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_EXAMPLE1_SIZE sets sizes for a sample triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2004
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
        node_num = 48;
        triangle_num = 16;
        hole_num = 1;
    }
}