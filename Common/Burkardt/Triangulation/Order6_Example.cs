namespace Burkardt.TriangulationNS;

public static class Order6_Example
{
    public static void triangulation_order6_example2(int node_num, int triangle_num,
            ref double[] node_xy, ref int[] triangle_node, ref int[] triangle_neighbor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_EXAMPLE2 sets up a sample triangulation.
        //
        //  Discussion:
        //
        //    This triangulation is actually a Delaunay triangulation.
        //
        //    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
        //    determined by calling TRIANGULATION_ORDER6_EXAMPLE2_SIZE first.
        //
        //  Diagram:
        //
        //   21-22-23-24-25
        //    |\  6 |\  8 |
        //    | \   | \   |
        //   16 17 18 19 20
        //    |   \ |   \ |
        //    | 5  \| 7  \|
        //   11-12-13-14-15
        //    |\  2 |\  4 |
        //    | \   | \   |
        //    6  7  8  9 10
        //    | 1 \ | 3 \ |
        //    |    \|    \|
        //    1--2--3--4--5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2007
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
        const int DIM_NUM = 2;
        const int NODE_NUM = 25;
        const int TRIANGLE_NUM = 8;
        const int TRIANGLE_ORDER = 6;

        int i;
        int j;
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
        int[] triangle_node_save =
        {
            1, 3, 11, 2, 7, 6,
            13, 11, 3, 12, 7, 8,
            3, 5, 13, 4, 9, 8,
            15, 13, 5, 14, 9, 10,
            11, 13, 21, 12, 17, 16,
            23, 21, 13, 22, 17, 18,
            13, 15, 23, 14, 19, 18,
            25, 23, 15, 24, 19, 20
        };
        int[] triangle_neighbor_save =
        {
            -1, 2, -1,
            5, 1, 3,
            -1, 4, 2,
            7, 3, -1,
            2, 6, -1,
            -1, 5, 7,
            4, 8, 6,
            -1, 7, -1
        };

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
                triangle_node[i + j * TRIANGLE_ORDER] = triangle_node_save[i + j * TRIANGLE_ORDER];
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

    public static void triangulation_order6_example2_size(ref int node_num, ref int triangle_num,
            ref int hole_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_EXAMPLE2_SIZE sets sizes for a sample triangulation.
        //
        //  Diagram:
        //
        //   21-22-23-24-25
        //    |\  6 |\  8 |
        //    | \   | \   |
        //   16 17 18 19 20
        //    |   \ |   \ |
        //    | 5  \| 7  \|
        //   11-12-13-14-15
        //    |\  2 |\  4 |
        //    | \   | \   |
        //    6  7  8  9 10
        //    | 1 \ | 3 \ |
        //    |    \|    \|
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
        triangle_num = 8;
        hole_num = 0;
    }
}