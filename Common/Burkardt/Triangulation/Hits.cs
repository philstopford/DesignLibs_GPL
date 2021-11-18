namespace Burkardt.TriangulationNS;

public static class Hits
{
    public static int triangulation_hits(int node_num, double[] node_xy, int triangle_order,
            int triangle_num, int[] triangle_node, int data_num, double[] data_xy,
            ref int[] triangle_hit)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_HITS counts "data hits" in triangles and a triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes in the 
        //    triangulation.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_ORDER, the order of triangles in 
        //    the triangulation.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in 
        //    the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], 
        //    the nodes making up each triangle.
        //
        //    Input, int DATA_NUM, the number of data points.
        //
        //    Input, double DATA_XY[2*DATA_NUM], the data points.
        //
        //    Output, int TRIANGLE_HIT[TRIANGLE_NUM], the number of
        //    hits in each triangle.
        //
        //    Output, int TRIANGULATION_HITS, the number of hits in
        //    the triangulation.
        //
    {
        int triangle;

        int triangulation_hit = 0;

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            triangle_hit[triangle] = 0;

            int d;
            for (d = 0; d < data_num; d++)
            {
                bool inside = true;

                int j;
                for (j = 0; j < 3; j++)
                {
                    int k = (j + 1) % 3;

                    int nj = triangle_node[j + triangle * triangle_order];
                    int nk = triangle_node[k + triangle * triangle_order];

                    if (0.0 < (data_xy[0 + d * 2] - node_xy[0 + nj * 2])
                        * (node_xy[1 + nk * 2] - node_xy[1 + nj * 2])
                        - (data_xy[1 + d * 2] - node_xy[1 + nj * 2])
                        * (node_xy[0 + nk * 2] - node_xy[0 + nj * 2]))
                    {
                        inside = false;
                        break;
                    }
                }

                switch (inside)
                {
                    case true:
                        triangle_hit[triangle] += 1;
                        break;
                }
            }

            triangulation_hit += triangle_hit[triangle];
        }

        return triangulation_hit;
    }
}