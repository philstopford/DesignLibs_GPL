using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class VertexCount
    {
        public static int triangulation_order6_vertex_count(int tri_num, int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER6_VERTEX_COUNT counts vertex nodes in a triangulation.
            //
            //  Discussion:
            //
            //    In a triangulation of order 6, some nodes are midside nodes and some
            //    nodes are vertex nodes.
            //
            //    Especially when an order 6 triangulation is used to handle the
            //    Navier Stokes equations, it is useful to know the number of
            //    vertex and midside nodes.
            //
            //    Note that the number of midside nodes is simple NODE_NUM - VERTEX_NUM.
            //
            //  Diagram:
            //
            //       3
            //    s  |.
            //    i  | .
            //    d  |  .
            //    e  6   5  side 2
            //       |    .
            //    3  |     .
            //       |      .
            //       1---4---2
            //
            //         side 1
            //
            //    The local node numbering.  Local nodes 1, 2 and 3 are "vertex nodes",
            //    while nodes 4, 5 and 6 are "midside nodes".
            //
            //
            //   21-22-23-24-25
            //    |.    |.    |
            //    | .   | .   |
            //   16 17 18 19 20
            //    |   . |   . |
            //    |    .|    .|
            //   11-12-13-14-15
            //    |.    |.    |
            //    | .   | .   |
            //    6  7  8  9 10
            //    |   . |   . |
            //    |    .|    .|
            //    1--2--3--4--5
            //
            //    A sample grid, which contains 9 vertex nodes and 16 midside nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters
            //
            //    Input, int TRI_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[6*TRI_NUM], lists the nodes that
            //    make up each triangle.  The first three nodes are the vertices,
            //    in counterclockwise order.  The fourth value is the midside
            //    node between nodes 1 and 2; the fifth and sixth values are
            //    the other midside nodes in the logical order.
            //
            //    Output, int TRIANGULATION_ORDER6_VERTEX_COUNT, the number of nodes
            //    which are vertices.
            //
        {
            int j;
            int vertex_num;
            int[] vertices;

            vertices = new int[3 * tri_num];

            for (j = 0; j < tri_num; j++)
            {
                vertices[j] = triangle_node[0 + j * 6];
            }

            for (j = 0; j < tri_num; j++)
            {
                vertices[tri_num + j] = triangle_node[1 + j * 6];
            }

            for (j = 0; j < tri_num; j++)
            {
                vertices[2 * tri_num + j] = triangle_node[2 + j * 6];
            }

            typeMethods.i4vec_sort_heap_a(3 * tri_num, ref vertices);

            vertex_num = typeMethods.i4vec_sorted_unique(3 * tri_num, ref vertices);
            
            return vertex_num;
        }
    }
}