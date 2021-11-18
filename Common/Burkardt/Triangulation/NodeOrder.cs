using System;

namespace Burkardt.TriangulationNS;

public static class NodeOrder
{
    public static int[] triangulation_node_order ( int triangle_order, int triangle_num,
            int[] triangle_node, int node_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_NODE_ORDER determines the order of nodes in a triangulation.
        //
        //  Discussion:
        //
        //    The order of a node is the number of triangles that use that node
        //    as a vertex.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer TRIANGLE_ORDER, the order of the triangulation.
        //
        //    Input, integer TRIANGLE_NUM, the number of triangles.
        //
        //    Input, integer TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes
        //    that make up the triangles.
        //
        //    Input, integer NODE_NUM, the number of nodes.
        //
        //    Output, integer TRIANGULATION_NODE_ORDER[NODE_NUM], the order of
        //    each node.
        //
    {
        int node;
        int triangle;

        int[] node_order = new int[node_num];

        for ( node = 0; node < node_num; node++ )
        {
            node_order[node] = 0;
        }

        for ( triangle = 0; triangle < triangle_num; triangle++ )
        {
            int i;
            for ( i = 0; i < triangle_order; i++ )
            {
                node = triangle_node[i+triangle*triangle_order];

                if ( node < 1 || node_num < node )
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_NODE_ORDER - Fatal error!");
                    Console.WriteLine("  Illegal entry in TRIANGLE_NODE.");
                    return null;
                }

                node_order[node-1] += 1;
            }
        }

        return node_order;
    }
}