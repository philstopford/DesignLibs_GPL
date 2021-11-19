using System;
using Burkardt.Types;

namespace Burkardt.QuadMesh;

public static class NodeOrder
{
    public static int[] node_order_q4_mesh ( int element_num, int[] element_node, int node_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NODE_ORDER_Q4_MESH determines the order of nodes in a Q4 mesh.
        //
        //  Discussion:
        //
        //    The order of a node is the number of elements that use that node
        //    as a vertex.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], 
        //    the nodes that make up the elements.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Output, int NODE_ORDER_Q4_MESH[NODE_NUM], the order of each node.
        //
    {
        int element;

        int[] node_order = typeMethods.i4vec_zero_new ( node_num );

        for ( element = 0; element < element_num; element++ )
        {
            int i;
            for ( i = 0; i < 4; i++ )
            {
                int node = element_node[i+element*4];
                if ( node < 0 || node_num <= node )
                {
                    Console.WriteLine("");
                    Console.WriteLine("NODE_ORDER_Q4_MESH - Fatal error!");
                    Console.WriteLine("  Illegal entry in ELEMENT_NODE.");
                    return null;
                }

                node_order[node] += 1;
            }
        }
        return node_order;
    }
}