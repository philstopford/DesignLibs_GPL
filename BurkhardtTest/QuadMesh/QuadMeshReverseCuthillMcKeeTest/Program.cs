using System;
using Burkardt.MatrixNS;
using Burkardt.MeshNS;
using Burkardt.QuadMesh;
using Burkardt.Table;
using Burkardt.Types;

namespace QuadMeshReverseCuthillMcKeeTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for QUAD_MESH_RCM.
        //
        //  Discussion:
        //
        //    QUAD_MESH_RCM applies the RCM reordering to a quadrilateral mesh.
        //
        //    The user supplies a node file and an element file, containing
        //    the coordinates of the nodes, and the indices of the nodes that
        //    make up each element.  
        //
        //    The program reads the data, computes the adjacency information,
        //    carries out the RCM algorithm to get the permutation, applies
        //    the permutation to the nodes and elements, and writes out
        //    new node and element files that correspond to the RCM permutation.
        //
        //    Note that node data is normally two dimensional, that is,
        //    each node has an X and Y coordinate.  In some applications, it
        //    may be desirable to specify more information.  This program
        //    will accept node data that includes DIM_NUM entries on each line,
        //    as long as DIM_NUM is the same for each entry.  
        //
        //  Usage:
        //
        //    quad_mesh_rcm prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txt contains the element definitions.
        //    * prefix_rcm_nodes.txt will contain the RCM node coordinates,
        //    * prefix_rcm_elements.txt will contain the RCM element definitions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        string prefix;

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("QUAD_MESH_RCM");
        Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
        Console.WriteLine("  Read an associated quad mesh dataset of ELEMENT_NUM");
        Console.WriteLine("  4 node quaderilaterals.");
        Console.WriteLine("");
        Console.WriteLine("  Apply the RCM reordering (Reverse Cuthill-McKee).");
        Console.WriteLine("");
        Console.WriteLine("  Reorder the data and write it out to files.");
        Console.WriteLine("");
        //
        //  Get the filename prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("QUAD_MESH_RCM");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string node_rcm_filename = prefix + "_rcm_nodes.txt";
        string element_rcm_filename = prefix + "_rcm_elements.txt";
        //
        //  Read the data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int dim_num = h.m;
        int node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of nodes NODE_NUM =  " + node_num + "");

        double[] node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1,
            dim_num, 5, "  Coordinates of first 5 nodes:");

        h = typeMethods.i4mat_header_read(element_filename);
        int element_order = h.m;
        int element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order  = " + element_order + "");
        Console.WriteLine("  Number of elements = " + element_num + "");

        if (element_order != 4)
        {
            Console.WriteLine("");
            Console.WriteLine("QUAD_MESH_RCM - Fatal error!");
            Console.WriteLine("  Data is not for 4-node quadrilaterals.");
            return;
        }

        int[] element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node, 1, 1,
            element_order, 5, "  First 5 elements:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Create the element neighbor array.
        //
        int[] element_neighbor = Neighbors.neighbor_elements_q4_mesh(element_num, element_node);
        //
        //  Count the number of adjacencies, and set up the ADJ_ROW 
        //  adjacency pointer array.
        //
        int[] adj_row = new int [node_num + 1];

        int adj_num = Adjacency.adj_size_q4_mesh(node_num, element_num,
            element_node, element_neighbor, ref adj_row);
        //
        //  Set up the ADJ adjacency array.
        //
        int[] adj = Adjacency.adj_set_q4_mesh(node_num, element_num, element_node,
            element_neighbor, adj_num, adj_row);

        int bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

        Console.WriteLine("");
        Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");
        //
        //  Compute the RCM permutation.
        //
        int[] perm = new int[node_num];
        //
        //  For now, add 1 to ADJ and ADJ_ROW since GENRCM assumes 1-based indexing.
        //
        for (i = 0; i < node_num + 1; i++)
        {
            adj_row[i] += 1;
        }

        for (i = 0; i < adj_num; i++)
        {
            adj[i] += 1;
        }

        perm = Burkardt.Graph.GenRCM.genrcm(node_num, adj_num, adj_row, adj);
        //
        //  On return, subtract 1 from ADJ, ADJROW and PERM.
        //
        for (i = 0; i < node_num; i++)
        {
            perm[i] -= 1;
        }

        for (i = 0; i < node_num + 1; i++)
        {
            adj_row[i] -= 1;
        }

        for (i = 0; i < adj_num; i++)
        {
            adj[i] -= 1;
        }

        int[] perm_inv = new int[node_num];

        typeMethods.perm_inverse3(node_num, perm, ref perm_inv);

        bandwidth = AdjacencyMatrix.adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
            perm, perm_inv);

        Console.WriteLine("");
        Console.WriteLine("  Permuted ADJ bandwidth = " + bandwidth + "");
        //
        //  Permute the nodes according to the permutation vector.
        //
        typeMethods.r8col_permute(dim_num, node_num, perm, ref node_xy);
        //
        //  Permute the node indices in the element array.
        //
        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < element_order; i++)
            {
                int node = element_node[i + j * element_order];
                element_node[i + j * element_order] = perm_inv[node];
            }
        }

        //
        //  Write out the new data.
        //
        typeMethods.r8mat_write(node_rcm_filename, dim_num, node_num, node_xy);

        Console.WriteLine("");
        Console.WriteLine("  Created the node file \"" + node_rcm_filename + "\".");

        typeMethods.i4mat_write(element_rcm_filename, element_order,
            element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Created the element file \"" +
                          element_rcm_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("QUAD_MESH_RCM:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}