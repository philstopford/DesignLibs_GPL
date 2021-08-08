using System;
using Burkardt;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshRCMTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TET_MESH_RCM.
            //
            //  Discussion:
            //
            //    TET_MESH_RCM applies the RCM reordering to a tet mesh.
            //
            //    The user supplies a node file and a tetrahedron file, containing
            //    the coordinates of the nodes, and the indices of the nodes that
            //    make up each tetrahedron.  Either 4-node or 10-node tetrahedrons may
            //    be used.
            //
            //    The program reads the data, computes the adjacency information,
            //    carries out the RCM algorithm to get the permutation, applies
            //    the permutation to the nodes and tetrahedrons, and writes out
            //    new node and tetrahedron files that correspond to the RCM permutation.
            //
            //    Note that node data is normally three dimensional, that is,
            //    each node has an X, Y and Z coordinate.  In some applications, it
            //    may be desirable to specify more information.  This program
            //    will accept node data that includes DIM_NUM entries on each line,
            //    as long as DIM_NUM is the same for each entry.  
            //
            //    Thanks to Xingxing Zhang for pointing out some problems with a
            //    previous version of this program, 10 May 2011.
            //
            //  Usage:
            //
            //    tet_mesh_rcm prefix
            //
            //    where prefix is the common file prefix:
            //
            //    * prefix_nodes.txt,           the node coordinates (input);
            //    * prefix_elements.txt,        the element definitions (input).
            //    * prefix_rcm_nodes.txt,       the new node coordinates (output);
            //    * prefix_rcm_elements.txt,    the new element definitions (output).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] adj = new int[1];
            int adj_num = 0;
            int[] adj_row = new int[1];
            int bandwidth;
            int base_internal;
            int base_user;
            bool debug = false;
            int dim_num;
            int[] element_node;
            int element_num;
            int element_order;
            int i;
            string node_filename;
            string element_filename;
            int j;
            int node;
            int node_num;
            double[] node_xyz;
            string node_rcm_filename;
            string element_rcm_filename;
            int[] perm;
            int[] perm_inv;
            string prefix;

            Console.WriteLine("");

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_RCM");
            Console.WriteLine("  Read a node dataset of NODE_NUM points in 3 dimensions.");
            Console.WriteLine("  Read an associated tet mesh dataset of TETRA_NUM");
            Console.WriteLine("  tetrahedrons using 4 or 10 nodes.");
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
                Console.WriteLine("TET_MESH_RCM:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            node_rcm_filename = prefix + "_rcm_nodes.txt";
            element_rcm_filename = prefix + "_rcm_elements.txt";
            //
            //  Read the node data.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_filename);
            dim_num = h.m;
            node_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension DIM_NUM =  " + dim_num + "");
            Console.WriteLine("  Number of points NODE_NUM  = " + node_num + "");

            node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xyz, 1, 1, dim_num, 5,
                "  Coordinates of first 5 nodes:");
            //
            //  Read the tet mesh data.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

            if (element_order != 4 && element_order != 10)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_RCM - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 4 or order 10.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron order = " + element_order + "");
            Console.WriteLine("  Number of tetras  = " + element_num + "");

            element_node = typeMethods.i4mat_data_read(element_filename, element_order,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(element_order, element_num,
                element_node, 1, 1, element_order, 5, "  First 5 tetrahedrons:");
            //
            //  If the element information is 1-based, make it 0-based.
            //
            base_user = TetMesh.tet_mesh_base_zero(node_num, element_order,
                element_num, ref element_node);

            if (base_user != 0 && base_user != 1)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_RCM - Fatal error!");
                Console.WriteLine("  The input data does not seem to be 0-based or 1-based.");
            }

            //
            //  Following code depends on the element order.
            //
            if (element_order == 4)
            {
                //
                //  Count the number of adjacencies.
                //  Set up the ADJ_ROW adjacency pointer array.
                //
                adj_row = new int[node_num + 1];

                TetMesh.tet_mesh_order4_adj_count(node_num, element_num, element_node,
                    ref adj_num, ref adj_row);

                if (debug || node_num < 10)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  ADJ_NUM = " + adj_num + "");

                    typeMethods.i4vec_print(node_num + 1, adj_row, "  ADJ_ROW:");
                }

                //
                //  Set up the ADJ adjacency array.
                //
                adj = TetMesh.tet_mesh_order4_adj_set(node_num, element_num, element_node,
                    adj_num, adj_row);

                if (node_num < 10)
                {
                    AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ");
                }
            }
            else if (element_order == 10)
            {
                //
                //  Count the number of adjacencies.
                //  Set up the ADJ_ROW adjacency pointer array.
                //
                adj_row = new int[node_num + 1];

                TetMesh.tet_mesh_order10_adj_count(node_num, element_num, element_node,
                    ref adj_num, ref adj_row);

                if (debug || node_num < 10)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  ADJ_NUM = " + adj_num + "");

                    typeMethods.i4vec_print(node_num + 1, adj_row, "  ADJ_ROW:");
                }

                //
                //  Set up the ADJ adjacency array.
                //
                adj = TetMesh.tet_mesh_order10_adj_set(node_num, element_num, element_node,
                    ref adj_num, ref adj_row);

                if (node_num < 10)
                {
                    AdjacencyMatrix.adj_print(node_num, adj_num, adj_row, adj, "  ADJ");
                }
            }

            //
            //  Compute the bandwidth.
            //
            bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

            Console.WriteLine("");
            Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");
            //
            //  GENRCM computes the RCM permutation.
            //
            perm = GenRCM.genrcm(node_num, adj_num, adj_row, adj);
            //
            //  Compute the inverse permutation.
            //
            perm_inv = typeMethods.perm_inverse3(node_num, perm);

            if (node_num < 10)
            {
                Console.WriteLine("");
                Console.WriteLine("         I   PERM[I] INVERSE[I]");
                Console.WriteLine("");
                for (i = 0; i < node_num; i++)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + perm[i].ToString().PadLeft(8)
                                           + "  " + perm_inv[i].ToString().PadLeft(8) + "");
                }
            }

            //
            //  Compute the bandwidth of the permuted array.
            //
            bandwidth = AdjacencyMatrix.adj_perm_bandwidth(node_num, adj_num, adj_row, adj,
                perm, perm_inv);

            Console.WriteLine("");
            Console.WriteLine("  ADJ bandwidth after RCM permutation = " + bandwidth + "");
            //
            //  Permute the nodes in NODE_XYZ.
            //
            base_internal = 0;
            typeMethods.r8col_permute(dim_num, node_num, perm, base_internal, node_xyz);
            //
            //  Permute the node indices in ELEMENT_NODE.
            //
            for (j = 0; j < element_num; j++)
            {
                for (i = 0; i < element_order; i++)
                {
                    node = element_node[i + j * element_order];
                    element_node[i + j * element_order] = perm_inv[node];
                }
            }

            //
            //  If the user base was 1, restore it!
            //
            if (base_user == 1)
            {
                for (i = 0; i < element_num; i++)
                {
                    for (j = 0; j < element_order; j++)
                    {
                        element_node[i + j * element_order] = element_node[i + j * element_order] + 1;
                    }
                }

                Console.WriteLine("");
                Console.WriteLine("  Output files will use the same 1-based ordering used by the input.");
            }

            //
            //  Write the node and element data.
            //
            typeMethods.r8mat_write(node_rcm_filename, dim_num, node_num, node_xyz);

            Console.WriteLine("");
            Console.WriteLine("  Created the file \"" + node_rcm_filename + "\".");

            typeMethods.i4mat_write(element_rcm_filename, element_order, element_num, element_node);

            Console.WriteLine("  Created the file \"" + element_rcm_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_RCM:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}