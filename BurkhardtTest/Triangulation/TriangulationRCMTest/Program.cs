﻿using System;
using Burkardt.Graph;
using Burkardt.MatrixNS;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationRCMTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGULATION_RCM.
            //
            //  Discussion:
            //
            //    TRIANGULATION_RCM applies the RCM reordering to a triangulation.
            //
            //    The user supplies a node file and a triangle file, containing
            //    the coordinates of the nodes, and the indices of the nodes that
            //    make up each triangle.  Either 3-node or 6-node triangles may
            //    be used.
            //
            //    The program reads the data, computes the adjacency information,
            //    carries out the RCM algorithm to get the permutation, applies
            //    the permutation to the nodes and triangles, and writes out
            //    new node and triangle files that correspond to the RCM permutation.
            //
            //    Note that node data is normally two dimensional, that is,
            //    each node has an X and Y coordinate.  In some applications, it
            //    may be desirable to specify more information.  This program
            //    will accept node data that includes DIM_NUM entries on each line,
            //    as long as DIM_NUM is the same for each entry.  
            //
            //  Usage:
            //
            //    triangulation_rcm prefix
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
            //    28 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] adj = null;
            int adj_num = 0;
            int[] adj_row = null;
            int bandwidth;
            int dim_num;
            int i;
            string node_filename;
            string element_filename;
            int j;
            int node;
            int node_num;
            double[] node_xy = null;
            string node_rcm_filename;
            string element_rcm_filename;
            int[] perm;
            int[] perm_inv;
            string prefix;
            int triangle_num;
            int triangle_order;
            int[] triangle_neighbor;
            int[] triangle_node;

            Console.WriteLine("");

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_RCM");
            Console.WriteLine("  C++ version:");
            Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
            Console.WriteLine("  Read an associated triangulation dataset of TRIANGLE_NUM");
            Console.WriteLine("  triangles using 3 or 6 nodes.");
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
                Console.WriteLine("TRIANGULATION_RCM");
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
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("  Number of nodes NODE_NUM =  " + node_num + "");

            node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1,
                dim_num, 5, "  Coordinates of first 5 nodes:");
            //
            //  Read the element data.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            triangle_order = h.m;
            triangle_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Triangle order TRIANGLE_ORDER = " + triangle_order + "");
            Console.WriteLine("  Number of triangles TRIANGLE_NUM = " + triangle_num + "");

            if (triangle_order != 3 && triangle_order != 6)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_RCM - Fatal error!");
                Console.WriteLine("  Data is not for a 3-node or 6-node triangulation.");
                return;
            }

            triangle_node = typeMethods.i4mat_data_read(element_filename,
                triangle_order, triangle_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node, 1, 1,
                triangle_order, 10, "  First 5 triangles::");
            //
            //  Detect and correct 1-based node indexing.
            //
            Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
            //
            //  Create the element neighbor array.
            //
            triangle_neighbor = NeighborElements.triangulation_neighbor_triangles(triangle_order,
                triangle_num, triangle_node);
            //
            //  Set up the information needed for the RCM computation.
            //
            if (triangle_order == 3)
            {
                //
                //  Count the number of adjacencies, and set up the ADJ_ROW 
                //  adjacency pointer array.
                //
                adj_row = new int [node_num + 1];

                adj_num = Adjacency.triangulation_order3_adj_count(node_num, triangle_num,
                    triangle_node, triangle_neighbor, adj_row);
                //
                //  Set up the ADJ adjacency array.
                //
                adj = Adjacency.triangulation_order3_adj_set(node_num, triangle_num, triangle_node,
                    triangle_neighbor, adj_num, adj_row);
            }
            else if (triangle_order == 6)
            {
                //
                //  Count the number of adjacencies, and set up the ADJ_ROW 
                //  adjacency pointer array.
                //
                adj_row = new int [node_num + 1];

                adj_num = Adjacency.triangulation_order6_adj_count(node_num, triangle_num,
                    triangle_node, ref triangle_neighbor, ref adj_row);
                //
                //  Set up the ADJ adjacency array.
                //
                adj = Adjacency.triangulation_order6_adj_set(node_num, triangle_num, triangle_node,
                    triangle_neighbor, adj_num, adj_row);
            }

            bandwidth = AdjacencyMatrix.adj_bandwidth(node_num, adj_num, adj_row, adj);

            Console.WriteLine("");
            Console.WriteLine("  ADJ bandwidth = " + bandwidth + "");
            //
            //  Compute the RCM permutation.
            //
            perm = new int[node_num];

            perm = GenRCM.genrcm(node_num, adj_num, adj_row, adj);

            perm_inv = new int[node_num];

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
            //  Permute the node indices in the triangle array.
            //
            for (j = 0; j < triangle_num; j++)
            {
                for (i = 0; i < triangle_order; i++)
                {
                    node = triangle_node[i + j * triangle_order];
                    triangle_node[i + j * triangle_order] = perm_inv[node - 1];
                }
            }

            //
            //  Write out the new data.
            //
            typeMethods.r8mat_write(node_rcm_filename, dim_num, node_num, node_xy);

            Console.WriteLine("");
            Console.WriteLine("  Created the node file \"" + node_rcm_filename + "\".");

            typeMethods.i4mat_write(element_rcm_filename, triangle_order,
                triangle_num, triangle_node);

            Console.WriteLine("");
            Console.WriteLine("  Created the triangulation file \"" +
                              element_rcm_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_RCM:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}