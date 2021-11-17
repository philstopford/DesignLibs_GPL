using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationBoundaryNodesTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_BOUNDARY_NODES.
        //
        //  Discussion:
        //
        //    TRIANGULATION_BOUNDARY_NODES outputs boundary nodes of a triangulation.
        //
        //    Each connected segment of the boundary is written to a separate file.
        //
        //  Usage:
        //
        //    triangulation_boundary_nodes  prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txtcontains the element definitions.
        //    * prefix_nodes_boundary_nodes.txt contains a 0 for internal nodes,
        //      and a 1 for boundary nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 December 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] boundarynode;
        string boundarynode_filename;
        int dim_num;
        string element_filename;
        int j;
        int node;
        bool[] node_boundary = null;
        int node_boundary_num;
        string node_filename;
        int node_num;
        double[] node_xy;
        int node2;
        string prefix;
        int[] triangle_node;
        int triangle_num;
        int triangle_order;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_BOUNDARY_NODES");
        Console.WriteLine("  Identify triangulation nodes that lie on the boundary.");
        Console.WriteLine("");
        Console.WriteLine("* Read a node dataset of NODE_NUM points in 2 dimensions;");
        Console.WriteLine("");
        Console.WriteLine("* Read an associated triangulation dataset of ");
        Console.WriteLine("  triangles using 3 or 6 nodes;");
        Console.WriteLine("");
        Console.WriteLine("* Determine which nodes are on the boundary;");
        Console.WriteLine("");
        Console.WriteLine("* Write a file which is 1 for each boundary node.");
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
            Console.WriteLine("TRIANGULATION_BOUNDARY_NODES:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        boundarynode_filename = prefix + "_boundary_nodes.txt";
        //
        //  Read the node coordinates.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        dim_num = h.m;
        node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of nodes NODE_NUM  = " + node_num + "");

        node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, dim_num, 5,
            "  Portion of coordinate data from file:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        triangle_order = h.m;
        triangle_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Triangle order TRIANGLE_ORDER = " + triangle_order + "");
        Console.WriteLine("  Number of triangles TRIANGLE_NUM  = " + triangle_num + "");

        triangle_node = typeMethods.i4mat_data_read(element_filename,
            triangle_order, triangle_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node, 1, 1,
            triangle_order, 5, "  Portion of data read from file:");
        //
        //  Detect and correct 1-based indexing.
        //
        Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
        node_boundary = triangle_order switch
        {
            //
            //  Determine which nodes lie on the boundary.
            //
            3 => Boundary.triangulation_order3_boundary_node(node_num, triangle_num, triangle_node),
            6 => Boundary.triangulation_order6_boundary_node(node_num, triangle_num, triangle_node),
            _ => node_boundary
        };

        //
        //  Print the data
        //
        node_boundary_num = 0;
        for (node = 0; node < node_num; node++)
        {
            switch (node_boundary[node])
            {
                case true:
                    node_boundary_num += 1;
                    break;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of boundary nodes is " + node_boundary_num + "");

        switch (node_boundary_num)
        {
            case < 100:
            {
                Console.WriteLine("");
                Console.WriteLine("  Boundary nodes:");
                Console.WriteLine("");
                Console.WriteLine("     New     Old");
                Console.WriteLine("   Index   Index         X and Y Coordinates");
                Console.WriteLine("");

                node2 = 0;
                for (node = 0; node < node_num; node++)
                {
                    switch (node_boundary[node])
                    {
                        case true:
                            Console.WriteLine("  " + node2.ToString().PadLeft(8)
                                                   + "  " + node.ToString().PadLeft(8)
                                                   + "  " + node_xy[0 + node * 2].ToString().PadLeft(14)
                                                   + "  " + node_xy[1 + node * 2].ToString().PadLeft(14) + "");
                            node2 += 1;
                            break;
                    }
                }

                break;
            }
        }

        //
        //  Write the output file.
        //
        boundarynode = new int[node_num];
        for (j = 0; j < node_num; j++)
        {
            boundarynode[j] = node_boundary[j] switch
            {
                false => 0,
                _ => 1
            };
        }

        typeMethods.i4mat_write(boundarynode_filename, 1, node_num, boundarynode);

        Console.WriteLine("");
        Console.WriteLine("  Created the file \""
                          + boundarynode_filename + "\"");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_BOUNDARY_NODES:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}