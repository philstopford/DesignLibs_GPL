using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationTriangleNeighborsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_TRIANGLE_NEIGHBORS.
        //
        //  Discussion:
        //
        //    TRIANGULATION_TRIANGLE_NEIGHBORS determines the neighbor triangles 
        //    of each triangle in a triangulation.
        //
        //    The user supplies a node file and a triangle file, containing
        //    the coordinates of the nodes, and the indices of the nodes that
        //    make up each triangle.  Either 3-node or 6-node triangles may
        //    be used.
        //
        //    The program reads the node and triangle data, computes the triangle
        //    neighbor information, and writes it to a file.
        //
        //  Usage:
        //
        //    triangulation_triangle_neighbors prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txt contains the element definitions.
        //    * prefix_element_neighbors.txt will contain the triangle neighbors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        string neighbor_filename;
        string node_filename;
        int node_num;
        double[] node_xy;
        string prefix;
        string element_filename;
        int[] triangle_neighbor;
        int[] triangle_node;
        int triangle_num;
        int triangle_order;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_TRIANGLE_NEIGHBORS.");
        Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangulation dataset of ");
        Console.WriteLine("  TRIANGLE_NUM triangles using 3 or 6 nodes.");
        Console.WriteLine("");
        Console.WriteLine("  For each triangle, determine the indices of the");
        Console.WriteLine("  triangles opposite vertices 1, 2 and 3.");
        Console.WriteLine("");
        Console.WriteLine("  Write this triangle neighbor data to files.");
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
            Console.WriteLine("TRIANGULATION_TRIANGLE_NEIGHBORS:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        neighbor_filename = prefix + "_neighbors.txt";
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
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
        //
        //  Create the triangle neighbors.
        //
        triangle_neighbor = NeighborElements.triangulation_neighbor_triangles(triangle_order,
            triangle_num, triangle_node);
        //
        //  Write the output file.
        //
        typeMethods.i4mat_write(neighbor_filename, 3, triangle_num, triangle_neighbor);

        Console.WriteLine("");
        Console.WriteLine("  Created the triangle neighbor file \""
                          + neighbor_filename + "\"");


        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_TRIANGLE_NEIGHBORS.:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}