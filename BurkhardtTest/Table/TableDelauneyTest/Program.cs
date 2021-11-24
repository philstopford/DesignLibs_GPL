using System;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TableDelauneyTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TABLE_DELAUNAY.
        //
        //  Discussion:
        //
        //    TABLE_DELAUNAY computes the Delaunay triangulation of a TABLE dataset.
        //
        //    The dataset is simply a set of points in the plane.
        //
        //    Thus, given a set of points V1, V2, ..., VN, we apply a standard 
        //    Delaunay triangulation.  The Delaunay triangulation is an organization 
        //    of the data into triples, forming a triangulation of the data, with
        //    the property that the circumcircle of each triangle never contains
        //    another data point.  
        //
        //  Usage:
        //
        //    table_delaunay prefix
        //
        //    where:
        //
        //    'prefix' is the common prefix for the node and triangle files:
        //
        //    * prefix_nodes.txt,     the node coordinates (input).
        //    * prefix_elements.txt,  the nodes that make up each triangle (output).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int base_ = 1;
        string prefix;
        int triangle_num = 0;
        int triangle_order = 0;

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("TABLE_DELAUNAY");
        Console.WriteLine("");
        Console.WriteLine("  Read a TABLE dataset of N points in 2 dimensions,");
        Console.WriteLine("  Compute the Delaunay triangulation.");
        Console.WriteLine("  Write an integer TABLE dataset of the triangulation.");
        //
        //  First argument is the file prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TABLE_DELAUNAY:");
            Console.WriteLine("  Please enter the file prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string triangle_filename = prefix + "_elements.txt";
        //
        //  Read the point coordinates.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int node_dim = h.m;
        int node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\"");
        Console.WriteLine("");
        Console.WriteLine("  Node dimension NODE_DIM = " + node_dim + "");
        Console.WriteLine("  Node number    NODE_NUM = " + node_num + "");

        if (node_dim != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("TABLE_DELAUNAY - Fatal error!");
            Console.WriteLine("  The node dimension is not 2.");
            return;
        }

        double[] node_xy = typeMethods.r8mat_data_read(node_filename, node_dim, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data of \"" + node_filename + "\"");

        typeMethods.r8mat_transpose_print_some(node_dim, node_num, node_xy, 1, 1, node_dim, 5,
            "  Initial portion of node data:");
        //
        //  Determine the Delaunay triangulation.
        //
        triangle_order = 3;

        int[] triangle_node = new int[triangle_order * 3 * node_num];
        int[] triangle_neighbor = new int[triangle_order * 3 * node_num];

        Delauney.dtris2(node_num, base_, ref node_xy, ref triangle_num, ref triangle_node,
            ref triangle_neighbor);
        //
        //  Print a portion of the triangulation.
        //
        Console.WriteLine("");
        Console.WriteLine("  Computed the triangulation.");
        Console.WriteLine("  Number of triangles is " + triangle_num + "");

        typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node,
            1, 1, 3, 5, "  Initial portion of triangulation data:");
        //
        //  Write the triangulation to a file.
        //
        typeMethods.i4mat_write(triangle_filename, triangle_order, triangle_num,
            triangle_node);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the triangulation data to \""
                          + triangle_filename + "\".");
        //
        //  Terminate execution.
        //
        Console.WriteLine("");
        Console.WriteLine("TABLE_DELAUNAY:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}