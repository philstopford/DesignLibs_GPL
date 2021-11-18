using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangulationSVGTest;

using Plot = Burkardt.TriangulationNS.Plot;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_SVG.
        //
        //  Discussion:
        //
        //    TRIANGULATION_SVG plots a triangulated set of nodes.
        //
        //  Usage:
        //
        //    triangulation_svg prefix
        //
        //    where:
        //
        //    'prefix' is the common prefix for the node and element files:
        //
        //    * prefix_nodes.txt,     the node coordinates.
        //    * prefix_elements.txt,  the nodes that make up each element.
        //    * prefix.svg,           the plot of the triangulation (output).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        string element_filename;
        int[] element_node;
        int element_num;
        int element_order;
        string node_filename;
        int node_num;
        double[] node_xy;
        string plot_filename;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_SVG");
        Console.WriteLine("  Make an SVG plot of triangulated data.");
        Console.WriteLine("");
        Console.WriteLine("  This program expects two files:");
        Console.WriteLine("  * prefix_nodes.txt,    node coordinates,");
        Console.WriteLine("  * prefix_elements.txt, indices of nodes forming elements,");
        Console.WriteLine("  and creates:");
        Console.WriteLine("  * prefix.svg, an SVG image of the triangulation.");
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
            Console.WriteLine("TRIANGULATION_SVG:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        plot_filename = prefix + ".svg";
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

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, 5, 5,
            "  5 by 5 portion of data read from file:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order ELEMENT_ORDER = " + element_order + "");
        Console.WriteLine("  Number of elements ELEMENT_NUM  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename,
            element_order, element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, 5, 5, "  5 by 5 portion of data read from file:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Create the output file.
        //
        Plot.triangulation_plot(plot_filename, node_num, node_xy, element_order,
            element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_SVG:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}