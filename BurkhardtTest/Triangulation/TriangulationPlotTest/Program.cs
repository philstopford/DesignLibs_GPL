using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationPlotTest;

using Plot = Plot;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_PLOT.
        //
        //  Discussion:
        //
        //    TRIANGULATION_PLOT plots a triangulated set of nodes.
        //
        //  Usage:
        //
        //    triangulation_plot prefix node_vis element_vis
        //
        //    where:
        //
        //    'prefix' is the common prefix for the node and element files:
        //
        //    * prefix_nodes.txt,     the node coordinates.
        //    * prefix_elements.txt,  the nodes that make up each element.
        //    * prefix.eps, the plot of the triangulation (output).
        //
        //    'node_vis' indicates the node visibility:
        //
        //    0: do not show the nodes;
        //    1:        show the nodes;
        //    2:        show the nodes, and label them.
        //
        //    'element_vis' indicates the element visibility:
        //
        //    0: do not show the elements;
        //    1:        show the elements;
        //    2:        show the elements, and label them.
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
        string element_filename;
        int[] element_node;
        int element_num;
        int element_order;
        int element_show;
        string node_filename;
        int node_num;
        int node_show;
        double[] node_xy;
        string plot_filename;
        string prefix;

        Console.WriteLine("TRIANGULATION_PLOT");
        Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangulation dataset of ");
        Console.WriteLine("  ELEMENT_NUM element in a triangulation of order");
        Console.WriteLine("  ELEMENT_ORDER = 3, 4, or 6 .");
        Console.WriteLine("");
        Console.WriteLine("  Make an EPS plot of the triangulated data.");
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
            Console.WriteLine("TRIANGULATION_PLOT:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Get the node visibility.
        //
        try
        {
            node_show = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the option for showing the nodes:");
            Console.WriteLine("  0: do not show the nodes;");
            Console.WriteLine("  1:        show the nodes;");
            Console.WriteLine("  2:        show the nodes, and label them.");
            node_show = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get the element visibility.
        //
        try
        {
            element_show = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the option for showing the elements:");
            Console.WriteLine("  0: do not show the elements;");
            Console.WriteLine("  1:        show the elements;");
            Console.WriteLine("  2:        show the elements, and label them.");
            element_show = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        plot_filename = prefix + ".eps";
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
        switch (element_order)
        {
            //
            //  Create the output file.
            //
            case 3:
                Plot.triangulation_order3_plot(plot_filename, node_num, node_xy,
                    element_num, element_node, node_show, element_show);
                break;
            case 4:
                Plot.triangulation_order4_plot(plot_filename, node_num, node_xy,
                    element_num, element_node, node_show, element_show);
                break;
            case 6:
                Plot.triangulation_order6_plot(plot_filename, node_num, node_xy,
                    element_num, element_node, node_show, element_show);
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Created the EPS file \"" + plot_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_PLOT:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}