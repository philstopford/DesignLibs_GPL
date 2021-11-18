using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationHistogramTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_HISTOGRAM.
        //
        //  Discussion:
        //
        //    TRIANGULATION_HISTOGRAM reads the definition of a triangulation and
        //    a list of N points, most of which are presumed to lie somewhere within
        //    the triangulation.
        //
        //    The program determines 
        //    * Ai, the area of the Ith triangle;
        //    * Atotal, the area of the triangulation;
        //    * Ni, the number of data points inside the Ith triangle;
        //    * Ntotal, the number of data points inside the triangulation.
        //
        //    It then reports the values of Ai/Atotal and Ni/Ntotal.
        //
        //    This allows an estimation of whether the set of data points represents
        //    a uniform sampling of the triangulation.
        //
        //  Usage:
        //
        //    triangulation_histogram prefix data_filename
        //
        //    where
        //
        //    * prefix is the common filename prefix for a triangulation:
        //      'prefix_nodes.txt' contains the node coordinates;
        //      'prefix_elements.txt' contains lists of node indices forming elements;
        //
        //    * data_filename is the name of the file containing the sample points,
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string data_filename;
        int data_num = 0;
        double[] data_xy;
        int dim_num = 0;
        int element;
        double[] element_area;
        int[] element_hit;
        string element_filename;
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        string node_filename;
        int node_num;
        double[] node_xy;
        string prefix;
        double triangulation_area;
        int triangulation_hit;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_HISTOGRAM:");
        Console.WriteLine("  Compute a histogram for datapoints in a triangulation.");
        //
        //  Command line argument #1 is the filename prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_HISTOGRAM:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Command line argument #2 is the data filename.
        //
        try
        {
            data_filename = args[1];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_HISTOGRAM:");
            Console.WriteLine("  Please enter the data filename.");

            data_filename = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
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
            "  First 5 nodes:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Triangle order ELEMENT_ORDER = " + element_order + "");
        Console.WriteLine("  Number of triangles ELEMENT_NUM  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 10, "  First 10 triangles:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Read the data points.
        //
        h = typeMethods.r8mat_header_read(data_filename);
        dim_num = h.m;
        data_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + data_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of data points DATA_NUM  = " + data_num + "");

        data_xy = typeMethods.r8mat_data_read(data_filename, dim_num, data_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + data_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, data_num, data_xy, 1, 1, dim_num, 5,
            "  First 5 data points:");
        //
        //  Compute the triangle areas.
        //
        element_area = new double[element_num];

        triangulation_area = Triangulation.triangulation_areas(node_num, node_xy, element_order,
            element_num, element_node, element_area);
        //
        //  Count the hits.
        //
        element_hit = new int[element_num];

        triangulation_hit = Hits.triangulation_hits(node_num, node_xy, element_order,
            element_num, element_node, data_num, data_xy, ref element_hit);
        //
        //  Print the histogram.
        //
        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_HISTOGRAM:");
        Console.WriteLine("  Histogram report:");
        Console.WriteLine("");
        Console.WriteLine("  Node data from         \"" + node_filename + "\".");
        Console.WriteLine("  Element data from      \"" + element_filename + "\".");
        Console.WriteLine("  Sample point data from \"" + data_filename + "\".");
        Console.WriteLine("  Number of sample points = " + data_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Triangle        Area          Hits  Area Ratio       Hit Ratio");
        Console.WriteLine("");
        for (element = 0; element < element_num; element++)
        {
            Console.WriteLine("  " + element.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + element_area[element].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + element_hit[element].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + (element_area[element] / triangulation_area).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (element_hit[element]
                                             / (double) triangulation_hit).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  " + "   Total"
                               + "  " + triangulation_area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + triangulation_hit.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + 1.0.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + 1.0.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_HISTOGRAM:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}