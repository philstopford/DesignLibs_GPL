using System;
using Burkardt;
using Burkardt.MeasureNS;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangulationQualityTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_QUALITY.
        //
        //  Discussion:
        //
        //    TRIANGULATION_QUALITY determines quality measures for a triangulation.
        //
        //    The code has been modified to 'allow' 6-node triangulations.
        //    However, no effort is made to actually process the midside nodes.
        //    Only information from the vertices is used.
        //
        //    The three quality measures are:
        //
        //      ALPHA_MEASURE
        //      AREA_MEASURE
        //      Q_MEASURE
        //
        //    In each case, the ideal value of the quality measure is 1, and
        //    the worst possible value is 0.
        //
        //    The program also prints out the geometric bandwidth, which is the
        //    bandwidth of the adjacency matrix of the nodes.
        //
        //  Usage:
        //
        //    triangulation_quality prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txt contains the element definitions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local parameters:
        //
        //    Local, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    lists the nodes that make up each element.
        //
        //    Local, int ELEMENT_NUM, the number of elements.
        //
        //    Local, int ELEMENT_ORDER, the order of the elements, either 3 or 6.
        //
        //    Local, int NODE_DIM, the spatial dimension.
        //
        //    Local, int NODE_NUM, the number of nodes.
        //
        //    Local, double NODE_XY[NODE_DIM*NODE_NUM], the point set.
        //
    {
        double alpha_area = 0;
        double alpha_ave = 0;
        double alpha_min = 0;
        double area_ave = 0;
        double area_max = 0;
        double area_min = 0;
        int area_negative = 0;
        double area_ratio = 0;
        double area_std = 0;
        int area_zero = 0;
        string element_filename;
        int[] element_node;
        int element_num;
        int element_order;
        int m = 0;
        int ml = 0;
        int mu = 0;
        int node_dim;
        string node_filename;
        int node_num;
        double[] node_xy;
        string prefix;
        double q_area = 0;
        double q_ave = 0;
        double q_max = 0;
        double q_min = 0;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_QUALITY:");
        Console.WriteLine("  Compute triangulation quality measures.");
        Console.WriteLine("");
        //
        //  The first commandline argument is the filename prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_QUALITY:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
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
        node_dim = h.m;
        node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension NODE_DIM = " + node_dim + "");
        Console.WriteLine("  Number of nodes NODE_NUM  = " + node_num + "");

        node_xy = typeMethods.r8mat_data_read(node_filename, node_dim, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(node_dim, node_num, node_xy, 1, 1, node_dim, 5,
            "  First 5 nodes:");
        //
        //  Read the triangulation data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order      = " + element_order + "");
        Console.WriteLine("  Number of elements = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename,
            element_order, element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 10, "  First 10 elements:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Compute the quality measures.
        //
        Measure.alpha_measure(node_num, node_xy, element_order, element_num,
            element_node, ref alpha_min, ref alpha_ave, ref alpha_area);

        Console.WriteLine("");
        Console.WriteLine("  ALPHA compares the smallest angle against 60 degrees.");
        Console.WriteLine("  Values of ALPHA range from 0 (extremely poor) to 1 (excellent).");
        Console.WriteLine("  (The second figure is the same number in degrees.)");
        Console.WriteLine("");
        Console.WriteLine("  ALPHA_MIN  : minimum over all triangles = " + alpha_min
                                                                         + "  " + 60.0 * alpha_min + "");
        Console.WriteLine("  ALPHA_AVE  : average over all triangles = " + alpha_ave
                                                                         + "  " + 60.0 * alpha_ave + "");
        Console.WriteLine("  ALPHA_AREA : average weighted by area   = " + alpha_area
                                                                         + "  " + 60.0 * alpha_area + "");

        Measure.area_measure(node_num, node_xy, element_order, element_num,
            element_node, ref area_min, ref area_max, ref area_ratio, ref area_ave, ref area_std,
            ref area_negative, ref area_zero);

        Console.WriteLine("");
        Console.WriteLine("  AREA compares the areas of the triangles.");
        Console.WriteLine("  Values of AREA_RATIO range from 0 (extremely poor) to 1 (excellent).");
        Console.WriteLine("");
        Console.WriteLine("  AREA_MIN   : minimum area         = " + area_min + "");
        Console.WriteLine("  AREA_MAX   : maximum area         = " + area_max + "");
        Console.WriteLine("  AREA_RATIO : minimum/maximum area = " + area_ratio + "");
        Console.WriteLine("  AREA_AVE   : average area         = " + area_ave + "");
        Console.WriteLine("  AREA_STD   : standard deviation   = " + area_std + "");
        Console.WriteLine("  AREA_NEG   : area < 0             = " + area_negative + "");
        Console.WriteLine("  AREA_ZERO  : area = 0             = " + area_zero + "");

        Measure.q_measure(node_num, node_xy, element_order, element_num, element_node,
            ref q_min, ref q_max, ref q_ave, ref q_area);

        Console.WriteLine("");
        Console.WriteLine("  Q is the ratio of 2 * inradius to outradius.");
        Console.WriteLine("  Values of Q range from 0 (extremely poor) to 1 (excellent).");
        Console.WriteLine("");
        Console.WriteLine("  Q_MIN  : minimum Q                  = " + q_min + "");
        Console.WriteLine("  Q_MAX  : maximum Q                  = " + q_max + "");
        Console.WriteLine("  Q_AVE  : average Q                  = " + q_ave + "");
        Console.WriteLine("  Q_AREA : average Q weighted by area = " + q_area + "");

        Mesh.bandwidth_mesh(element_order, element_num, element_node, ref ml, ref mu, ref m);

        Console.WriteLine("");
        Console.WriteLine("  The geometric bandwidth          M = " + m + "");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_QUALITY:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}