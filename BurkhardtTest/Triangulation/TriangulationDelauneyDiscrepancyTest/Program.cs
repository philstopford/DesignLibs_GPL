using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationDelauneyDiscrepancyTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_DELAUNAY_DISCREPANCY.
        //
        //  Discussion:
        //
        //    TRIANGULATION_DELAUNAY_DISCREPANCY measures the amount (possibly zero)
        //    by which a triangulation fails the local Delaunay test.
        //
        //    The local Delaunay considers pairs of neighboring triangles.
        //    The two triangles form a quadrilateral, and their common edge is one
        //    diagonal of that quadrilateral.  The program considers the effect of
        //    replacing that diagonal with the other one.  If the minimum angle
        //    of the original configuration is smaller than in the new configuration,
        //    then the pair of triangles has failed the local Delaunay test.
        //    The amount by which the minimum angle would have increased is the
        //    local discrepancy.
        //
        //    This program searches all pairs of triangles, and records the maximum
        //    discrepancy found.  If this discrepancy is essentially zero, then the
        //    triangulation is a Delaunay triangulation.  Otherwise, it is not a
        //    Delaunay triangulation.  
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
        //    triangulation_delaunay_discrepancy prefix
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
        //    06 May 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle_max = 0;
        int angle_max_triangle = 0;
        double angle_min = 0;
        int angle_min_triangle = 0;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_DELAUNAY_DISCREPANCY:");
        Console.WriteLine("");
        Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangulation dataset of ");
        Console.WriteLine("  TRIANGLE_NUM triangles using 3 or 6 nodes.");
        Console.WriteLine("");
        Console.WriteLine("  Determine the Delaunay discrepancy, that is, the amount");
        Console.WriteLine("  by which the minimum angle in the triangulation could be");
        Console.WriteLine("  changed by a single adjustment of a pair of triangles.");
        Console.WriteLine("");
        Console.WriteLine("  If this discrepancy is negative, ");
        Console.WriteLine("  then the triangulation is not a Delaunay triangulation.");
        Console.WriteLine("");
        Console.WriteLine("  If this discrepancy is 0 or essentially so, ");
        Console.WriteLine("  then the triangulation is a Delaunay triangulation.");
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
            Console.WriteLine("TRIANGULATION_DELAUNAY_DISCREPANCY:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename );
        int dim_num = h.m;
        int node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of nodes NODE_NUM  = " + node_num + "");

        double[] node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, dim_num, 5,
            "  First 5 nodes:");
        //
        //  Read the triangulation data.
        //
        h = typeMethods.i4mat_header_read(element_filename );
        int triangle_order = h.m;
        int triangle_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Triangle order TRIANGLE_ORDER = " + triangle_order + "");
        Console.WriteLine("  Number of triangles TRIANGLE_NUM  = " + triangle_num + "");

        int[] triangle_node = typeMethods.i4mat_data_read(element_filename, triangle_order, triangle_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node,
            1, 1, triangle_order, 10, "  First 10 triangles:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
        //
        //  Create the triangle neighbors.
        //
        int[] triangle_neighbor = NeighborElements.triangulation_neighbor_triangles(triangle_order,
            triangle_num, triangle_node);

        //
        //  Convert to 0-based index.
        //
        for (int j = 0; j < triangle_num; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                triangle_neighbor[i + 3 * j] -= 1;
            }
        }
            
        typeMethods.i4mat_transpose_print_some(3, triangle_num, triangle_neighbor,
            1, 1, 3, 10, "  First 10 triangle neighbors:");
        //
        //  Now we are ready to check.
        //
        double discrepancy = Delauney.triangulation_delaunay_discrepancy_compute(node_num, node_xy,
            triangle_order, triangle_num, triangle_node, triangle_neighbor, ref angle_min,
            ref angle_min_triangle, ref angle_max, ref angle_max_triangle);

        Console.WriteLine("");
        Console.WriteLine("  Discrepancy (degrees) =   " + discrepancy + "");
        Console.WriteLine("  Minimum angle (degrees) = " + angle_min + "");
        Console.WriteLine("  occurred in triangle      " + angle_min_triangle + "");
        Console.WriteLine("  Maximum angle (degrees) = " + angle_max + "");
        Console.WriteLine("  occurred in triangle      " + angle_max_triangle + "");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_DELAUNAY_DISCREPANCY:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}