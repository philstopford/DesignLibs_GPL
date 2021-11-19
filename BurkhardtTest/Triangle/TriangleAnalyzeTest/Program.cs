using System;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangleAnalyzeTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_ANALYZE.
        //
        //  Discussion:
        //
        //    TRIANGLE_ANALYZE reports properties of a triangle.
        //
        //  Usage:
        //
        //    triangle_analyze filename
        //
        //    where "filename" is a file containing the coordinates of the vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] angles = new double[3];
        double area;
        double[] centroid;
        double[] circum_center = new double[2];
        double circum_radius = 0;
        int dim_num;
        double[] edge_length;
        bool flag = false;
        int i;
        double[] in_center = new double[2];
        double in_radius = 0;
        string node_filename;
        int node_num;
        double[] node_xy;
        int orientation;
        double[] ortho_center = new double[2];
        double r8_pi = 3.141592653589793;
        double quality;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_ANALYZE:");
        Console.WriteLine("  Determine properties of a triangle.");

        try
        {
            node_filename = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Please enter the name of the node coordinate file.");
            node_filename = Console.ReadLine();
        }

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
        Console.WriteLine("  Number of points NODE_NUM = " + node_num + "");

        if (dim_num != 2)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_ANALYZE - Fatal error!");
            Console.WriteLine("  Dataset must have spatial dimension 2.");
            return;
        }

        if (node_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_ANALYZE - Fatal error!");
            Console.WriteLine("  Dataset must have 3 nodes.");
            return;
        }

        node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print(dim_num, node_num, node_xy, "  Node coordinates:");
        //
        //  ANGLES
        //
        typeMethods.triangle_angles_2d(node_xy, angles);

        typeMethods.r8vec_print(3, angles, "  ANGLES (radians):");

        for (i = 0; i < 3; i++)
        {
            angles[i] = angles[i] * 180.0 / r8_pi;
        }

        typeMethods.r8vec_print(3, angles, "  ANGLES (degrees):");
        //
        //  AREA
        //
        area = typeMethods.triangle_area_2d(node_xy);

        Console.WriteLine("");
        Console.WriteLine("  AREA: " + area + "");
        //
        //  CENTROID
        //
        centroid = typeMethods.triangle_centroid_2d(node_xy);

        Console.WriteLine("");
        Console.WriteLine("  CENTROID: "
                          + "  " + centroid[0].ToString().PadLeft(14)
                          + "  " + centroid[1].ToString().PadLeft(14) + "");
        //
        //  CIRCUM_CIRCLE
        //
        typeMethods.triangle_circumcircle_2d(node_xy, ref circum_radius, ref circum_center);

        Console.WriteLine("");
        Console.WriteLine("  CIRCUM_RADIUS: " + circum_radius + "");
        Console.WriteLine("  CIRCUM_CENTER: "
                          + "  " + circum_center[0].ToString().PadLeft(14)
                          + "  " + circum_center[1].ToString().PadLeft(14) + "");
        //
        //  EDGE LENGTHS
        //
        edge_length = typeMethods.triangle_edge_length_2d(node_xy);

        typeMethods.r8vec_print(3, edge_length, "  EDGE_LENGTHS:");
        //
        //  IN_CIRCLE
        //
        typeMethods.triangle_incircle_2d(node_xy, ref in_center, ref in_radius);

        Console.WriteLine("");
        Console.WriteLine("  IN_RADIUS: " + in_radius + "");
        Console.WriteLine("  IN_CENTER: "
                          + "  " + in_center[0].ToString().PadLeft(14)
                          + "  " + in_center[1].ToString().PadLeft(14) + "");
        //
        //  ORIENTATION
        //
        orientation = typeMethods.triangle_orientation_2d(node_xy);

        Console.WriteLine("");
        switch (orientation)
        {
            case 0:
                Console.WriteLine("  ORIENTATION: CounterClockwise.");
                break;
            case 1:
                Console.WriteLine("  ORIENTATION: Clockwise.");
                break;
            case 2:
                Console.WriteLine("  ORIENTATION: Degenerate Distinct Colinear Points.");
                break;
            case 3:
                Console.WriteLine("  ORIENTATION: Degenerate, at least two points identical.");
                break;
        }

        //
        //  ORTHO_CENTER
        //
        typeMethods.triangle_orthocenter_2d(node_xy, ortho_center, ref flag);

        switch (flag)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  ORTHO_CENTER: Could not be computed.");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  ORTHO_CENTER: "
                                  + "  " + ortho_center[0].ToString().PadLeft(14)
                                  + "  " + ortho_center[1].ToString().PadLeft(14) + "");
                break;
        }

        //
        //  QUALITY
        //
        quality = typeMethods.triangle_quality_2d(node_xy);

        Console.WriteLine("");
        Console.WriteLine("  QUALITY: " + quality + "");

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_ANALYZE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}