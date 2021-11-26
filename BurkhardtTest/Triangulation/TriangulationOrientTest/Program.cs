﻿using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangulationOrientTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_ORIENT.
        //
        //  Discussion:
        //
        //    TRIANGULATION_ORIENT forces a triangulation to have positive orientation.
        //
        //    The user supplies a node file and a triangle file, containing
        //    the coordinates of the nodes, and the indices of the nodes that
        //    make up each triangle.  Either 3-node or 6-node triangles may
        //    be used.
        //
        //    The program reads the data, and for each triangle, determines
        //    whether the triangle has positive orientation.  This essentially
        //    means that the vertices are listed in counter clockwise order.
        //    If the vertices are listed in the wrong order, they are reordered.
        //    The reordered triangle file is written out.
        //
        //    Note that for an order 6 triangulation, the vertices are listed
        //    in the first three positions.
        //
        //  Usage:
        //
        //    triangulation_orient prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txtcontains the element definitions.
        //    * prefix_orient_elements.txt will contain the oriented element definitions.
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
        string prefix;
        double[] t3 = new double[2 * 3];
        int triangle;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_ORIENT");
        Console.WriteLine("  Read a node dataset of NODE_NUM points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangle file of TRIANGLE_NUM");
        Console.WriteLine("  triangles using 3 or 6 nodes.");
        Console.WriteLine(" ");
        Console.WriteLine("  Ensure that every triangle has positive orientation..");
        Console.WriteLine("");
        Console.WriteLine("  Write the reoriented triangle file.");
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
            Console.WriteLine("TRIANGULATION_ORIENT:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string element_orient_filename = prefix + "_orient_elements.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
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
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        int triangle_order = h.m;
        int triangle_num = h.n;

        if (triangle_order != 3 && triangle_order != 6)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_ORIENT - Fatal error!");
            Console.WriteLine("  Data is not for a 3-node or 6-node triangulation.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Triangle order TRIANGLE_ORDER = " + triangle_order + "");
        Console.WriteLine("  Number of triangles TRIANGLE_NUM  = " + triangle_num + "");

        int[] triangle_node = typeMethods.i4mat_data_read(element_filename,
            triangle_order, triangle_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node, 1, 1,
            triangle_order, 5, "  First 5 triangles:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
        //
        //  Compute the area, and reorient if necessary.
        //
        int triangle_negative_num = 0;
        int triangle_zero_num = 0;

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int j;
            int i;
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t3[(i + j * 2 + t3.Length ) % t3.Length] = node_xy[(i + (triangle_node[(j + triangle * triangle_order + triangle_node.Length) % triangle_node.Length] - 1) * dim_num + node_xy.Length) % node_xy.Length];
                }
            }

            double area = typeMethods.triangle_area_2d(t3);

            switch (area)
            {
                case < 0.0:
                {
                    triangle_negative_num += 1;

                    int node = triangle_node[1 + triangle * triangle_order];
                    triangle_node[1 + triangle * triangle_order] = triangle_node[2 + triangle * triangle_order];
                    triangle_node[2 + triangle * triangle_order] = node;

                    switch (triangle_order)
                    {
                        case 6:
                            node = triangle_node[3 + triangle * triangle_order];
                            triangle_node[3 + triangle * triangle_order] = triangle_node[5 + triangle * triangle_order];
                            triangle_node[5 + triangle * triangle_order] = node;
                            break;
                    }

                    //
                    //  As a check, repeat the area calculation.
                    //  Now, we expect to get a positive result.
                    //
                    for (j = 0; j < 3; j++)
                    {
                        for (i = 0; i < 2; i++)
                        {
                            t3[i + j * 2] = node_xy[i + (triangle_node[j + triangle * triangle_order] - 1) * dim_num];
                        }
                    }

                    area = typeMethods.triangle_area_2d(t3);

                    switch (area)
                    {
                        case < 0.0:
                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_ORIENT - Fatal error!");
                            Console.WriteLine("  I thought I fixed the area but I did not!");
                            return;
                    }

                    break;
                }
                case 0.0:
                    triangle_zero_num += 1;
                    break;
            }
        }

        switch (triangle_zero_num)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORIENT - Warning!");
                Console.WriteLine("  You have " + triangle_zero_num + " triangles with");
                Console.WriteLine("  area equal to zero.");
                break;
        }

        switch (triangle_negative_num)
        {
            case > 0:
                typeMethods.i4mat_write(element_orient_filename, triangle_order,
                    triangle_num, triangle_node);

                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORIENT - Warning!");
                Console.WriteLine("  You have " + triangle_negative_num
                                                + " triangles with negative area.");
                Console.WriteLine("");
                Console.WriteLine("  We have reoriented these triangles to have positive");
                Console.WriteLine("  area, and written the new triangle data to");
                Console.WriteLine("  the triangle file \"" + element_orient_filename + "\".");

                typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num,
                    triangle_node, 1, 1, triangle_order, 5, "  First 5 triangles:");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORIENT:");
                Console.WriteLine("  None of your triangles had negative area.");
                Console.WriteLine("  Therefore, no new triangle file was written.");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_ORIENT:");
        Console.WriteLine("  Normal end of execution.");

    }
}