﻿using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationCornerTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_CORNER.
        //
        //  Discussion:
        //
        //    TRIANGULATION_CORNER refines a triangulation by doubling.
        //
        //  Usage:
        //
        //    triangulation_corner prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txt contains the element definitions.
        //    * prefix_corner_nodes.txt will contain the revised node coordinates,
        //    * prefix_corner_elements.txt will contain the revised element definitions.
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
        int[] negative_total = new int[4];
        string prefix;
        int t1_to_t2 = 0;
        double[] t3 = new double[2 * 3];

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_CORNER");
        Console.WriteLine("");
        Console.WriteLine("  Read a node file of NODE_NUM point coordinates in 2 dimensions.");
        Console.WriteLine("  Read an associated triangle file of");
        Console.WriteLine("  TRIANGLE_NUM triangles, listing 3 or 6 node indices.");
        Console.WriteLine("");
        Console.WriteLine("  Any triangle which has exactly two sides on the boundary");
        Console.WriteLine("  is a corner triangle.");
        Console.WriteLine("");
        Console.WriteLine("  If there are any corner triangles this program tries to");
        Console.WriteLine("  eliminate them.");
        Console.WriteLine("");
        Console.WriteLine("  The \"repaired\" triangle file is written out.");
        Console.WriteLine("");

        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_CORNER:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();

            //
            //  Create the filenames.
            //
            string node_filename = prefix + "_nodes.txt";
            string element_filename = prefix + "_elements.txt";
            string node_corner_filename = prefix + "_corner_nodes.txt";
            string element_corner_filename = prefix + "_corner_triangles.txt";
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
                Console.WriteLine("TRIANGULATION_CORNER - Fatal error!");
                Console.WriteLine("  Data is not for a 3-node or 6-node triangulation.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Triangle order TRIANGLE_ORDER     = " + triangle_order + "");
            Console.WriteLine("  Number of triangles TRIANGLE_NUM  = " + triangle_num + "");

            int[] triangle_node = typeMethods.i4mat_data_read(element_filename,
                triangle_order, triangle_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num, triangle_node,
                1, 1, triangle_order, 5, "  First 5 triangles:");
            //
            //  Detect and correct 1-based node indexing.
            //
            Mesh.mesh_base_zero(node_num, triangle_order, triangle_num, ref triangle_node);
            //
            //  Create the triangle neighbor array.
            //
            int[] triangle_neighbor = new int[3 * triangle_num];

            switch (triangle_order)
            {
                case 3:
                    Neighbor.triangulation_order3_neighbor_triangles_a(
                        triangle_num, triangle_node, ref triangle_neighbor);
                    break;
                case 6:
                    Neighbor.triangulation_order6_neighbor_triangles_a(
                        triangle_num, triangle_node, ref triangle_neighbor);
                    break;
            }

            //
            //  Examine the triangle neighbor array.
            //
            typeMethods.i4vec_zero(4, ref negative_total);

            int negative;
            int neighbor;
            int triangle;
            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                negative = 0;
                for (neighbor = 0; neighbor < 3; neighbor++)
                {
                    switch (triangle_neighbor[neighbor + triangle * 3])
                    {
                        case < 0:
                            negative += 1;
                            break;
                    }
                }

                negative_total[negative] += 1;
            }

            Console.WriteLine("");
            Console.WriteLine("  Number of boundary sides     Number of triangles");
            Console.WriteLine("");
            int i;
            for (i = 0; i < 4; i++)
            {
                Console.WriteLine("            " + i.ToString().PadLeft(8)
                                                 + "            " + negative_total[i].ToString().PadLeft(8) + "");
            }

            switch (negative_total[3])
            {
                //
                //  Try to patch problems.
                //
                case > 0:
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_CORNER - Fatal error!");
                    Console.WriteLine("  There is at least one triangle with all sides");
                    Console.WriteLine("  on the boundary.");
                    break;
                default:
                {
                    switch (negative_total[2])
                    {
                        case 0:
                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_CORNER:");
                            Console.WriteLine("  No corner triangles were found.");
                            Console.WriteLine("  No corrections need to be made.");
                            break;
                        default:
                        {
                            //
                            //  We need the triangles to be oriented properly.
                            //
                            negative = 0;

                            int node;
                            for (triangle = 0; triangle < triangle_num; triangle++)
                            {
                                int j;
                                for (j = 0; j < 3; j++)
                                {
                                    for (i = 0; i < 2; i++)
                                    {
                                        t3[(t3.Length + i + j * 2) % t3.Length] =
                                            node_xy[( node_xy.Length + i + (triangle_node[( triangle_node.Length + j + triangle * triangle_order) % triangle_node.Length] - 1) * dim_num) % node_xy.Length];
                                    }
                                }

                                double area = typeMethods.triangle_area_2d(t3);

                                switch (area)
                                {
                                    case < 0.0:
                                    {
                                        negative += 1;

                                        node = triangle_node[1 + triangle * triangle_order];
                                        triangle_node[1 + triangle * triangle_order] = triangle_node[2 + triangle * triangle_order];
                                        triangle_node[2 + triangle * triangle_order] = node;

                                        switch (triangle_order)
                                        {
                                            case 6:
                                                node = triangle_node[3 + triangle * triangle_order];
                                                triangle_node[3 + triangle * triangle_order] =
                                                    triangle_node[5 + triangle * triangle_order];
                                                triangle_node[5 + triangle * triangle_order] = node;
                                                break;
                                        }

                                        neighbor = triangle_neighbor[0 + triangle * 3];
                                        triangle_neighbor[0 + triangle * 3] = triangle_neighbor[2 + triangle * 3];
                                        triangle_neighbor[2 + triangle * 3] = neighbor;
                                        break;
                                    }
                                }
                            }

                            switch (negative)
                            {
                                case > 0:
                                    Console.WriteLine("");
                                    Console.WriteLine("TRIANGULATION_CORNER:");
                                    Console.WriteLine("  Reoriented " + negative + " triangles.");
                                    Console.WriteLine("");
                                    break;
                                default:
                                    Console.WriteLine("");
                                    Console.WriteLine("TRIANGULATION_CORNER:");
                                    Console.WriteLine("  Triangles were already properly oriented.");
                                    Console.WriteLine("");
                                    break;
                            }

                            //
                            //  Now consider each triangle that has exactly two boundary sides.
                            //
                            int triangle1;
                            for (triangle1 = 0; triangle1 < triangle_num; triangle1++)
                            {
                                negative = 0;
                                for (neighbor = 0; neighbor < 3; neighbor++)
                                {
                                    switch (triangle_neighbor[neighbor + triangle1 * 3])
                                    {
                                        case < 0:
                                            negative += 1;
                                            break;
                                    }
                                }

                                switch (negative)
                                {
                                    case 2:
                                    {
                                        int triangle2 = -1;

                                        for (neighbor = 0; neighbor < 3; neighbor++)
                                        {
                                            switch (triangle_neighbor[neighbor + triangle1 * 3])
                                            {
                                                case > 0:
                                                    triangle2 = triangle_neighbor[neighbor + triangle1 * 3] - 1;
                                                    t1_to_t2 = neighbor;
                                                    break;
                                            }
                                        }

                                        Console.WriteLine("  Adjusting triangle " + triangle1 + 1
                                                          + " using triangle " + triangle2 + 1 + "");

                                        int t2_to_t1 = -1;
                                        for (neighbor = 0; neighbor < 3; neighbor++)
                                        {
                                            if (triangle_neighbor[neighbor + triangle2 * 3] - 1 == triangle1)
                                            {
                                                t2_to_t1 = neighbor;
                                            }
                                        }

                                        int n1_index = typeMethods.i4_wrap(t1_to_t2 - 1, 0, 2);
                                        node = triangle_node[n1_index + triangle1 * triangle_order];

                                        n1_index = typeMethods.i4_wrap(t1_to_t2 + 1, 0, 2);
                                        int n2_index = typeMethods.i4_wrap(t2_to_t1 - 1, 0, 2);
                                        triangle_node[n1_index + triangle1 * triangle_order]
                                            = triangle_node[n2_index + triangle2 * triangle_order];

                                        n1_index = typeMethods.i4_wrap(t1_to_t2 - 1, 0, 2);
                                        n2_index = typeMethods.i4_wrap(t2_to_t1 + 1, 0, 2);
                                        triangle_node[n2_index + triangle2 * triangle_order] = node;

                                        switch (triangle_order)
                                        {
                                            case 6:
                                                //
                                                //  Adjust coordinates of the new midside node.
                                                //
                                                int j1 = typeMethods.i4_wrap(t1_to_t2 - 1, 0, 2);
                                                int j2 = typeMethods.i4_wrap(t1_to_t2 + 3, 3, 5);
                                                int j3 = typeMethods.i4_wrap(t2_to_t1 - 1, 0, 2);

                                                int node1 = triangle_node[j1 + triangle1 * triangle_order] - 1;
                                                int node2 = triangle_node[j2 + triangle1 * triangle_order] - 1;
                                                int node3 = triangle_node[j3 + triangle2 * triangle_order] - 1;

                                                node_xy[0 + node2 * 2] = 0.5 * (node_xy[0 + node1 * 2] + node_xy[0 + node3 * 2]);
                                                node_xy[1 + node2 * 2] = 0.5 * (node_xy[1 + node1 * 2] + node_xy[1 + node3 * 2]);
                                                //
                                                //  Update the triangle array.
                                                //
                                                j1 = typeMethods.i4_wrap(t1_to_t2 + 4, 3, 5);
                                                j2 = typeMethods.i4_wrap(t1_to_t2 + 3, 3, 5);
                                                j3 = typeMethods.i4_wrap(t2_to_t1 + 4, 3, 5);
                                                int j4 = typeMethods.i4_wrap(t2_to_t1 + 3, 3, 5);

                                                node
                                                    = triangle_node[j1 + triangle1 * triangle_order];
                                                triangle_node[j1 + triangle1 * triangle_order]
                                                    = triangle_node[j2 + triangle1 * triangle_order];
                                                triangle_node[j2 + triangle1 * triangle_order]
                                                    = triangle_node[j3 + triangle2 * triangle_order];
                                                triangle_node[j3 + triangle2 * triangle_order]
                                                    = triangle_node[j4 + triangle2 * triangle_order];
                                                triangle_node[j4 + triangle2 * triangle_order] = node;
                                                break;
                                        }

                                        //
                                        //  Update the neighbor array.
                                        //
                                        n2_index = typeMethods.i4_wrap(t2_to_t1 + 1, 0, 2);
                                        triangle = triangle_neighbor[n2_index + triangle2 * 3] - 1;
                                        triangle_neighbor[n2_index + triangle2 * 3] = triangle1 + 1;
                                        triangle_neighbor[t2_to_t1 + triangle2 * 3] = -1;

                                        n1_index = typeMethods.i4_wrap(t1_to_t2 + 1, 0, 2);
                                        triangle_neighbor[n1_index + triangle1 * 3] = triangle2 + 1;
                                        triangle_neighbor[t1_to_t2 + triangle1 * 3] = triangle + 1;
                                        break;
                                    }
                                }
                            }

                            //
                            //  Write out the corrected triangle file.
                            //
                            typeMethods.i4mat_write(element_corner_filename, triangle_order, triangle_num,
                                triangle_node);

                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_CORNER:");
                            Console.WriteLine("  New triangle file with repaired corners written to");
                            Console.WriteLine("  \"" + element_corner_filename + "\".");
                            //
                            //  Write out the corrected node coordinate file.
                            //  This may only differ for quadratic elements.
                            //
                            typeMethods.r8mat_write(node_corner_filename, dim_num, node_num, node_xy);

                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_CORNER:");
                            Console.WriteLine("  New node coordinate file with adjusted midside nodes");
                            Console.WriteLine("  written to \"" + node_corner_filename + "\".");

                            Console.WriteLine("");
                            Console.WriteLine("TRIANGULATION_CORNER:");
                            Console.WriteLine("  Normal end of execution.");

                            Console.WriteLine("");
                            break;
                        }
                    }

                    break;
                }
            }
        }
    }
}