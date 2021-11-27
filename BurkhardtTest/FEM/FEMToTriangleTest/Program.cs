﻿using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;
using Triangle = Burkardt.FEM.Triangle;

namespace FEMToTriangleTest;

internal static class Program
{
    private static void Main(string[] args)
    {
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM_TO_TRIANGLE.
        //
        //  Discussion:
        //
        //    FEM_TO_TRIANGLE converts a mesh of triangles from FEM to TRIANGLE format.
        //
        //  Usage:
        //
        //    fem_to_triangle prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'_nodes.txt contains the FEM node coordinates,
        //    * 'prefix'_elements.txt contains the FEM element connectivities.
        //    * 'prefix'.node will contains the TRIANGLE node coordinates,
        //    * 'prefix'.ele will contains the TRIANGLE element connectivities.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            string prefix;

            Console.WriteLine("");
            Console.WriteLine("FEM_TO_TRIANGLE");
            Console.WriteLine("  Convert a 2D mesh from FEM to TRIANGLE format.");
            Console.WriteLine("");
            Console.WriteLine("  Read:");
            Console.WriteLine("  * \"prefix\"_nodes.txt, FEM node coordinates.");
            Console.WriteLine("  * \"prefix\"_elements.txt, FEM element connectivities.");
            Console.WriteLine("");
            Console.WriteLine("  Create:");
            Console.WriteLine("  * \"prefix\".node, TRIANGLE node coordinates.");
            Console.WriteLine("  * \"prefix\".ele, TRIANGLE element connectivities.");
            switch (args.Length)
            {
                //
                //  Get the filename prefix.
                //
                case < 1:
                    Console.WriteLine("");
                    Console.WriteLine("FEM_TO_TRIANGLE:");
                    Console.WriteLine("  Please enter the filename prefix.");

                    prefix = Console.ReadLine();
                    break;
                default:
                    prefix = args[0];
                    break;
            }

            //
            //  Create the filenames.
            //
            string fem_node_filename = prefix + "_nodes.txt";
            string fem_element_filename = prefix + "_elements.txt";
            string triangle_node_filename = prefix + ".node";
            string triangle_element_filename = prefix + ".ele";
            //
            //  Read the node data.
            //
            TableHeader h = typeMethods.r8mat_header_read(fem_node_filename);
            int node_num = h.n;
            int m = h.m;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + fem_node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + m + "");
            Console.WriteLine("  Number of nodes  = " + node_num + "");

            if (m != 2)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_TO_TRIANGLE - Fatal error!");
                Console.WriteLine("  Spatial dimension must be 2.");
                return;
            }

            double[] node_x = typeMethods.r8mat_data_read(fem_node_filename, m, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + fem_node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(m, node_num, node_x, 1, 1, m, 5,
                "  Portion of node coordinate data:");
            //
            //  Read the element data.
            //
            h = typeMethods.i4mat_header_read(fem_element_filename);
            int element_order = h.m;
            int element_num = h.n;


            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + fem_element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Element order = " + element_order + "");
            Console.WriteLine("  Number of elements  = " + element_num + "");

            if (element_order != 3)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_TO_TRIANGLE - Fatal error!");
                Console.WriteLine("  Element order must be 3.");
                return;
            }

            int[] element_node = typeMethods.i4mat_data_read(fem_element_filename, element_order,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + fem_element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
                1, 1, element_order, 10, "  Initial portion of element data:");
            //
            //  Force 1-based indexing.
            //
            Mesh.mesh_base_one(node_num, element_order, element_num, ref element_node);
            //
            //  Write out the TRIANGLE version of the data.
            //
            int element_att_num = 0;
            double[] element_att = new double[1];

            Triangle.triangle_element_write(triangle_element_filename, element_num,
                element_order, element_att_num, element_node, element_att);

            Console.WriteLine("");
            Console.WriteLine("  Created the TRIANGLE element file \"" + triangle_element_filename + "\".");

            int node_att_num = 0;
            double[] node_att = new double[1];
            int node_marker_num = 0;
            int[] node_marker = new int[1];

            Triangle.triangle_node_write(triangle_node_filename, node_num, m,
                node_att_num, node_marker_num, node_x, node_att, node_marker);

            Console.WriteLine("  Created the TRIANGLE node file \"" + triangle_node_filename + "\".");
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("FEM_TO_TRIANGLE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}