﻿using System;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace FEM2DSampleTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FEM2D_SAMPLE.
            //
            //  Discussion:
            //
            //    FEM2D_SAMPLE reads files defining a 2D FEM representation of data,
            //    and a set of sample points, and writes out a file containing the 
            //    value of the finite element function at the sample points.
            //
            //  Usage:
            //
            //    fem2d_sample fem_prefix sample_prefix
            //
            //    where 'fem_prefix' is the common prefix for the FEM files:
            //
            //    * fem_prefix_nodes.txt,    the node coordinates.
            //    * fem_prefix_elements.txt, the nodes that make up each element;
            //    * fem_prefix_values.txt,   the values defined at each node.
            //
            //    and 'sample_prefix' is the common prefix for the SAMPLE files.
            //    (the node file is input, and the values file is created by the program.)
            //
            //    * sample_prefix_nodes.txt,  the node coordinates where samples are desired.
            //    * sample_prefix_values.txt, the values computed at each sample node.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string fem_element_filename;
            int[] fem_element_neighbor = new int[1];
            int[] fem_element_node;
            int fem_element_num;
            int fem_element_order;
            int fem_node_dim;
            string fem_node_filename;
            int fem_node_num;
            double[] fem_node_xy;
            string fem_prefix;
            double[] fem_value;
            int fem_value_dim;
            string fem_value_filename;
            int fem_value_num;
            int sample_node_dim;
            string sample_node_filename;
            int sample_node_num;
            double[] sample_node_xy;
            string sample_prefix;
            double[] sample_value;
            int sample_value_dim;
            string sample_value_filename;
            int sample_value_num;

            Console.WriteLine("");
            Console.WriteLine("FEM2D_SAMPLE");
            Console.WriteLine("");
            Console.WriteLine("  Read files defining an FEM function of 2 arguments.");
            Console.WriteLine("  Read a file of sample arguments.");
            Console.WriteLine("  Write a file of function values at the arguments.");
            //
            //  Get the number of command line arguments.
            //
            try
            {
                fem_prefix = args[0];
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("Enter the FEM file prefix:");
                fem_prefix = Console.ReadLine();
            }

            try
            {
                sample_prefix = args[1];
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("Enter the sample file prefix:");
                sample_prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            fem_node_filename = fem_prefix + "_nodes.txt";
            fem_element_filename = fem_prefix + "_elements.txt";
            fem_value_filename = fem_prefix + "_values.txt";

            sample_node_filename = sample_prefix + "_nodes.txt";
            sample_value_filename = sample_prefix + "_values.txt";
            //
            //  Read the FEM data.
            //
            TableHeader h = typeMethods.r8mat_header_read(fem_node_filename);
            fem_node_dim = h.m;
            fem_node_num = h.n;

            fem_node_xy = typeMethods.r8mat_data_read(fem_node_filename, fem_node_dim, fem_node_num);

            Console.WriteLine("");
            Console.WriteLine("  The FEM node dimension is        " + fem_node_dim + "");
            Console.WriteLine("  The FEM node number is           " + fem_node_num + "");

            if (fem_node_dim != 2)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_SAMPLE - Fatal error!");
                Console.WriteLine("  Spatial dimension of the nodes is not 2.");
                return;
            }

            h = typeMethods.i4mat_header_read(fem_element_filename);
            fem_element_order = h.m;
            fem_element_num = h.n;

            fem_element_node = typeMethods.i4mat_data_read(fem_element_filename, fem_element_order,
                fem_element_num);

            Console.WriteLine("  The FEM element order is         " + fem_element_order + "");
            Console.WriteLine("  The FEM element number is        " + fem_element_num + "");

            h = typeMethods.r8mat_header_read(fem_value_filename);
            fem_value_dim = h.m;
            fem_value_num = h.n;

            Console.WriteLine("  The FEM value order is           " + fem_value_dim + "");
            Console.WriteLine("  the FEM value number is          " + fem_value_num + "");

            if (fem_value_num != fem_node_num)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_SAMPLE - Fatal error!");
                Console.WriteLine("  Number of values and nodes differ.");
                return;
            }

            fem_value = typeMethods.r8mat_data_read(fem_value_filename, fem_value_dim, fem_value_num);
            //
            //  Create the element neighbor array.
            //
            if (fem_element_order == 3)
            {
                Neighbor.triangulation_order3_neighbor_triangles(
                    fem_element_num, fem_element_node, ref fem_element_neighbor);
            }
            else if (fem_element_order == 6)
            {
                Neighbor.triangulation_order6_neighbor_triangles(
                    fem_element_num, fem_element_node, ref fem_element_neighbor);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_SAMPLE - Fatal error!");
                Console.WriteLine("  The element order must be 3 or 6.");
                Console.WriteLine("  But this data has element order = " + fem_element_order + "");
                return;
            }

            Console.WriteLine("  The element neighbor array has been computed.");
            //
            //  Read the SAMPLE node data.
            //
            h = typeMethods.r8mat_header_read(sample_node_filename);
            sample_node_dim = h.m;
            sample_node_num = h.n;

            sample_node_xy = typeMethods.r8mat_data_read(sample_node_filename, sample_node_dim,
                sample_node_num);

            Console.WriteLine("");
            Console.WriteLine("  Sample node spatial dimension is " + sample_node_dim + "");
            Console.WriteLine("  Sample node number is            " + sample_node_num + "");

            if (sample_node_dim != 2)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_SAMPLE - Fatal error!");
                Console.WriteLine("  Spatial dimension of the sample nodes is not 2.");
                return;
            }

            //
            //  Compute the SAMPLE values.
            //
            sample_value_dim = fem_value_dim;
            sample_value_num = sample_node_num;

            sample_value = FEM_2D_Evaluate.fem2d_evaluate(fem_node_num, fem_node_xy, fem_element_order,
                fem_element_num, fem_element_node, fem_element_neighbor, fem_value_dim,
                fem_value, sample_node_num, sample_node_xy);
            //
            //  Write the sample values.
            //
            typeMethods.r8mat_write(sample_value_filename, sample_value_dim, sample_value_num,
                sample_value);

            Console.WriteLine("");
            Console.WriteLine("  Interpolated FEM data written to \"" + sample_value_filename + "\"");

            Console.WriteLine("");
            Console.WriteLine("FEM2D_SAMPLE");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}