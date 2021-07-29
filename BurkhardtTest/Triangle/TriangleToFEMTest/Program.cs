using System;
using Burkardt;
using Burkardt.FEM;
using Burkardt.Types;

namespace TriangleToFEMTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGLE_TO_FEM.
            //
            //  Discussion:
            //
            //    The TRIANGLE program creates "node" and "element" files that define
            //    a triangular mesh.  A typical pair of such files might have the names
            //    "suv.node" and "suv.ele".
            //
            //    This program reads this pair of files and creates a pair of FEM files, whose
            //    names might be "suv_nodes.txt" and "suv_elements.txt".
            //
            //  Usage:
            //
            //    triangle_to_fem prefix
            //
            //    where 'prefix' is the common filename prefix so that:
            //
            //    * prefix.node contains the coordinates of the nodes;
            //    * prefix.ele contains the indices of nodes forming each element.
            //    * prefix_nodes.txt will be the FEM node file created by the program;
            //    * prefix_elements.txt will be the FEM element file created by the program.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim_num;
            string element_filename;
            int element_att_num = 0;
            double[] element_att;
            int[] element_node;
            int element_num = 0;
            int element_order = 0;
            string fem_element_filename;
            string fem_node_filename;
            double[] node_att;
            int node_att_num = 0;
            int node_marker_num = 0;
            string node_filename;
            int[] node_marker;
            double[] node_coord;
            int node_dim = 0;
            int node_num = 0;
            string prefix;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_TO_FEM:");
            Console.WriteLine("  Read a pair of NODE and ELE files created by TRIANGLE.");
            Console.WriteLine("  Write a corresponding pair of FEM node and element files.");
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
                Console.WriteLine("  Please enter the filename prefix:");

                prefix = Console.ReadLine();
            }

            //
            //  Create the file names.
            //
            node_filename = prefix + ".node";
            element_filename = prefix + ".ele";
            fem_node_filename = prefix + "_nodes.txt";
            fem_element_filename = prefix + "_elements.txt";

            Console.WriteLine("");
            Console.WriteLine("  Read Node file \"" + node_filename + "\"");
            Console.WriteLine("    and Element file \"" + element_filename + "\".");
            Console.WriteLine("  Create FEM node file \"" + fem_node_filename + "\"");
            Console.WriteLine("    and FEM element file \"" + fem_element_filename + "\".");
            //
            //  Read the TRIANGLE NODE data.
            //
            Node.node_size_read(node_filename, ref node_num, ref node_dim, ref node_att_num,
                ref node_marker_num);

            node_coord = new double[node_dim * node_num];
            node_att = new double[node_att_num * node_num];
            node_marker = new int[node_marker_num * node_num];

            Node.node_data_read(node_filename, node_num, node_dim, node_att_num,
                node_marker_num, ref node_coord, ref node_att, ref node_marker);
            //
            //  Read the TRIANGLE ELE data.
            //
            Element.element_size_read(element_filename, ref element_num, ref element_order,
                ref element_att_num);

            element_node = new int[element_order * element_num];
            element_att = new double[element_att_num * element_num];

            Element.element_data_read(element_filename, element_num, element_order,
                element_att_num, ref element_node, ref element_att);
            //
            //  Write the FEM NODE and ELEMENT data.
            //
            dim_num = 2;

            typeMethods.r8mat_write(fem_node_filename, dim_num, node_num, node_coord);

            typeMethods.i4mat_write(fem_element_filename, element_order, element_num, element_node);

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_TO_FEM:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}