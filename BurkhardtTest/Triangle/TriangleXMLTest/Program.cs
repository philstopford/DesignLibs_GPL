using System;
using Burkardt.IO;
using Burkardt.Types;

namespace TriangleXMLTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_TO_XML.
        //
        //  Discussion:
        //
        //    The TRIANGLE program creates "node" and "element" files that define
        //    a triangular mesh.  A typical pair of such files might have the names
        //    "suv.node" and "suv.ele".
        //
        //    This program reads this pair of files and creates an equivalent XML file
        //    whose name might be "suv.xml".
        //
        //  Usage:
        //
        //    triangle_to_xml prefix
        //
        //    where 'prefix' is the common filename prefix so that:
        //
        //    * prefix.node contains the coordinates of the nodes;
        //    * prefix.ele contains the indices of nodes forming each element.
        //    * prefix.xml will be the XML file created by the program.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_att_num = 0;
        double[] element_att;
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        int m = 0;
        double[] node_att;
        int node_att_num = 0;
        int node_marker_num = 0;
        int[] node_marker;
        double[] node_x;
        int node_num = 0;
        string prefix;
        string triangle_element_filename;
        string triangle_node_filename;
        string xml_filename;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_TO_XML:");
        Console.WriteLine("  Read a pair of NODE and ELE files created by TRIANGLE.");
        Console.WriteLine("  Write a corresponding XML mesh file.");
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
        triangle_node_filename = prefix + ".node";
        triangle_element_filename = prefix + ".ele";
        xml_filename = prefix + ".xml";

        Console.WriteLine("");
        Console.WriteLine("  Read:");
        Console.WriteLine("  * TRIANGLE node file \"" + triangle_node_filename + "\"");
        Console.WriteLine("  * TRIANGLE element file \"" + triangle_element_filename + "\".");
        Console.WriteLine("  Create:");
        Console.WriteLine("  * XML file \"" + xml_filename + "\".");
        //
        //  Read the TRIANGLE NODE data.
        //
        typeMethods.triangle_node_size_read(triangle_node_filename, ref node_num, ref m, ref node_att_num,
            ref node_marker_num);

        node_x = new double[m * node_num];
        node_att = new double[node_att_num * node_num];
        node_marker = new int[node_marker_num * node_num];

        typeMethods.triangle_node_data_read(triangle_node_filename, node_num, m, node_att_num,
            node_marker_num, ref node_x, ref node_att, ref node_marker);
        //
        //  Read the TRIANGLE ELE data.
        //
        typeMethods.triangle_element_size_read(triangle_element_filename, ref element_num, ref element_order,
            ref element_att_num);

        element_node = new int[element_order * element_num];
        element_att = new double[element_att_num * element_num];

        typeMethods.triangle_element_data_read(triangle_element_filename, element_num, element_order,
            element_att_num, ref element_node, ref element_att);
        //
        //  Write the XML file.
        //
        XML.xml_mesh2d_write(xml_filename, m, node_num, node_x, element_order,
            element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Created XML file \"" + xml_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_TO_XML:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}