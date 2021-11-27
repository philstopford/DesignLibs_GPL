﻿using System;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.Types;

namespace FEM1DProjectTest;

internal static class Program
{
    private static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_PROJECT.
//
//  Discussion:
//
//    FEM1D_PROJECT reads files defining a sampling of a (scalar or vector)
//    function of 1 argument, and a list of nodes and elements to use for
//    a finite element representation of the data.
//
//    It computes a set of finite element coefficients to be associated with
//    the given finite element mesh, and writes that information to a file
//    so that an FEM representation is formed by the node, element and value 
//    files.
//
//  Usage:
//
//    fem1d_project sample_prefix fem_prefix
//
//    where 'sample_prefix' is the common prefix for the SAMPLE files:
//
//    * sample_prefix_nodes.txt,  the node coordinates where samples were taken,
//    * sample_prefix_values.txt, the sample values.
//
//    and 'fem_prefix' is the common prefix for the FEM files:
//
//    * fem_prefix_nodes.txt,    the node coordinates.
//    * fem_prefix_elements.txt, the nodes that make up each element;
//    * fem_prefix_values.txt,   the values defined at each node.
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
        string fem_prefix;
        string sample_prefix;

        Console.WriteLine("");
        Console.WriteLine("FEM1D_PROJECT");
        Console.WriteLine("");
        Console.WriteLine("  Read files defining a sampling of a function of 1 argument.");
        Console.WriteLine("  Read files defining a finite element mesh.");
        Console.WriteLine("  Project the sample data onto the mesh, and");
        Console.WriteLine("  write a file of FEM coefficient values.");
//
//  Get the number of command line arguments.
//
        try
        {
            sample_prefix = args[0];
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("Enter the sample file prefix:");
            sample_prefix = Console.ReadLine();
        }

        try
        {
            fem_prefix = args[1];
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("Enter the FEM file prefix:");
            fem_prefix = Console.ReadLine();
        }

//
//  Create the filenames.
//
        string sample_node_filename = sample_prefix + "_nodes.txt";
        string sample_value_filename = sample_prefix + "_values.txt";

        string fem_node_filename = fem_prefix + "_nodes.txt";
        string fem_element_filename = fem_prefix + "_elements.txt";
        string fem_value_filename = fem_prefix + "_values.txt";
//
//  Read the SAMPLE data.
//
        TableHeader h = typeMethods.r8mat_header_read(sample_node_filename);

        int sample_node_dim = h.m;
        int sample_node_num = h.n;

        double[] sample_node_x = typeMethods.r8mat_data_read(sample_node_filename, sample_node_dim,
            sample_node_num);

        Console.WriteLine("");
        Console.WriteLine("  Sample node spatial dimension is " + sample_node_dim + "");
        Console.WriteLine("  Sample node number is            " + sample_node_num + "");

        if (sample_node_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM1D_PROJECT - Fatal error!");
            Console.WriteLine("  Spatial dimension of the sample nodes is not 1.");
            return;
        }

        h = typeMethods.r8mat_header_read(sample_value_filename);

        int sample_value_dim = h.m;
        int sample_value_num = h.n;

        Console.WriteLine("  The SAMPLE value dimension is    " + sample_value_dim + "");
        Console.WriteLine("  the SAMPLE value number is        " + sample_value_num + "");

        if (sample_value_num != sample_node_num)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM1D_PROJECT - Fatal error!");
            Console.WriteLine("  Number of sample values and nodes differ.");
            return;
        }

        double[] sample_value = typeMethods.r8mat_data_read(sample_value_filename, sample_value_dim,
            sample_value_num);
//
//  Read the FEM data.
//
        h = typeMethods.r8mat_header_read(fem_node_filename);

        int fem_node_dim = h.m;
        int fem_node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  The FEM node dimension is        " + fem_node_dim + "");
        Console.WriteLine("  The FEM node number is           " + fem_node_num + "");

        if (fem_node_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM1D_PROJECT - Fatal error!");
            Console.WriteLine("  Spatial dimension of the nodes is not 1.");
            return;
        }

        double[] fem_node_x = typeMethods.r8mat_data_read(fem_node_filename, fem_node_dim, fem_node_num);

        h = typeMethods.i4mat_header_read(fem_element_filename);
        int fem_element_order = h.m;
        int fem_element_num = h.n;

        Console.WriteLine("  The FEM element order is         " + fem_element_order + "");
        Console.WriteLine("  The FEM element number is        " + fem_element_num + "");

        typeMethods.i4mat_data_read(fem_element_filename, fem_element_order,
            fem_element_num);
//
//  Compute the FEM values.
//

        double[] fem_value = FEM_1D_Project.fem1d_approximate(sample_node_num, sample_value_dim, sample_node_x,
            sample_value, fem_node_num, fem_node_x, fem_element_order,
            fem_element_num, sample_value_dim, fem_node_num);
//
//  Write the FEM values.
//
        typeMethods.r8mat_write(fem_value_filename, sample_value_dim, fem_node_num,
            fem_value);

        Console.WriteLine("");
        Console.WriteLine("  Projected FEM values written to \"" + fem_value_filename + "\"");

        Console.WriteLine("");
        Console.WriteLine("FEM1D_PROJECT");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}