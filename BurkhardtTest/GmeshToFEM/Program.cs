using System;
using Burkardt.GMesh;
using Burkardt.Types;

namespace GmeshToFEM;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for GMSH_TO_FEM.
        //
        //  Discussion:
        //
        //    GMSH_TO_FEM converts mesh data from GMSH to FEM format.
        //
        //  Usage:
        //
        //    gmsh_to_fem prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'.msh contains the GMSH mesh file.
        //    * 'prefix'_nodes.txt will contain the node coordinates.
        //    * 'prefix'_elements.txt will contain the element node connectivity.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        string fem_element_filename;
        string fem_node_filename;
        string gmsh_filename;
        int m = 0;
        int node_num = 0;
        double[] node_x;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("GMSH_TO_FEM");
        Console.WriteLine("  Read a mesh description created by GMSH:");
        Console.WriteLine("  * \"prefix\".msh, the GMSH mesh file.");
        Console.WriteLine("  Write out two simple FEM files listing nodes and elements.");
        Console.WriteLine("  * \"prefix\"_nodes.txt, node coordinates.");
        Console.WriteLine("  * \"prefix\"_elements.txt, element connectivity.");
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
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        gmsh_filename = prefix + ".msh";
        fem_node_filename = prefix + "_nodes.txt";
        fem_element_filename = prefix + "_elements.txt";
        //
        //  Read GMSH sizes.
        //
        IO.gmsh_size_read(gmsh_filename, ref node_num, ref m, ref element_num, ref element_order);
        //
        //  Report sizes.
        //
        Console.WriteLine("");
        Console.WriteLine("  Size information from GMSH:");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Number of nodes NODE_NUM = " + node_num + "");
        Console.WriteLine("  Number of elements ELEMENT_NUM = " + element_num + "");
        Console.WriteLine("  Element order ELEMENT_ORDER = " + element_order + "");
        //
        //  Allocate memory.
        //
        node_x = new double[m * node_num];
        element_node = new int[element_order * element_num];
        //
        //  Read GMSH data.
        //
        IO.gmsh_data_read(gmsh_filename, m, node_num, node_x, element_order,
            element_num, element_node);
        //
        //  Write FEM data.
        //
        typeMethods.r8mat_write(fem_node_filename, m, node_num, node_x);

        typeMethods.i4mat_write(fem_element_filename, element_order, element_num,
            element_node);

        Console.WriteLine("");
        Console.WriteLine("GMSH_TO_FEM:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}