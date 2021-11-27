using System;
using Burkardt.Table;
using Burkardt.FEM;
using Burkardt.Types;

namespace FEMGMeshTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM_TO_GMSH.
        //
        //  Discussion:
        //
        //    FEM_TO_GMSH converts a 1D, 2D or 3D mesh from FEM to GMSH format.
        //
        //  Usage:
        //
        //    fem_to_gmsh prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'_nodes.txt contains the node coordinates,
        //    * 'prefix'_elements.txt contains the element node connectivity.
        //    * 'prefix'.msh will contain the Gmsh version of the data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string prefix;

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("FEM_TO_GMSH");
        Console.WriteLine("  Convert a 1D, 2D or 3D mesh from FEM to GMSH format.");
        Console.WriteLine("");
        Console.WriteLine("  Read \"prefix\"_nodes.txt, node coordinates.");
        Console.WriteLine("  Read \"prefix\"_elements.txt, element node connectivity.");
        Console.WriteLine("");
        Console.WriteLine("  Create \"prefix\".msh, a corresponding Gmsh mesh file.");

        switch (args.Length)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("FEM_TO_GMSH:");
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
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string gmsh_filename = prefix + ".msh";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int node_num = h.n;
        int m = h.m;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + m + "");
        Console.WriteLine("  Number of nodes  = " + node_num + "");

        double[] node_x = typeMethods.r8mat_data_read(node_filename, m, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(m, node_num, node_x, 1, 1, m, 5,
            "  Portion of node coordinate data:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        int element_order = h.m;
        int element_num = h.n;

        switch (m)
        {
            case 1 when element_order == 2:
                break;
            case 1:
                break;
            case 2 when element_order == 3:
                break;
            case 2 when element_order == 6:
                break;
            case 2:
                Console.Write("");
                Console.Write("FEM_TO_GMSH - Fatal error!");
                Console.Write("2D mesh data must use 3 or 6 nodes.");
                return;
            case 3 when element_order == 4:
                break;
            case 3 when element_order == 10:
                break;
            case 3 when element_order == 20:
                break;
            case 3:
                Console.Write("");
                Console.Write("FEM_TO_GMSH - Fatal error!");
                Console.Write("  3D mesh data must use 4, 10, or 20 nodes.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order = " + element_order + "");
        Console.WriteLine("  Number of elements  = " + element_num + "");

        int[] element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 10, "  Initial portion of element data:");
        switch (m)
        {
            //
            //  Write out the Gmsh version of the data.
            //
            case 1:
                GMesh.gmsh_mesh1d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
                break;
            case 2:
                GMesh.gmsh_mesh2d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
                break;
            case 3:
                GMesh.gmsh_mesh3d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Created the GMSH file \"" + gmsh_filename + "\".");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("FEM_TO_GMSH:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}