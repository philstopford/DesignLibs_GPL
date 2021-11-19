using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshGMeshTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_TO_GMSH.
        //
        //  Discussion:
        //
        //    TET_MESH_TO_GMSH converts a tet mesh, using 4 or 10 node 
        //    tetrahedral elements, to a Gmsh file.
        //
        //  Usage:
        //
        //    tet_mesh_to_gmsh prefix
        //
        //    where prefix is the common file prefix:
        //
        //    * prefix_nodes.txt,       the node coordinates;
        //    * prefix_elements.txt,    the linear element definitions.
        //    * prefix.msh,             the Gmsh msh file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Christophe Geuzaine, Jean-Francois Remacle,
        //    Gmsh: a three-dimensional finite element mesh generator with
        //    built-in pre- and post-processing facilities,
        //    International Journal for Numerical Methods in Engineering,
        //    Volume 79, Number 11, pages 1309-1331, 2009.
        //
    {
        int dim_num;
        string element_filename;
        int[] element_node;
        int element_num;
        int element_order;
        string gmsh_filename;
        string node_filename;
        int node_num;
        double[] node_xyz;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TO_GMSH;");
        Console.WriteLine("  Convert a linear or quadratic tet mesh to Gmsh format");
        Console.WriteLine("");
        Console.WriteLine("  Read \"prefix\"_nodes.txt, node coordinates.");
        Console.WriteLine("  Read \"prefix\"_elements.txt, 4 or 10 node element definitions.");
        Console.WriteLine("");
        Console.WriteLine("  Create \"prefix\".msh, a corresponding Gmsh mesh file.");
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
            Console.WriteLine("TET_MESH_TO_GMSH;:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        gmsh_filename = prefix + ".msh";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        dim_num = h.m;
        node_num = h.n;
        if (dim_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_TO_GMSH; - Fatal error!");
            Console.WriteLine("  The spatial dimension must be 3.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of nodes   = " + node_num + "");

        node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num,
            node_xyz, 1, 1, dim_num, 5, "  First 5 nodes:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;
        switch (element_order)
        {
            case 4:
                break;
            case 10:
                break;
            case 20:
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_TO_GMSH; - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 4, 10, or 20.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order = " + element_order + "");
        Console.WriteLine("  Number of elements  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num,
            element_node, 1, 1, element_order, 5, "  First 5 elements:");
        //
        //  Use 1-based indexing.
        //
        TetMesh.tet_mesh_base_one(node_num, element_order, element_num, ref element_node);
        //
        //  Write the Gmsh file.
        //
        Burkardt.GMesh.IO.gmsh_write(gmsh_filename, dim_num, node_num, node_xyz, element_order,
            element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Created GMSH file \"" + gmsh_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TO_GMSH;");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}