using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshNeighborsTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_TET_NEIGHBORS.
        //
        //  Discussion:
        //
        //    TET_MESH_TET_NEIGHBORS manages the tet mesh neighbor calculation.
        //
        //    A tet mesh of order 4 or order 10 may be used.
        //
        //  Usage:
        //
        //    tet_mesh_tet_neighbors prefix
        //
        //    where prefix is the common file prefix:
        //
        //    * prefix_nodes.txt,    the node coordinates (not needed by this program);
        //    * prefix_elements.txt,    the linear element definitions.
        //    * prefix_element_neighbors.txt, the element neighbors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_neighbor;
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        int node_num = 0;
        string element_filename;
        string neighbor_filename;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TET_NEIGHBORS");
        Console.WriteLine("  Read a tet mesh dataset of TETRA_NUM");
        Console.WriteLine("  tetrahedrons, using 4 or 10 nodes.");
        Console.WriteLine("");
        Console.WriteLine("  Compute the tet mesh neighbors, and write this");
        Console.WriteLine("  information to a file");
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
            Console.WriteLine("TET_MESH_TET_NEIGHBORS:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        element_filename = prefix + "_elements.txt";
        neighbor_filename = prefix + "_element_neighbors.txt";
        //
        //  Read the tet mesh data.
        //
        TableHeader h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;

        if (element_order != 4 && element_order != 10)
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_TET_NEIGHBORS - Fatal error!");
            Console.WriteLine("  The tet mesh must have order 4 or order 10.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron order = " + element_order + "");
        Console.WriteLine("  Number of tetras  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num,
            element_node, 1, 1, element_order, 5, "  First 5 tetrahedrons:");
        //
        //  Detect and correct 1-based node indexing.
        //
        TetMesh.tet_mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Compute the neighbor information.
        //
        element_neighbor = TetMesh_Neighbors.tet_mesh_neighbor_tets(element_order, element_num,
            element_node);

        typeMethods.i4mat_transpose_print_some(4, element_num,
            element_neighbor, 1, 1, 4, 5, "  First 5 neighbor sets:");
        //
        //  Write the neighbor information to a file.
        //
        typeMethods.i4mat_write(neighbor_filename, 4, element_num, element_neighbor);

        Console.WriteLine("");
        Console.WriteLine("  Created the file \"" + neighbor_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_TET_NEIGHBORS:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}