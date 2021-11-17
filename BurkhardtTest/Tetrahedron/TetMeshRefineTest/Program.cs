using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshRefineTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_REFINE.
        //
        //  Discussion:
        //
        //    TET_MESH_REFINE refines a tetrahedral mesh of order 4 (linear).
        //
        //  Usage:
        //
        //    tet_mesh_refine prefix
        //
        //    where prefix is the common file prefix:
        //
        //    * prefix_nodes.txt,    the node coordinates;
        //    * prefix_elements.txt,    the element definitions.
        //    * prefix_ref_nodes.txt,    the new node coordinates;
        //    * prefix_ref_elements.txt,    the new element definitions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 0;
        int[] edge_data;
        int[] element_node1 = new int[1];
        int[] element_node2 = new int[1];
        int element_num1 = 0;
        int element_num2 = 0;
        int element_order = 0;
        string input_node_filename;
        string input_element_filename;
        int node_num1 = 0;
        int node_num2 = 0;
        double[] node_xyz1 = new double[1];
        double[] node_xyz2 = new double[1];
        string output_node_filename;
        string output_element_filename;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_REFINE");
        Console.WriteLine("  READ a tet mesh, REFINE it, and WRITE the new data.");
        Console.WriteLine("");
        Console.WriteLine("  READ:");
        Console.WriteLine("    a node dataset of NODE_NUM1 points in 3 dimensions.");
        Console.WriteLine("    a tet mesh of TETRA_NUM1 tets of order TET_ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  REFINE:");
        Console.WriteLine("    compute a new set of nodes and tets, which is an");
        Console.WriteLine("    eightfold refinement of the input mesh.");
        Console.WriteLine("");
        Console.WriteLine("  WRITE:");
        Console.WriteLine("    a node dataset of NODE_NUM2 points in 3 dimensions.");
        Console.WriteLine("    a tet mesh of 8*TETRA_NUM1 tets of order TET_ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  At the moment, this program only works for a linear");
        Console.WriteLine("  mesh (TET_ORDER=4).");
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
            Console.WriteLine("TET_MESH_REFINE:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        input_node_filename = prefix + "_nodes.txt";
        input_element_filename = prefix + "_elements.txt";
        output_node_filename = prefix + "_ref_nodes.txt";
        output_element_filename = prefix + "_ref_elements.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(input_node_filename);
        dim_num = h.m;
        node_num1 = h.n;

        if (dim_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_REFINE - Fatal error!");
            Console.WriteLine("  The spatial dimension must be 3.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + input_node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of nodes   = " + node_num1 + "");

        node_xyz1 = typeMethods.r8mat_data_read(input_node_filename, dim_num,
            node_num1);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + input_node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num1,
            node_xyz1, 1, 1, dim_num, 5, "  First 5 input nodes:");
        //
        //  Read the tet mesh data.
        //
        h = typeMethods.i4mat_header_read(input_element_filename);
        element_order = h.m;
        element_num1 = h.n;

        switch (element_order)
        {
            case 4:
                break;
            case 10:
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_REFINE - Fatal error!");
                Console.WriteLine("  The program cannot yet handel the 10-node case.");
                Console.WriteLine("  Try using the sequence:");
                Console.WriteLine("    TET_MESH_Q2L --> TET_MESH_REFINE --> TET_MESH_L2Q.");
                return;
            default:
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_REFINE - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 4 or order 10.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + input_element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron order = " + element_order + "");
        Console.WriteLine("  Number of tetras  = " + element_num1 + "");

        element_node1 = typeMethods.i4mat_data_read(input_element_filename, element_order,
            element_num1);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + input_element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num1,
            element_node1, 1, 1, element_order, 5, "  First 5 input tetrahedrons:");
        //
        //  Check for 1-based node-indexing, and convert it to 0-based.
        //
        TetMesh.tet_mesh_base_zero(node_num1, element_order, element_num1, ref element_node1);
        switch (element_order)
        {
            //
            //  Compute the refined mesh.
            //
            case 4:
                edge_data = new int[5 * 6 * element_num1];

                TetMesh_Refine.tet_mesh_order4_refine_size(node_num1, element_num1, element_node1,
                    ref node_num2, ref element_num2, ref edge_data);

                Console.WriteLine("  Number of refined nodes =  " + node_num2 + "");
                Console.WriteLine("  Number of refined tetras = " + element_num2 + "");

                node_xyz2 = new double[dim_num * node_num2];
                element_node2 = new int[element_order * element_num2];

                TetMesh_Refine.tet_mesh_order4_refine_compute(node_num1, element_num1, node_xyz1,
                    element_node1, node_num2, element_num2, edge_data, ref node_xyz2, ref element_node2);
                break;
            case 10:
                break;
        }

        //
        //  Print a small amount of the refined data.
        //
        typeMethods.r8mat_transpose_print_some(dim_num, node_num2, node_xyz2,
            1, 1, dim_num, 5, "  First 5 output nodes:");

        typeMethods.i4mat_transpose_print_some(element_order, element_num2, element_node2,
            1, 1, element_order, 5, "  First 5 output tetras");
        //
        //  Write out the node and tetra data for the refined mesh
        //
        typeMethods.r8mat_write(output_node_filename, dim_num, node_num2, node_xyz2);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the file \"" + output_node_filename + "\".");

        typeMethods.i4mat_write(output_element_filename, element_order, element_num2, element_node2);

        Console.WriteLine("  Wrote the file \"" + output_element_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_REFINE:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}