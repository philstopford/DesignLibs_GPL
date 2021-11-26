using System;
using System.Globalization;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshQualityTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_QUALITY.
        //
        //  Discussion:
        //
        //    TET_MESH_QUALITY determines quality measures for a tet mesh.
        //
        //  Usage:
        //
        //    tet_mesh_quality prefix
        //
        //    where prefix is the common file prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates;
        //    * prefix_elements.txt contains the element definitions.
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
        int i;
        string prefix;
        double value_max = 0;
        double value_mean = 0;
        double value_min = 0;
        double value_var = 0;

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_QUALITY");
        Console.WriteLine("  Compute tet mesh quality measures.");
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
            Console.WriteLine("TET_MESH_QUALITY:");
            Console.WriteLine("  Please enter the file prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int dim_num = h.m;
        int node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM =  " + dim_num + "");
        Console.WriteLine("  Number of points NODE_NUM  = " + node_num + "");

        if (dim_num != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_QUALITY - Fatal error!");
            Console.WriteLine("  Dataset must have spatial dimension 3.");
            return;
        }

        double[] node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xyz, 1, 1, 5, 5,
            "  5 by 5 portion of data read from file:");
        //
        //  Read the tetra data.
        //
        h = typeMethods.i4mat_header_read(element_filename);

        int element_order = h.m;
        int element_num = h.n;

        if (element_order != 4)
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_QUALITY - Fatal error!");
            Console.WriteLine("  Data is not for a 4 node tet mesh.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order      = " + element_order + "");
        Console.WriteLine("  Number of elements = " + element_num + "");

        int[] element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num,
            element_node, 1, 1, element_order, 10, "  Portion of ELEMENT_NODE:");
        //
        //  If the element information is 1-based, make it 0-based.
        //
        TetMesh.tet_mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Compute and print the quality measures.
        //
        Console.WriteLine("");
        Console.WriteLine("                           Minimum    Mean         Maximum    Variance");
        Console.WriteLine("");

        TetMesh.tet_mesh_quality1(node_num, node_xyz, element_order, element_num,
            element_node, ref value_min, ref value_mean, ref value_max, ref value_var);

        Console.WriteLine("  Quality measure 1 = "
                          + "  " + value_min.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_mean.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_max.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_var.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        TetMesh.tet_mesh_quality2(node_num, node_xyz, element_order, element_num,
            element_node, ref value_min, ref value_mean, ref value_max, ref value_var);

        Console.WriteLine("  Quality measure 2 = "
                          + "  " + value_min.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_mean.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_max.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_var.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        TetMesh.tet_mesh_quality3(node_num, node_xyz, element_order, element_num,
            element_node, ref value_min, ref value_mean, ref value_max, ref value_var);

        Console.WriteLine("  Quality measure 3 = "
                          + "  " + value_min.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_mean.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_max.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_var.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        TetMesh.tet_mesh_quality4(node_num, node_xyz, element_order, element_num,
            element_node, ref value_min, ref value_mean, ref value_max, ref value_var);

        Console.WriteLine("  Quality measure 4 = "
                          + "  " + value_min.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_mean.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_max.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_var.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        TetMesh.tet_mesh_quality5(node_num, node_xyz, element_order, element_num,
            element_node, ref value_min, ref value_mean, ref value_max, ref value_var);

        Console.WriteLine("  Quality measure 5 = "
                          + "  " + value_min.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_mean.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_max.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + value_var.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        int[] node_order = TetMesh.tet_mesh_node_order(element_order, element_num, element_node,
            node_num);

        int histo_num = typeMethods.i4vec_max(node_num, node_order);

        int[] histo_gram = typeMethods.i4vec_histogram(node_num, node_order, histo_num);

        Console.WriteLine("");
        Console.WriteLine("  Here is a numerical histogram of the order of");
        Console.WriteLine("  each node in the mesh, that is, the number of");
        Console.WriteLine("  tetrahedrons that include that node as a vertex.");
        Console.WriteLine("");
        Console.WriteLine("  For a regular mesh, most nodes would have the");
        Console.WriteLine("  same order.");
        Console.WriteLine("");
        Console.WriteLine("   Order  Number of Nodes");
        Console.WriteLine("");

        for (i = 0; i <= histo_num; i++)
        {
            if (histo_gram[i] != 0)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + histo_gram[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
            }
        }

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_QUALITY:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}