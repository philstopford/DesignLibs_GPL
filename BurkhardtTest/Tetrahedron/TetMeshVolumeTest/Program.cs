using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshVolumeTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_VOLUMES.
        //
        //  Discussion:
        //
        //    TET_MESH_VOLUMES determines the element volumes of a tet mesh.
        //
        //  Usage:
        //
        //    tet_mesh_volumes prefix
        //
        //    where
        //
        //    * prefix_nodes.txt contains nodal coordinates;
        //    * prefix_elements.txt contains the element definitions;
        //    * prefix_volumes.txt will contain the element volumes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element;
        string prefix;
        double[] tetra = new double[3 * 4];

        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_VOLUMES");
        Console.WriteLine("  Compute volume of each tetrahedron in a tet mesh..");
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
            Console.WriteLine("TET_MESH_VOLUMES:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string volume_filename = prefix + "_volumes.txt";
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
            Console.WriteLine("TET_MESH_VOLUMES - Fatal error!");
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
            Console.WriteLine("TET_MESH_VOLUMES - Fatal error!");
            Console.WriteLine("  Data is not for a 4 node tet mesh.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron order = " + element_order + "");
        Console.WriteLine("  Number of tetras  = " + element_num + "");

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
        double[] volume = new double[element_num];

        for (element = 0; element < element_num; element++)
        {
            int j;
            for (j = 0; j < 4; j++)
            {
                int i;
                for (i = 0; i < 3; i++)
                {
                    int node = element_node[j + element * element_order];
                    tetra[i + j * 3] = node_xyz[i + node * 3];
                }
            }

            volume[element] = Tetrahedron.tetrahedron_volume(tetra);
        }

        double volume_max = typeMethods.r8vec_max(element_num, volume);
        double volume_min = typeMethods.r8vec_min(element_num, volume);
        double volume_ave = typeMethods.r8vec_mean(element_num, volume);
        double volume_tot = typeMethods.r8vec_sum(element_num, volume);
        double volume_var = typeMethods.r8vec_variance(element_num, volume);

        Console.WriteLine("");
        Console.WriteLine("  Minimum:   " + volume_min + "");
        Console.WriteLine("  Average:   " + volume_ave + "");
        Console.WriteLine("  Maximum:   " + volume_max + "");
        Console.WriteLine("  Total:     " + volume_tot + "");
        Console.WriteLine("  Variance:  " + volume_var + "");

        typeMethods.r8mat_write(volume_filename, 1, element_num, volume);

        Console.WriteLine("");
        Console.WriteLine("  Full list of volumes written to \"" + volume_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_VOLUMES:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}