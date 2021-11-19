using System;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace FEM3DProjectTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM3D_PROJECT.
        //
        //  Discussion:
        //
        //    FEM3D_PROJECT reads files defining a sampling of a (scalar or vector)
        //    function of 3 arguments, and a list of nodes and tetrahedral elements 
        //    to use for a finite element representation of the data.
        //
        //    It computes a set of finite element coefficients to be associated with
        //    the given finite element mesh, and writes that information to a file
        //    so that an FEM representation is formed by the node, element and value 
        //    files.
        //
        //  Usage:
        //
        //    fem3d_project sample_prefix fem_prefix
        //
        //    where 'sample_prefix' is the common prefix for the SAMPLE files:
        //
        //    * sample_prefix_nodes.txt,     the node coordinates where samples were taken,
        //    * sample_prefix_elements.txt,  the 4 nodes that make up each element;
        //    * sample_prefix_values.txt,    the sample values.
        //
        //    and 'fem_prefix' is the common prefix for the FEM files:
        //
        //    * fem_prefix_nodes.txt,    the node coordinates.
        //    * fem_prefix_elements.txt, the 4 nodes that make up each element;
        //    * fem_prefix_values.txt,   the values defined at each node (output).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_min;
        string fem_element_filename;
        int[] fem_element_node;
        int fem_element_num;
        int fem_element_order;
        int fem_node_dim;
        string fem_node_filename;
        int fem_node_num;
        double[] fem_node_xyz;
        string fem_prefix;
        double[] fem_value;
        int fem_value_dim;
        string fem_value_filename;
        int fem_value_num;
        int i;
        int j;
        string sample_element_filename;
        int[] sample_element_neighbor;
        int[] sample_element_node;
        int sample_element_num;
        int sample_element_order;
        string sample_prefix;
        int sample_node_dim;
        string sample_node_filename;
        int sample_node_num;
        double[] sample_node_xyz;
        int sample_value_dim;
        int sample_value_num;
        double[] sample_value;
        string sample_value_filename;

        Console.WriteLine("");
        Console.WriteLine("FEM3D_PROJECT");
        Console.WriteLine("");
        Console.WriteLine("  Read files defining a sampling of a function of 3 arguments.");
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
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("Enter the sample file prefix:");
            sample_prefix = Console.ReadLine();
        }

        try
        {
            fem_prefix = args[1];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("Enter the FEM file prefix:");
            fem_prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        sample_node_filename = sample_prefix + "_nodes.txt";
        sample_element_filename = sample_prefix + "_elements.txt";
        sample_value_filename = sample_prefix + "_values.txt";

        fem_node_filename = fem_prefix + "_nodes.txt";
        fem_element_filename = fem_prefix + "_elements.txt";
        fem_value_filename = fem_prefix + "_values.txt";
        //
        //  Read the SAMPLE NODE, ELEMENT and VALUE data.
        //
        TableHeader h = typeMethods.r8mat_header_read(sample_node_filename);

        sample_node_dim = h.m;
        sample_node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Sample node spatial dimension is " + sample_node_dim + "");
        Console.WriteLine("  Sample node number is            " + sample_node_num + "");

        if (sample_node_dim != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM3D_PROJECT - Fatal error!");
            Console.WriteLine("  Spatial dimension of the sample nodes is not 3.");
            return;
        }

        sample_node_xyz = typeMethods.r8mat_data_read(sample_node_filename, sample_node_dim,
            sample_node_num);

        h = typeMethods.i4mat_header_read(sample_element_filename);
        sample_element_order = h.m;
        sample_element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Sample element order is  " + sample_element_order + "");
        Console.WriteLine("  Sample element number is " + sample_element_num + "");

        if (sample_element_order != 4)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM3D_PROJECT - Fatal error!");
            Console.WriteLine("  The sample element order must be 4.");
            return;
        }
            
        sample_element_node = typeMethods.i4mat_data_read(sample_element_filename,
            sample_element_order, sample_element_num);

        element_min = typeMethods.i4mat_min(sample_element_order, sample_element_num,
            sample_element_node);

        switch (element_min)
        {
            case 1:
            {
                Console.WriteLine("");
                Console.WriteLine("  Converting 1-based sample element array to 0 base.");
                for (j = 0; j < sample_element_num; j++)
                {
                    for (i = 0; i < sample_element_order; i++)
                    {
                        sample_element_node[i + j * sample_element_order] -= 1;
                    }
                }

                break;
            }
        }

        h = typeMethods.r8mat_header_read(sample_value_filename);
        sample_value_dim = h.m;
        sample_value_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  The sample value dimension is    " + sample_value_dim + "");
        Console.WriteLine("  The sample value number is        " + sample_value_num + "");

        if (sample_value_num != sample_node_num)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM3D_PROJECT - Fatal error!");
            Console.WriteLine("  Number of sample values and nodes differ.");
            return;
        }

        sample_value = typeMethods.r8mat_data_read(sample_value_filename, sample_value_dim,
            sample_value_num);
        //
        //  Create the sample element neighbor array.
        //
        sample_element_neighbor = TetMesh.tet_mesh_neighbor_tets(sample_element_order,
            sample_element_num, sample_element_node);

        Console.WriteLine("");
        Console.WriteLine("  The element neighbor array has been computed.");
        //
        //  Read the FEM NODE and ELEMENT data.
        //
        h = typeMethods.r8mat_header_read(fem_node_filename);
        fem_node_dim = h.m;
        fem_node_num = h.n;
        Console.WriteLine("");
        Console.WriteLine("  The FEM node dimension is        " + fem_node_dim + "");
        Console.WriteLine("  The FEM node number is           " + fem_node_num + "");

        if (fem_node_dim != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM3D_PROJECT - Fatal error!");
            Console.WriteLine("  Spatial dimension of the nodes is not 3.");
            return;
        }

        fem_node_xyz = typeMethods.r8mat_data_read(fem_node_filename, fem_node_dim, fem_node_num);

        h = typeMethods.i4mat_header_read(fem_element_filename);
        fem_element_order = h.m;
        fem_element_num = h.n;

        Console.WriteLine("  The FEM element order is         " + fem_element_order + "");
        Console.WriteLine("  The FEM element number is        " + fem_element_num + "");

        if (fem_element_order != 4)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM3D_PROJECT - Fatal error!");
            Console.WriteLine("  The FEM element order is not 4.");
            return;
        }

        fem_element_node = typeMethods.i4mat_data_read(fem_element_filename, fem_element_order,
            fem_element_num);

        element_min = typeMethods.i4mat_min(fem_element_order, fem_element_num,
            fem_element_node);

        switch (element_min)
        {
            case 1:
            {
                Console.WriteLine("");
                Console.WriteLine("  Converting 1-based FEM element array to 0 base.");
                for (j = 0; j < fem_element_num; j++)
                {
                    for (i = 0; i < fem_element_order; i++)
                    {
                        fem_element_node[i + j * fem_element_order] -= 1;
                    }
                }

                break;
            }
        }

        //
        //  Compute the FEM values.
        //
        fem_value_dim = sample_value_dim;
        fem_value_num = fem_node_num;

        fem_value = FEM_3D_Transfer.fem3d_transfer(sample_node_num, sample_element_order,
            sample_element_num, sample_value_dim, sample_value_num,
            sample_node_xyz, sample_element_node, sample_element_neighbor, sample_value,
            fem_node_num, fem_element_order,
            fem_element_num, fem_value_dim, fem_value_num,
            fem_node_xyz, fem_element_node);
        //
        //  Write the FEM values.
        //
        typeMethods.r8mat_write(fem_value_filename, fem_value_dim, fem_value_num,
            fem_value);

        Console.WriteLine("");
        Console.WriteLine("  FEM value data written to \"" + fem_value_filename + "\"");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("FEM3D_PROJECT");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}