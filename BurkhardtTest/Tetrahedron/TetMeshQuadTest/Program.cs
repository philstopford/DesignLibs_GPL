﻿using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TetMeshQuadTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TET_MESH_QUAD.
        //
        //  Discussion:
        //
        //    TET_MESH_QUAD estimates the integral of a scalar or vector function
        //    over a region defined by a tetrahedral mesh, using only the values
        //    of the function at the nodes of the mesh.
        //
        //  Usage:
        //
        //    tet_mesh_quad prefix
        //
        //    where
        //
        //    * 'prefix'_nodes.txt contains the node coordinates;
        //    * 'prefix'_elements.txt contains the element definitions.
        //    * 'prefix'_values.txt contains the function component values at each node.
        //
        //    The program will compute and print the values of QUAD.
        //    QUAD will be a scalar or vector of numbers, containing the estimates
        //    for the integral of each component of the function over the tet mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local parameters:
        //
        //    Local, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    lists the nodes that make up each element.
        //
        //    Local, int ELEMENT_NUM, the number of elements.
        //
        //    Local, int ELEMENT_ORDER, the order of the elements, either 4 or 10.
        //
        //    Local, int NODE_DIM, the spatial dimension.
        //
        //    Local, int NODE_NUM, the number of nodes.
        //
        //    Local, double NODE_XYZ[NODE_DIM*NODE_NUM], the point set.
        //
        //    Local, double VALUE[VALUE_DIM*VALUE_NUM], the VALUE_DIM components of a
        //    vector function evaluated at each of the NODE_NUM nodes.
        //
        //    Local, int VALUE_DIM, the number of components of the function.
        //
        //    Local, int VALUE_NUM, the number of points at which function values
        //    were supplied.  This must be equal to NODE_NUM.
        //
    {
        int element;
        int i;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_QUAD:");
        Console.WriteLine("  Estimate an integral over a region defined by a tetrahedral mesh.");
        Console.WriteLine("");
        //
        //  The first commandline argument is the filename prefix.
        //
        try
        {
            prefix = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_QUAD:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string value_filename = prefix + "_values.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int node_dim = h.m;
        int node_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension NODE_DIM = " + node_dim + "");
        Console.WriteLine("  Number of nodes NODE_NUM  = " + node_num + "");

        double[] node_xyz = typeMethods.r8mat_data_read(node_filename, node_dim, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(node_dim, node_num, node_xyz, 1, 1, node_dim, 5,
            "  First 5 nodes:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        int element_order = h.m;
        int element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order      = " + element_order + "");
        Console.WriteLine("  Number of elements = " + element_num + "");

        int[] element_node = typeMethods.i4mat_data_read(element_filename,
            element_order, element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 5, "  First 5 elements:");
        //
        //  Read the vaue data.
        //
        h = typeMethods.r8mat_header_read(value_filename);
        int value_dim = h.m;
        int value_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + value_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  VALUE dimension VALUE_DIM = " + value_dim + "");
        Console.WriteLine("  Number of nodes VALUE_NUM = " + value_num + "");

        if (value_num != node_num)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_QUAD - Fatal error!");
            Console.WriteLine("  Number of sets of values must equal number of  nodes.");
            return;
        }

        double[] value = typeMethods.r8mat_data_read(value_filename, value_dim, value_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + value_filename + "\".");

        typeMethods.r8mat_transpose_print_some(value_dim, value_num, value, 1, 1, value_dim, 5,
            "  First 5 values:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Compute the tetrahedron volumes.
        //
        double[] volume = new double[element_num];

        for (element = 0; element < element_num; element++)
        {
            volume[element] = 1.0 / 6.0 * Math.Abs(
                (node_xyz[0 + element_node[1 + element * element_order] * 3] -
                 node_xyz[0 + element_node[0 + element * element_order] * 3])
                * ((node_xyz[1 + element_node[2 + element * element_order] * 3] -
                    node_xyz[1 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[2 + element_node[3 + element * element_order] * 3] -
                      node_xyz[2 + element_node[0 + element * element_order] * 3])
                   - (node_xyz[2 + element_node[2 + element * element_order] * 3] -
                      node_xyz[2 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[1 + element_node[3 + element * element_order] * 3] -
                      node_xyz[1 + element_node[0 + element * element_order] * 3]))
                + (node_xyz[1 + element_node[1 + element * element_order] * 3] -
                   node_xyz[1 + element_node[0 + element * element_order] * 3])
                * ((node_xyz[2 + element_node[2 + element * element_order] * 3] -
                    node_xyz[2 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[0 + element_node[3 + element * element_order] * 3] -
                      node_xyz[0 + element_node[0 + element * element_order] * 3])
                   - (node_xyz[0 + element_node[2 + element * element_order] * 3] -
                      node_xyz[0 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[2 + element_node[3 + element * element_order] * 3] -
                      node_xyz[2 + element_node[0 + element * element_order] * 3]))
                + (node_xyz[2 + element_node[1 + element * element_order] * 3] -
                   node_xyz[2 + element_node[0 + element * element_order] * 3])
                * ((node_xyz[0 + element_node[2 + element * element_order] * 3] -
                    node_xyz[0 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[1 + element_node[3 + element * element_order] * 3] -
                      node_xyz[1 + element_node[0 + element * element_order] * 3])
                   - (node_xyz[1 + element_node[2 + element * element_order] * 3] -
                      node_xyz[1 + element_node[0 + element * element_order] * 3])
                   * (node_xyz[0 + element_node[3 + element * element_order] * 3] -
                      node_xyz[0 + element_node[0 + element * element_order] * 3])));
        }

        //
        //  Compute the integral estimate.
        //
        double[] quad = new double[value_dim];

        for (i = 0; i < value_dim; i++)
        {
            quad[i] = 0.0;
            for (element = 0; element < element_num; element++)
            {
                quad[i] += (value[i + element_node[0 + element * element_order] * value_dim]
                            + value[i + element_node[1 + element * element_order] * value_dim]
                            + value[i + element_node[2 + element * element_order] * value_dim]
                            + value[i + element_node[3 + element * element_order] * value_dim]) * volume[element] / 4.0;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Tet mesh volume = " + typeMethods.r8vec_sum(element_num, volume) + "");

        typeMethods.r8vec_print(value_dim, quad, "  Integral estimates:");

        Console.WriteLine("");
        Console.WriteLine("TET_MESH_QUAD:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}