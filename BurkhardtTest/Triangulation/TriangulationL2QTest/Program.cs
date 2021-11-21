using System;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationL2QTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_L2Q.
        //
        //  Discussion:
        //
        //    TRIANGULATION_L2Q makes a quadratic triangulation from a linear one.
        //
        //    Thanks to Zhu Wang for pointing out a problem caused by a change
        //    in the ordering of elements in the triangle neighbor array, 25 August 2010.
        //
        //  Usage:
        //
        //    triangulation_l2q prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the linear node coordinates,
        //    * prefix_elements.txt contains the linear element definitions.
        //    * prefix_l2q_nodes.txt will contain the quadratic node coordinates,
        //    * prefix_l2q_elements.txt will contain the quadratic element definitions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int boundary_num;
        int edge_num;
        int element;
        string element_filename;
        int[] element_node1;
        int[] element_node2;
        int element_num;
        int[] element_neighbor;
        int element_order1;
        int element_order2;
        int element2;
        int i;
        int ii;
        int iii;
        int interior_num;
        int ip1;
        int k1;
        int k2;
        int m;
        int node;
        string node_filename;
        int node_num1;
        int node_num2;
        double[] node_xy1;
        double[] node_xy2;
        string node_l2q_filename;
        string element_l2q_filename;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_L2Q");
        Console.WriteLine("  Read a \"linear\" triangulation and");
        Console.WriteLine("  write out a \"quadratic\" triangulation.");
        Console.WriteLine("");
        Console.WriteLine("  Read a dataset of NODE_NUM1 points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangulation dataset of ELEMENT_NUM ");
        Console.WriteLine("  elements which uses 3 nodes per triangular element.");
        Console.WriteLine("");
        Console.WriteLine("  Create new nodes which are triangle midpoints,");
        Console.WriteLine("  generate new node and triangulation data for");
        Console.WriteLine("  quadratic 6-node elements, and write them out.");
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
            Console.WriteLine("TRIANGULATION_L2Q:");
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        node_l2q_filename = prefix + "_l2q_nodes.txt";
        element_l2q_filename = prefix + "_l2q_elements.txt";
        //
        //  Read the data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        m = h.m;
        node_num1 = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Number of nodes NODE_NUM1  = " + node_num1 + "");

        node_xy1 = typeMethods.r8mat_data_read(node_filename, m, node_num1);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(m, node_num1, node_xy1, 1, 1, m, 10,
            "  Portion of node coordinate data:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order1 = h.m;
        element_num = h.n;

        if (element_order1 != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_L2Q - Fatal error!");
            Console.WriteLine("  Data is not for a 3-node triangulation.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order = " + element_order1 + "");
        Console.WriteLine("  Number of elements  = " + element_num + "");

        element_node1 = typeMethods.i4mat_data_read(element_filename, element_order1,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order1, element_num, element_node1,
            1, 1, element_order1, 10, "  First 10 elements:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num1, element_order1, element_num, ref element_node1);
        //
        //  Determine the number of midside nodes that will be added.
        //
        boundary_num = Boundary.triangulation_order3_boundary_edge_count(element_num,
            element_node1);

        interior_num = (3 * element_num - boundary_num) / 2;
        edge_num = interior_num + boundary_num;
        Console.WriteLine("");
        Console.WriteLine("  Number of midside nodes to add = " + edge_num + "");
        //
        //  Allocate space.
        //
        node_num2 = node_num1 + edge_num;
        element_order2 = 6;

        node_xy2 = new double[m * node_num2];
        element_node2 = new int[element_order2 * element_num];
        //
        //  Build the element neighbor array.
        //
        element_neighbor = NeighborElements.triangulation_neighbor_triangles(element_order1,
            element_num, element_node1);

        typeMethods.i4mat_transpose_print(3, element_num, element_neighbor,
            "  Element_neighbor:");
        //
        //  Create the midside nodes.
        //
        for (element = 0; element < element_num; element++)
        {
            for (i = 0; i < 3; i++)
            {
                element_node2[i + element * 6] = element_node1[i + element * 3];
            }

            for (i = 3; i < 6; i++)
            {
                element_node2[i + element * 6] = -1;
            }
        }

        for (i = 0; i < m; i++)
        {
            for (node = 0; node < node_num1; node++)
            {
                node_xy2[i + node * 2] = node_xy1[i + node * 2];
            }

            for (node = node_num1; node < node_num2; node++)
            {
                node_xy2[i + node * 2] = -1.0;
            }
        }

        node_num2 = node_num1;

        Console.WriteLine("");
        Console.WriteLine("  Generate midside nodes");
        Console.WriteLine("");

        for (element = 1; element <= element_num; element++)
        {
            for (i = 0; i < 3; i++)
            {
                //
                //  CORRECTION #1 because element neighbor definition was changed.
                //
                iii = typeMethods.i4_wrap(i + 2, 0, 2);
                element2 = element_neighbor[iii + (element - 1) * 3];

                switch (element2)
                {
                    case > 0 when element2 < element:
                        continue;
                }

                ip1 = typeMethods.i4_wrap(i + 1, 0, 2);
                //
                //  Temporary RETRO FIX!
                //
                k1 = element_node2[i + (element - 1) * 6] + 1;
                k2 = element_node2[ip1 + (element - 1) * 6] + 1;

                string cout = "  " + node_num2.ToString().PadLeft(8);
                for (ii = 0; ii < 2; ii++)
                {
                    node_xy2[(node_xy2.Length + ii + node_num2 * m) % node_xy2.Length] = 0.5
                        * (node_xy1[(ii + (k1 - 1) * m + node_xy1.Length) % node_xy1.Length] + node_xy1[(ii + (k2 - 1) * m + node_xy1.Length) % node_xy1.Length]);
                    cout += "  " + node_xy2[(ii + node_num2 * m + node_xy2.Length) % node_xy2.Length].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);

                element_node2[3 + i + (element - 1) * 6] = node_num2;

                switch (element2)
                {
                    case > 0:
                    {
                        for (ii = 0; ii < 3; ii++)
                        {
                            //
                            //  CORRECTION #2 because element neighbor definition changed.
                            //
                            iii = typeMethods.i4_wrap(ii + 2, 0, 2);
                            if (element_neighbor[iii + (element2 - 1) * 3] == element)
                            {
                                element_node2[ii + 3 + (element2 - 1) * 6] = node_num2;
                            }
                        }

                        break;
                    }
                }

                node_num2 += 1;
            }
        }

        typeMethods.i4mat_transpose_print(element_order2, element_num, element_node2,
            "  ELEMENT_NODE2:");
        //
        //  Write out the node and element data for the quadratic mesh.
        //
        typeMethods.r8mat_transpose_print(m, node_num2, node_xy2, "  NODE_XY2:");

        typeMethods.r8mat_write(node_l2q_filename, m, node_num2, node_xy2);

        typeMethods.i4mat_write(element_l2q_filename, element_order2, element_num,
            element_node2);

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_L2Q:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}