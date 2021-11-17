using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangulationT3toT4Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_T3_TO_T4.
        //
        //  Discussion:
        //
        //    TRIANGULATION_T3_TO_T4 converts a T3 mesh to a T4 mesh.
        //
        //  Usage:
        //
        //    triangulation_t3_to_t4 prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * prefix_nodes.txt contains the node coordinates,
        //    * prefix_elements.txt contains the element definitions.
        //    * prefix_t4_nodes.txt will contain the T4 node coordinates,
        //    * prefix_t4_elements.txt will contain the T4 element definitions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 October 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim;
        int dim_num;
        int element;
        string element_t3_filename;
        string element_t4_filename;
        int[] element_node_t3;
        int[] element_node_t4;
        int element_num;
        int element_order_t3;
        int element_order_t4;
        int i;
        int j;
        int node;
        string node_t3_filename;
        string node_t4_filename;
        int node_num_t3;
        int node_num_t4;
        double[] node_xy_t3;
        double[] node_xy_t4;
        string prefix;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_T3_TO_T4");
        Console.WriteLine("  Read a 3-node T3 triangulation and");
        Console.WriteLine("  write out a 4-node T4 triangulation.");
        Console.WriteLine("");
        Console.WriteLine("  Read a dataset of NODE_NUM_T3 points in 2 dimensions.");
        Console.WriteLine("  Read an associated triangulation dataset of ELEMENT_NUM ");
        Console.WriteLine("  triangles which uses 3 nodes per triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Create new nodes which are triangle centroids,");
        Console.WriteLine("  generate new node and triangulation data for");
        Console.WriteLine("  4-node elements, and write them out.");
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
            Console.WriteLine("TRIANGULATION_T3_TO_T4:");
            Console.WriteLine("  Please enter the file prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_t3_filename = prefix + "_nodes.txt";
        element_t3_filename = prefix + "_elements.txt";
        node_t4_filename = prefix + "_t4_nodes.txt";
        element_t4_filename = prefix + "_t4_elements.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_t3_filename);
        dim_num = h.m;
        node_num_t3 = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_t3_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of nodes NODE_NUM_T3  = " + node_num_t3 + "");

        node_xy_t3 = typeMethods.r8mat_data_read(node_t3_filename, dim_num, node_num_t3);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_t3_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num_t3, node_xy_t3, 1, 1,
            dim_num, 5, "  Portion of coordinate data from file:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_t3_filename);
        element_order_t3 = h.m;
        element_num = h.n;

        if (element_order_t3 != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_T3_TO_T4 - Fatal error!");
            Console.WriteLine("  Data is not for a 3-node triangulation.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_t3_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Triangle order ELEMENT_ORDER_T3 = " + element_order_t3 + "");
        Console.WriteLine("  Number of elements ELEMENT_NUM  = " + element_num + "");

        element_node_t3 = typeMethods.i4mat_data_read(element_t3_filename, element_order_t3,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_t3_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order_t3, element_num, element_node_t3,
            1, 1, element_order_t3, 10, "  Portion of data read from file:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num_t3, element_order_t3, element_num, ref element_node_t3);
        //
        //  Allocate space.
        //
        node_num_t4 = node_num_t3 + element_num;
        element_order_t4 = 4;

        node_xy_t4 = new double[dim_num * node_num_t4];
        element_node_t4 = new int[element_order_t4 * element_num];
        //
        //  Create the centroid nodes.
        //
        for (element = 0; element < element_num; element++)
        {
            for (i = 0; i < 3; i++)
            {
                element_node_t4[i + element * element_order_t4] = element_node_t3[i + element * element_order_t3];
            }

            element_node_t4[3 + element * element_order_t4] = -1;
        }

        for (node = 0; node < node_num_t3; node++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                node_xy_t4[dim + node * dim_num] = node_xy_t3[dim + node * dim_num];
            }
        }

        node = node_num_t3;
        for (element = 0; element < element_num; element++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                node_xy_t4[dim + node * dim_num] = 0.0;
                for (i = 0; i < 3; i++)
                {
                    j = element_node_t3[i + element * element_order_t3] - 1;
                    node_xy_t4[(dim + node * dim_num + node_xy_t4.Length ) % node_xy_t4.Length] += node_xy_t3[(dim + j * dim_num + node_xy_t3.Length ) % node_xy_t3.Length];
                }

                node_xy_t4[(dim + node * dim_num + node_xy_t4.Length ) % node_xy_t4.Length] /= 3.0;
                element_node_t4[(3 + element * element_order_t4 + element_node_t4.Length ) % element_node_t4.Length] = node;
            }

            node += 1;
        }

        typeMethods.i4mat_transpose_print(element_order_t4, element_num, element_node_t4,
            "  ELEMENT_NODE_T4");

        typeMethods.r8mat_transpose_print(dim_num, node_num_t4, node_xy_t4, "  NODE_XY_T4:");
        //
        //  Write out the node and triangle data for the quadratic mesh
        //
        typeMethods.r8mat_write(node_t4_filename, dim_num, node_num_t4, node_xy_t4);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the T4 node data to \"" + node_t4_filename + "\".");

        typeMethods.i4mat_write(element_t4_filename, element_order_t4, element_num,
            element_node_t4);

        Console.WriteLine("  Wrote the T4 element data to \"" + element_t4_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_T3_TO_T4:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}