using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace TriangulationNodeToElementTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGULATION_NODE_TO_ELEMENT.
        //
        //  Discussion:
        //
        //    TRIANGULATION_NODE_TO_ELEMENT averages node values to an element value.
        //
        //    This program is given a triangulation, along with the values of one or
        //    more data items associated with each node.
        //
        //    It produces a file containing the average of the nodal values for 
        //    each element.
        //
        //  Usage:
        //
        //    triangulation_node_to_element ( 'prefix' )
        //
        //    where
        //
        //    * 'prefix'_nodes.txt contains the node coordinates;
        //    * 'prefix'_elements.txt contains the element definitions 
        //      (this file is optional, and if missing, the elements will be generated
        //      by the program);
        //    * 'prefix'_values.txt contains the nodal values.
        //    * 'prefix'_element_values will contain the averaged element data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num;
        string element_filename;
        int[] element_node;
        int element_num;
        int element_order;
        double[] element_value;
        string element_value_filename;
        int i;
        int j;
        int k;
        int ni;
        string node_filename;
        int node_num;
        double[] node_xy;
        string prefix;
        int value_dim;
        string value_filename;
        int value_num;
        double[] value;

        Console.WriteLine("TRIANGULATION_NODE_TO_ELEMENT");
        Console.WriteLine("  Average nodal data to create element data.");
        Console.WriteLine("");
        Console.WriteLine("  This program expects three files:");
        Console.WriteLine("  * prefix_nodes.txt,    node coordinates,");
        Console.WriteLine("  * prefix_elements.txt, indices of nodes forming elements,");
        Console.WriteLine("  * prefix_values.txt,   data values at nodes,");
        Console.WriteLine("  and creates:");
        Console.WriteLine("  * prefix_element_values.txt, averaged data at elements.");
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
            Console.WriteLine("  Please enter the filename prefix.");

            prefix = Console.ReadLine();
        }

        //
        //  Create the filenames.
        //
        node_filename = prefix + "_nodes.txt";
        element_filename = prefix + "_elements.txt";
        element_value_filename = prefix + "element_values.txt";
        value_filename = prefix + "_values.txt";
        //
        //  Read the node data.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        dim_num = h.m;
        node_num = h.n;
        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + node_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of nodes NODE_NUM  = " + node_num + "");

        node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + node_filename + "\".");

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, 5, 5,
            "  5 by 5 portion of data read from file:");
        //
        //  Read the element data.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        element_order = h.m;
        element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine(" Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order ELEMENT_ORDER = " + element_order + "");
        Console.WriteLine("  Number of elements ELEMENT_NUM  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename,
            element_order, element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, 5, 5, "  5 by 5 portion of data read from file:");
        //
        //  Read the node value data.
        //
        h = typeMethods.r8mat_header_read(value_filename);
        value_dim = h.m;
        value_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + value_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Number of values per node VALUE_DIM = " + value_dim + "");
        Console.WriteLine("  Number of values VALUE_NUM  = " + value_num + "");

        value = typeMethods.r8mat_data_read(value_filename, value_dim, value_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + value_filename + "\".");

        typeMethods.r8mat_transpose_print_some(value_dim, value_num, value, 1, 1, 5, 5,
            "  Portion of data:");
        //
        //  Detect and correct 1-based node indexing.
        //
        Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
        //
        //  Create the element values data.
        //
        element_value = new double[value_dim * element_num];

        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < value_dim; i++)
            {
                element_value[i + j * value_dim] = 0.0;
            }
        }

        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < element_order; i++)
            {
                ni = element_node[i + j * element_order];
                for (k = 0; k < value_dim; k++)
                {
                    element_value[k + j * value_dim] += value[k + ni * value_dim];
                }
            }
        }

        for (j = 0; j < element_num; j++)
        {
            for (k = 0; k < value_dim; k++)
            {
                element_value[k + j * value_dim] /= element_order;
            }
        }

        //
        //  Write out the file.
        //
        typeMethods.r8mat_write(element_value_filename, value_dim, element_num, element_value);

        Console.WriteLine("");
        Console.WriteLine("  Element values written to '" + element_value_filename + "'");

        Console.WriteLine("");
        Console.WriteLine("TRIANGULATION_NODE_TO_ELEMENT:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}