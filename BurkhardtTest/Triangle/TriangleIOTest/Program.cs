using System;
using Burkardt.Types;

namespace TriangleIOTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGLE_IO_TEST.
            //
            //  Discussion:
            //
            //    TRIANGLE_IO_TEST tests the TRIANGLE_IO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("TRIANGLE_IO_TEST");
            Console.WriteLine("  Test the TRIANGLE_IO library.");

            test01();
            test02();
            test03();
            test04();

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_IO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        public static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 gets the example node data and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] node_att;
            int node_att_num = 0;
            double[] node_coord;
            int node_dim = 0;
            string node_file = "example.node";
            int[] node_marker;
            int node_marker_num = 0;
            int node_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Get example node data, write to a node file.");
            //
            //  Get node example size.
            //
            typeMethods.triangle_node_size_example(ref node_num, ref node_dim, ref node_att_num,
                ref node_marker_num);
            //
            //  Print the sizes.
            //
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + node_num + "");
            Console.WriteLine("  Spatial dimension =" + node_dim + "");
            Console.WriteLine("  Number of node attributes = " + node_att_num + "");
            Console.WriteLine("  Number of node markers = " + node_marker_num + "");
            //
            //  Allocate memory for node data.
            //
            node_coord = new double[node_dim * node_num];
            node_att = new double[node_att_num * node_num];
            node_marker = new int[node_marker_num * node_num];
            //
            //  Get the node data.
            //
            typeMethods.triangle_node_data_example(node_num, node_dim, node_att_num,
                node_marker_num, ref node_coord, ref node_att, ref node_marker);
            //
            //  Print some of the data.
            //
            typeMethods.r8mat_transpose_print_some(node_dim, node_num, node_coord,
                1, 1, node_dim, 10, "  Coordinates for first 10 nodes:");

            typeMethods.r8mat_transpose_print_some(node_att_num, node_num, node_att,
                1, 1, node_att_num, 10, "  Attributes for first 10 nodes:");

            typeMethods.i4mat_transpose_print_some(node_marker_num, node_num, node_marker,
                1, 1, node_marker_num, 10, "  Markers for first 10 nodes:");
            //
            //  Write the node information to node file.
            //
            typeMethods.triangle_node_write(node_file, node_num, node_dim, node_att_num,
                node_marker_num, ref node_coord, ref node_att, ref node_marker);

            Console.WriteLine("");
            Console.WriteLine("  Node data written to file \"" + node_file + "\"");
        }

        public static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 gets the example element data and writes it to a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] element_att;
            int element_att_num = 0;
            string element_file = "example.ele";
            int[] element_node;
            int element_num = 0;
            int element_order = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Get example element data, write to an element file.");
            //
            //  Get element example size.
            //
            typeMethods.triangle_element_size_example(ref element_num, ref element_order,
                ref element_att_num);
            //
            //  Print the sizes.
            //
            Console.WriteLine("");
            Console.WriteLine("  Number of elements = " + element_num + "");
            Console.WriteLine("  Order of elements = " + element_order + "");
            Console.WriteLine("  Number of element attributes = " + element_att_num + "");
            //
            //  Allocate memory.
            //
            element_node = new int[element_order * element_num];
            element_att = new double[element_att_num * element_num];
            //
            //  Get the data.
            //
            typeMethods.triangle_element_data_example(element_num, element_order, element_att_num,
                ref element_node, ref element_att);
            //
            //  Print some of the data.
            //
            typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
                1, 1, element_order, 10, "  Node connectivity of first 10 elements:");

            typeMethods.r8mat_transpose_print_some(element_att_num, element_num, element_att,
                1, 1, element_att_num, 10, "  Attributes for first 10 elements:");
            //
            //  Write the node information to node file.
            //
            typeMethods.triangle_element_write(element_file, element_num, element_order,
                element_att_num, element_node, element_att);

            Console.WriteLine("");
            Console.WriteLine("  Element data written to file \"" + element_file + "\"");
        }

        public static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 reads the example node data from a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] node_att;
            int node_att_num = 0;
            double[] node_coord;
            int node_dim = 0;
            string node_file = "example.node";
            int[] node_marker;
            int node_marker_num = 0;
            int node_num = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  Read node data from a node file.");
            //
            //  Get the data size.
            //
            typeMethods.triangle_node_size_read(node_file, ref node_num, ref node_dim, ref node_att_num,
                ref node_marker_num);
            //
            //  Print the sizes.
            //
            Console.WriteLine("");
            Console.WriteLine("  Node data read from file \"" + node_file + "\"");
            Console.WriteLine("");
            Console.WriteLine("  Number of nodes = " + node_num + "");
            Console.WriteLine("  Spatial dimension = " + node_dim + "");
            Console.WriteLine("  Number of node attributes = " + node_att_num + "");
            Console.WriteLine("  Number of node markers = " + node_marker_num + "");
            //
            //  Allocate memory.
            //
            node_coord = new double[node_dim * node_num];
            node_att = new double[node_att_num * node_num];
            node_marker = new int[node_marker_num * node_num];
            //
            //  Get the data.
            //
            typeMethods.triangle_node_data_read(node_file, node_num, node_dim, node_att_num,
                node_marker_num, ref node_coord, ref node_att, ref node_marker);
            //
            //  Print some of the data.
            //
            typeMethods.r8mat_transpose_print_some(node_dim, node_num, node_coord,
                1, 1, node_dim, 10, "  Coordinates for first 10 nodes:");

            typeMethods.r8mat_transpose_print_some(node_att_num, node_num, node_att,
                1, 1, node_att_num, 10, "  Attributes for first 10 nodes:");

            typeMethods.i4mat_transpose_print_some(node_marker_num, node_num, node_marker,
                1, 1, node_marker_num, 10, "  Markers for first 10 nodes:");
        }

        public static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 reads the example element data from a file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] element_att;
            int element_att_num = 0;
            string element_file = "example.ele";
            int[] element_node;
            int element_num = 0;
            int element_order = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST04:");
            Console.WriteLine("  Read element data from an element file.");
            //
            //  Get data size.
            //
            typeMethods.triangle_element_size_read(element_file, ref element_num, ref element_order,
                ref element_att_num);
            //
            //  Print the sizes.
            //
            Console.WriteLine("");
            Console.WriteLine("  Element data read from file \"" + element_file + "\"");
            Console.WriteLine("");
            Console.WriteLine("  Number of elements = " + element_num + "");
            Console.WriteLine("  Element order = " + element_order + "");
            Console.WriteLine("  Number of element attributes = " + element_att_num + "");
            //
            //  Allocate memory.
            //
            element_node = new int[element_order * element_num];
            element_att = new double[element_att_num * element_num];
            //
            //  Get the data.
            //
            typeMethods.triangle_element_data_read(element_file, element_num, element_order,
                element_att_num, ref element_node, ref element_att);
            //
            //  Print some of the data.
            //
            typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
                1, 1, element_order, 10, "  Connectivity for first 10 elements:");

            typeMethods.r8mat_transpose_print_some(element_att_num, element_num, element_att,
                1, 1, element_att_num, 10, "  Attributes for first 10 elements:");
        }
    }
}