using System;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.FEMIOTest
{
    class Program
    {
        static void Main(string[] args)
                    {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FEM_IO_TEST.
            //
            //  Discussion:
            //
            //    FEM_IO_TEST tests the FEM_IO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_IO_TEST:");
                Console.WriteLine("  Test the FEM_IO library.");

                test01();
                test02();

                Console.WriteLine("");
                Console.WriteLine("FEM_IO_TEST:");
                Console.WriteLine("  Normal end of execution.");
                Console.WriteLine("");
            }

            static void test01()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests FEM_READ.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                int dim_num = 0;
                string element_file_name = "ell_elements.txt";
                int[] element_node = new int[1];
                int element_num = 0;
                int element_order = 0;
                double[] node_coord = new double[1];
                string node_coord_file_name = "ell_nodes.txt";
                double[] node_data = new double[1];
                string node_data_file_name = "ell_values.txt";
                int node_data_num = 0;
                int node_num = 0;

                Console.WriteLine("");
                Console.WriteLine("TEST01");
                Console.WriteLine("  FEM_READ reads finite element data from files.");

                Console.WriteLine("");
                Console.WriteLine("  The node coordinate file name is \""
                     + node_coord_file_name + "\".");
                Console.WriteLine("  The element file name is \""
                     + element_file_name + "\".");
                Console.WriteLine("  The node data file name is \""
                     + node_data_file_name + "\".");

                IO.fem_header_read(node_coord_file_name, element_file_name,
                    node_data_file_name, ref dim_num, ref node_num, ref element_num,
                    ref element_order, ref node_data_num);

                IO.fem_header_print(dim_num, node_num, element_order, element_num,
                    node_data_num);

                IO.fem_data_read(node_coord_file_name, element_file_name,
                    node_data_file_name, dim_num, node_num, element_num,
                    element_order, node_data_num, ref node_coord, ref element_node, ref node_data);

                typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_coord, 1, 1,
                    dim_num, 10, "  First 10 node coordindates:");

                typeMethods.i4mat_transpose_print_some(element_order, element_num,
                    element_node, 1, 1, element_order, 10, "  First 10 elements");

                typeMethods.r8mat_transpose_print_some(node_data_num, node_num, node_data,
                    1, 1, node_data_num, 10, "  First 10 node data sets:");

            }

            static void test02()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEM_IO_TEST02 tests FEM_WRITE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                int DIM_NUM = 2;
                int NODE_NUM = 5;
                int ELEMENT_NUM = 3;
                int ELEMENT_ORDER = 3;
                int NODE_DATA_NUM = 2;

                string element_file_name = "tiny_elements.txt";
                int[] element_node =  {
                    1, 2, 4,
                    5, 4, 2,
                    2, 3, 5
                }
                ;
                double[] node_coord =  {
                    0.0, 0.0,
                    1.0, 0.0,
                    2.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0
                }
                ;
                string node_coord_file_name = "tiny_nodes.txt";
                double[] node_data =  {
                    1.0, 0.0,
                    0.8, 0.2,
                    0.6, 0.4,
                    0.9, 0.1,
                    0.5, 0.5
                }
                ;
                string node_data_file_name = "tiny_values.txt";

                Console.WriteLine("");
                Console.WriteLine("FEM_TEST02");
                Console.WriteLine("  Demonstrate the use of FEM_WRITE to write finite");
                Console.WriteLine("  element data to files.");

                Console.WriteLine("");
                Console.WriteLine("  The node coordinate file name is \""
                     + node_coord_file_name + "\".");
                Console.WriteLine("  The element file name is \""
                     + element_file_name + "\".");
                Console.WriteLine("  The node data file name is \""
                     + node_data_file_name + "\".");

                IO.fem_header_print(DIM_NUM, NODE_NUM, ELEMENT_ORDER, ELEMENT_NUM,
                    NODE_DATA_NUM);

                typeMethods.r8mat_transpose_print(DIM_NUM, NODE_NUM, node_coord,
                    "  Node coordindates:");

                typeMethods.i4mat_transpose_print(ELEMENT_ORDER, ELEMENT_NUM,
                    element_node, "  Elements");

                typeMethods.r8mat_transpose_print(NODE_DATA_NUM, NODE_NUM, node_data,
                    "  Node data sets:");

                IO.fem_write(node_coord_file_name, element_file_name,
                    node_data_file_name, DIM_NUM, NODE_NUM, ELEMENT_NUM,
                    ELEMENT_ORDER, NODE_DATA_NUM, node_coord, element_node, node_data);

            }
        }
    }
}