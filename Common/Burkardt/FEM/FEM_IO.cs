using System;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class IO
    {
        public static void fem_data_read(string node_coord_file_name, string element_file_name,
            string node_data_file_name, int dim_num, int node_num, int element_num,
            int element_order, int node_data_num, ref double[] node_coord, ref int[] element_node,
            ref double[] node_data)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM_DATA_READ reads data from a set of FEM files.
        //
        //  Discussion:
        //
        //    This program reads the node, element and data files that define
        //    a finite element geometry and data based on that geometry:
        //    * a set of nodes, 
        //    * a set of elements based on those nodes, 
        //    * a set of data values associated with each node.
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
        //  Parameters:
        //
        //    Input, string NODE_COORD_FILE_NAME, the name of the node
        //    coordinate file.  If this argument is not supplied, it will be requested.
        //    If the interactive response is blank, or otherwise defective, then the
        //    program terminates.
        //
        //    Input, string ELEMENT_FILE_NAME, the name of the element
        //    file.  If this argument is not supplied, it will be requested.  If the
        //    interactive response is blank, then the program will assume that no 
        //    element information is to be supplied.  (But the node coordinates must 
        //    be available and may be plotted.  And if a node data file is supplied, 
        //    then the data can be plotted against the node coordinates without using 
        //    any finite element structure.)
        //
        //    Input, string NODE_DATA_FILE_NAME, the name of the node 
        //    data file.  If this argument is not supplied, it will be requested.  If 
        //    the interactive response is blank, then the program will assume that 
        //    no node data information is to be supplied.  (But the node coordinates 
        //    will be available and may be plotted.
        //    And if an element file is supplied, then the elements can also be
        //    displayed.)
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int NODE_DATA_NUM, the number of data items per node.
        //
        //    Output, double **NODE_COORD. a pointer to a double[DIM_NUM*NODE_NUM], 
        //    the coordinates of nodes.
        //
        //    Output, int **ELEMENT_NODE, a point to int[ELEMENT_ORDER*ELEMENT_NUM]; 
        //    the global index of local node I in element J.
        //
        //    Output, double **NODE_DATA, a pointer to a double[NODE_DATA_NUM*NODE_NUM], 
        //    the data values associated with each node.
        //
        {
            if (typeMethods.s_len_trim(node_coord_file_name) <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_DATA_READ:");
                Console.WriteLine("  No node coordinate file name was supplied!");
                Console.WriteLine("  NO DATA WILL BE READ!");
                return;
            }

            node_coord = typeMethods.r8mat_data_read(node_coord_file_name, dim_num, node_num);

            if (0 < typeMethods.s_len_trim(element_file_name))
            {
                element_node = typeMethods.i4mat_data_read(element_file_name, element_order,
                    element_num);
            }

            if (0 < typeMethods.s_len_trim(node_data_file_name))
            {
                node_data = typeMethods.r8mat_data_read(node_data_file_name, node_data_num,
                    node_num);
            }

            return;
        }

        public static void fem_header_print(int dim_num, int node_num, int element_num,
            int element_order, int node_data_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM_HEADER_PRINT prints the header to set of FEM files.
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
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int NODE_DATA_NUM, the number of data items per node.
        //
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension         = " + dim_num + "");
            Console.WriteLine("  Number of nodes           = " + node_num + "");
            Console.WriteLine("  Number of elements        = " + element_num + "");
            Console.WriteLine("  Element order             = " + element_order + "");
            Console.WriteLine("  Number of node data items = " + node_data_num + "");

            return;
        }

        public static void fem_header_read(string node_coord_file_name, string element_file_name,
            string node_data_file_name, ref int dim_num, ref int node_num, ref int element_num,
            ref int element_order, ref int node_data_num)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM_HEADER_READ reads the sizes of arrays in a set of FEM files.
        //
        //  Discussion:
        //
        //    This program reads the node, element and data files that define
        //    a finite element geometry and data based on that geometry:
        //    * a set of nodes, 
        //    * a set of elements based on those nodes, 
        //    * a set of data values associated with each node.
        //    and returns the sizes DIM_NUM, NODE_NUM, ELEMENT_NUM, ELEMENT_ORDER,
        //    and NODE_DATA_NUM required to allocate space for these arrays.
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
        //  Parameters:
        //
        //    Input, string NODE_COORD_FILE_NAME, the name of the node
        //    coordinate file.  If this argument is not supplied, it will be requested.
        //    If the interactive response is blank, or otherwise defective, then the
        //    program terminates.
        //
        //    Input, string ELEMENT_FILE_NAME, the name of the element
        //    file.  If this argument is not supplied, it will be requested.  If the
        //    interactive response is blank, then the program will assume that no 
        //    element information is to be supplied.  (But the node coordinates must 
        //    be available and may be plotted.  And if a node data file is supplied, 
        //    then the data can be plotted against the node coordinates without using 
        //    any finite element structure.)
        //
        //    Input, string NODE_DATA_FILE_NAME, the name of the node 
        //    data file.  If this argument is not supplied, it will be requested.  If 
        //    the interactive response is blank, then the program will assume that 
        //    no node data information is to be supplied.  (But the node coordinates 
        //    will be available and may be plotted.
        //    And if an element file is supplied, then the elements can also be
        //    displayed.)
        //
        //    Output, int *DIM_NUM, the spatial dimension, inferred from the
        //    "shape" of the data in the node file.
        //
        //    Output, int *NODE_NUM, the number of nodes, inferred from the 
        //    number of lines of data in the node coordinate file.
        //
        //    Output, int *ELEMENT_NUM, the number of elements, inferred from the
        //    number of lines of data in the element file.
        //
        //    Output, int *ELEMENT_ORDER, the order of the elements, inferred from
        //    the number of items in the first line of the element file.
        //
        //    Output, int *NODE_DATA_NUM, the number of data items per node,
        //    inferred from the number of items in the first line of the node data file.
        //
        {
            int node_num2 = 0;

            if (typeMethods.s_len_trim(node_coord_file_name) <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_HEADER_READ - Error!");
                Console.WriteLine("  No node coordinate file name was supplied!");
                Console.WriteLine("  NO DATA WILL BE READ!");
                return;
            }

            if (typeMethods.s_len_trim(element_file_name) <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_HEADER_READ:");
                Console.WriteLine("  No element file name was supplied.");
                Console.WriteLine("  Therefore, no element data will be returned.");
            }

            if (typeMethods.s_len_trim(node_data_file_name) <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_HEADER_READ:");
                Console.WriteLine("  No node data file name was supplied!");
                Console.WriteLine("  Therefore, no node data will be returned.");
            }

            //
            //  Read the node coordinate file.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_coord_file_name);
            element_order = h.m;
            element_num = h.n;

            if (0 < typeMethods.s_len_trim(element_file_name))
            {
                typeMethods.i4mat_header_read(element_file_name);
            }
            else
            {
                element_order = 0;
                element_num = 0;
            }

            if (0 < typeMethods.s_len_trim(node_data_file_name))
            {
                h = typeMethods.r8mat_header_read(node_data_file_name);
                node_data_num = h.m;
                node_num2 = h.n;

                if (node_num2 != node_num)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  The number of nodes in the node coordinate");
                    Console.WriteLine("  file is " + node_num + " but the number of nodes");
                    Console.WriteLine("  in the node data file is " + node_num2 + "");
                    Console.WriteLine("  Because of this, no node data will be stored.");
                    node_data_num = 0;
                }
            }
            else
            {
                node_data_num = 0;
            }
        }

        public static void fem_write(string node_coord_file_name, string element_file_name,
            string node_data_file_name, int dim_num, int node_num, int element_num,
            int element_order, int node_data_num, double[] node_coord,
            int[] element_node, double[] node_data )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM_WRITE writes data files associated with a finite element solution.
        //
        //  Discussion:
        //
        //    This program writes the node, element and data files that define
        //    a finite element geometry and data based on that geometry:
        //    * a set of nodes, 
        //    * a set of elements based on those nodes, 
        //    * a set of data values associated with each node.
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
        //  Parameters:
        //
        //    Input, string NODE_COORD_FILE_NAME, the name of the node
        //    coordinate file.  If this argument is empty, no node coordinate file will
        //    be written.
        //
        //    Input, string ELEMENT_FILE_NAME, the name of the element
        //    file.  If this argument is empty, no element file will be written.
        //
        //    Input, string NODE_DATA_FILE_NAME, the name of the node 
        //    data file.  If this argument is empty, no node data file will be written.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int NODE_DATA_NUM, the number of data items per node.
        //
        //    Input, double NODE_COORD[DIM_NUM*NODE_NUM], the coordinates 
        //    of nodes.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
        //    the global index of local node I in element J.
        //
        //    Input, double NODE_DATA[NODE_DATA_NUM*NODE_NUM], the data 
        //    values associated with each node.
        //
        {
            //
            //  Write the node coordinate file.
            //
            if (0 < typeMethods.s_len_trim(node_coord_file_name))
            {
                typeMethods.r8mat_write(node_coord_file_name, dim_num, node_num, node_coord);

                Console.WriteLine("");
                Console.WriteLine("FEM_WRITE wrote node coordinates to \""
                     + node_coord_file_name + "\".");
            }

            //
            //  Write the element file.
            //
            if (0 < typeMethods.s_len_trim(element_file_name))
            {
                typeMethods.i4mat_write(element_file_name, element_order, element_num, element_node);
            }

            //
            //  Write the node data file.
            //
            if (0 < typeMethods.s_len_trim(node_data_file_name))
            {
                typeMethods.r8mat_write(node_data_file_name, node_data_num, node_num, node_data);

                Console.WriteLine("");
                Console.WriteLine("FEM_WRITE wrote node data to \""
                     + node_data_file_name + "\".");
            }
        }

        public static void mesh_base_zero(int node_num, int element_order, int element_num,
            ref int[] element_node)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_BASE_ZERO ensures that the element definition is zero-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
        //    definitions.
        //
        {
            const int i4_huge = 2147483647;

            int node_min = +i4_huge;
            int node_max = -i4_huge;

            for (int j = 0; j < element_num; j++)
            {
                for (int i = 0; i < element_order; i++)
                {
                    int t = element_node[i + j * element_order];
                    if (t < node_min)
                    {
                        node_min = t;
                    }

                    if (node_max < t)
                    {
                        node_max = t;
                    }
                }
            }

            if (node_min == 0 && node_max == node_num - 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ZERO:");
                Console.WriteLine("  The element indexing appears to be 0-based!");
                Console.WriteLine("  No conversion is necessary.");

            }
            else if (node_min == 1 && node_max == node_num)
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ZERO:");
                Console.WriteLine("  The element indexing appears to be 1-based!");
                Console.WriteLine("  This will be converted to 0-based.");
                for (int j = 0; j < element_num; j++)
                {
                    for (int i = 0; i < element_order; i++)
                    {
                        element_node[i + j * element_order] = element_node[i + j * element_order] - 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ZERO - Warning!");
                Console.WriteLine("  The element indexing is not of a recognized type.");
                Console.WriteLine("  NODE_MIN = " + node_min + "");
                Console.WriteLine("  NODE_MAX = " + node_max + "");
                Console.WriteLine("  NODE_NUM = " + node_num + "");
            }
        }
    }
}