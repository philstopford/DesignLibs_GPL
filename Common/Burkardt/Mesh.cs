using System;

namespace Burkardt
{
    public static class Mesh
    {
        public static void mesh_base_one(int node_num, int element_order, int element_num,
                ref int[] element_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MESH_BASE_ONE ensures that the element definition is 1-based.
            //
            //  Discussion:
            //
            //    The ELEMENT_NODE array contains nodes indices that form elements.
            //    The convention for node indexing might start at 0 or at 1.
            //
            //    If this function detects 0-based indexing, it converts to 1-based.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2014
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
            int element;
            const int i4_huge = 2147483647;
            int node;
            int node_max;
            int node_min;
            int order;

            node_min = +i4_huge;
            node_max = -i4_huge;
            for (element = 0; element < element_num; element++)
            {
                for (order = 0; order < element_order; order++)
                {
                    node = element_node[order + element * element_order];
                    if (node < node_min)
                    {
                        node_min = node;
                    }

                    if (node_max < node)
                    {
                        node_max = node;
                    }
                }
            }

            if (node_min == 0 && node_max == node_num - 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ONE:");
                Console.WriteLine("  The element indexing appears to be 0-based!");
                Console.WriteLine("  This will be converted to 1-based.");
                for (element = 0; element < element_num; element++)
                {
                    for (order = 0; order < element_order; order++)
                    {
                        element_node[order + element * element_order] =
                            element_node[order + element * element_order] + 1;
                    }
                }
            }
            else if (node_min == 1 && node_max == node_num)
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ONE:");
                Console.WriteLine("  The element indexing appears to be 1-based!");
                Console.WriteLine("  No conversion is necessary.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ONE - Warning!");
                Console.WriteLine("  The element indexing is not of a recognized type.");
                Console.WriteLine("  NODE_MIN = " + node_min + "");
                Console.WriteLine("  NODE_MAX = " + node_max + "");
                Console.WriteLine("  NODE_NUM = " + node_num + "");
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
            //  Discussion:
            //
            //    The ELEMENT_NODE array contains nodes indices that form elements.
            //    The convention for node indexing might start at 0 or at 1.
            //    Since a C++ program will naturally assume a 0-based indexing, it is
            //    necessary to check a given element definition and, if it is actually
            //    1-based, to convert it.
            //
            //    This function attempts to detect 1-based node indexing and correct it.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 October 2014
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
            int element;
            const int i4_huge = 2147483647;
            int node;
            int node_max;
            int node_min;
            int order;

            node_min = +i4_huge;
            node_max = -i4_huge;
            for (element = 0; element < element_num; element++)
            {
                for (order = 0; order < element_order; order++)
                {
                    node = element_node[order + element * element_order];
                    if (node < node_min)
                    {
                        node_min = node;
                    }

                    if (node_max < node)
                    {
                        node_max = node;
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
                for (element = 0; element < element_num; element++)
                {
                    for (order = 0; order < element_order; order++)
                    {
                        element_node[order + element * element_order] =
                            element_node[order + element * element_order] - 1;
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