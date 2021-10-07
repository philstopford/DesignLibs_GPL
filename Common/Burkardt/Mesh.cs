using System;

namespace Burkardt.MeshNS
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
            int node;
            int node_max;
            int node_min;
            int order;

            node_min = node_num + 1;
            node_max = -1;
            for (element = 0; element < element_num; element++)
            {
                for (order = 0; order < element_order; order++)
                {
                    node = element_node[order + element * element_order];
                    node_min = Math.Min ( node_min, node );
                    node_max = Math.Max ( node_max, node );                }
            }

            if (node_min == 1 && node_max == node_num)
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
            else if (node_min == 0 && node_max == node_num - 1)
            {
                Console.WriteLine("");
                Console.WriteLine("MESH_BASE_ZERO:");
                Console.WriteLine("  The element indexing appears to be 0-based!");
                Console.WriteLine("  No conversion is necessary.");
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

        public static void bandwidth_mesh(int element_order, int element_num, int[] element_node,
                ref int ml, ref int mu, ref int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
            //
            //  Discussion:
            //
            //    The quantity computed here is the "geometric" bandwidth determined
            //    by the finite element mesh alone.
            //
            //    If a single finite element variable is associated with each node
            //    of the mesh, and if the nodes and variables are numbered in the
            //    same way, then the geometric bandwidth is the same as the bandwidth
            //    of a typical finite element matrix.
            //
            //    The bandwidth M is defined in terms of the lower and upper bandwidths:
            //
            //      M = ML + 1 + MU
            //
            //    where 
            //
            //      ML = maximum distance from any diagonal entry to a nonzero
            //      entry in the same row, but earlier column,
            //
            //      MU = maximum distance from any diagonal entry to a nonzero
            //      entry in the same row, but later column.
            //
            //    Because the finite element node adjacency relationship is symmetric,
            //    we are guaranteed that ML = MU.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
            //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
            //
            //    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
            //
            //    Output, int *M, the bandwidth of the matrix.
            //
        {
            int element;
            int global_i;
            int global_j;
            int local_i;
            int local_j;

            ml = 0;
            mu = 0;

            for (element = 0; element < element_num; element++)
            {
                for (local_i = 0; local_i < element_order; local_i++)
                {
                    global_i = element_node[local_i + element * element_order];

                    for (local_j = 0; local_j < element_order; local_j++)
                    {
                        global_j = element_node[local_j + element * element_order];

                        mu = Math.Max(mu, global_j - global_i);
                        ml = Math.Max(ml, global_i - global_j);
                    }
                }
            }

            m = ml + 1 + mu;
        }
    }
}