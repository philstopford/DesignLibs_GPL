using System.Collections.Generic;
using System.IO;
using Burkardt.MeshNS;
using Burkardt.Types;

namespace Burkardt.GMesh
{
    public static class Mesh2D
    {
        public static int[] gmsh_mesh2d_element_data_example(int element_num, int element_order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GMSH_MESH2D_ELEMENT_DATA_EXAMPLE returns element information for the example.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Output, int GMSH_MESH2D_ELEMENT_DATA_EXAMPLE[ELEMENT_ORDER*ELEMENT_NUM], 
            //    the indices of the nodes that make up each element.
            //
        {
            int[] element_node;
            int[] element_node_save =
                {
                    1, 2, 6,
                    7, 6, 2,
                    2, 3, 7,
                    8, 7, 3,
                    3, 4, 8,
                    9, 8, 4,
                    4, 5, 9,
                    10, 9, 5,
                    6, 7, 11,
                    12, 11, 7,
                    7, 8, 12,
                    13, 12, 8,
                    8, 9, 13,
                    14, 13, 9,
                    9, 10, 14,
                    15, 14, 10,
                    11, 12, 16,
                    17, 16, 12,
                    12, 13, 17,
                    18, 17, 13,
                    16, 17, 19,
                    20, 19, 17,
                    17, 18, 20,
                    21, 20, 18
                }
                ;

            element_node = typeMethods.i4mat_copy_new(element_order, element_num,
                element_node_save);

            return element_node;
        }

        public static void gmsh_mesh2d_element_size_example(ref int element_num, ref int element_order)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GMSH_MESH2D_ELEMENT_SIZE_EXAMPLE returns element sizes for the example.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int &ELEMENT_NUM, the number of elements.
            //
            //    Output, int &ELEMENT_ORDER, the order of the elements.
            //
        {
            element_num = 24;
            element_order = 3;
        }

        public static double[] gmsh_mesh2d_node_data_example(int node_num, int node_dim)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GMSH_MESH2D_NODE_DATA_EXAMPLE returns node information for the example.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int NODE_DIM, the spatial dimension.
            //
            //    Output, double GMSH_MESH2D_NODE_DATA_EXAMPLE[NODE_DIM*NODE_NUM], 
            //    the nodal coordinates.
            //
        {
            double[] node_coord;
            double[] node_coord_save =
                {
                    0.0, 0.0,
                    1.0, 0.0,
                    2.0, 0.0,
                    3.0, 0.0,
                    4.0, 0.0,
                    0.0, 1.0,
                    1.0, 1.0,
                    2.0, 1.0,
                    3.0, 1.0,
                    4.0, 1.0,
                    0.0, 2.0,
                    1.0, 2.0,
                    2.0, 2.0,
                    3.0, 2.0,
                    4.0, 2.0,
                    0.0, 3.0,
                    1.0, 3.0,
                    2.0, 3.0,
                    0.0, 4.0,
                    1.0, 4.0,
                    2.0, 4.0
                }
                ;

            node_coord = typeMethods.r8mat_copy_new(2, 21, node_coord_save);

            return node_coord;
        }

        public static void gmsh_mesh2d_node_size_example(ref int node_num, ref int node_dim)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GMSH_MESH2D_NODE_SIZE_EXAMPLE returns the sizes of node information for the example.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 October 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int &NODE_NUM, the number of nodes.
            //
            //    Output, int &NODE_DIM, the spatial dimension.
            //
        {
            node_num = 21;
            node_dim = 2;
        }

        public static void gmsh_mesh2d_write(string gmsh_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GMSH_MESH2D_WRITE writes 2d mesh data as a Gmsh mesh file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Christophe Geuzaine, Jean-Francois Remacle,
        //    Gmsh: a three-dimensional finite element mesh generator with
        //    built-in pre- and post-processing facilities,
        //    International Journal for Numerical Methods in Engineering,
        //    Volume 79, Number 11, pages 1309-1331, 2009.
        //
        //  Parameters:
        //
        //    Input, string GMSH_FILENAME, the name of the Gmsh file.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_X[M*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
        {
            int element;
            int element_type = 0;
            List<string> gmsh = new List<string>();
            int i;
            int node;
            int tag_num;
            int tag1;
            //
            //  Detect and correct 0-based node indexing.
            //
            Mesh.mesh_base_one(node_num, element_order, element_num, ref element_node);

            //
            //  Write the data.
            //
            gmsh.Add("$MeshFormat");
            gmsh.Add("2.2 0 8");
            gmsh.Add("$EndMeshFormat");

            gmsh.Add("$Nodes");
            gmsh.Add(node_num + "");
            for (node = 0; node < node_num; node++)
            {
                gmsh.Add((node + 1)
                     + "  " + node_x[0 + node * m]
                     + "  " + node_x[1 + node * m]
                    + "  0.0");
            }

            gmsh.Add("$EndNodes");

            if (element_order == 3)
            {
                element_type = 2;
            }
            else if (element_order == 6)
            {
                element_type = 9;
            }

            tag_num = 2;
            tag1 = 0;
            gmsh.Add("$Elements");
            gmsh.Add(element_num + "");
            for (element = 0; element < element_num; element++)
            {
                string tmp = (element + 1)
                     + "  " + element_type
                     + "  " + tag_num
                     + "  " + tag1
                     + "  " + (element + 1);
                for (i = 0; i < element_order; i++)
                {
                    tmp += "  " + element_node[i + element * element_order];
                }

                gmsh.Add(tmp);
            }

            gmsh.Add("$EndElements");

            File.WriteAllLines(gmsh_filename, gmsh);
        }
    }
}