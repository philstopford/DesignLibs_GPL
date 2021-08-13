using System.Collections.Generic;
using System.IO;
using Burkardt.MeshNS;

namespace Burkardt.GMesh
{
    public static class Mesh3D
    {
        public static void gmsh_mesh3d_write(string gmsh_filename, int m, int node_num,
                double[] node_x, int element_order, int element_num, int[] element_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.
            //
            //  Discussion:
            //
            //    The node ordering for the 20 node element is not standard.
            //
            //    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:
            //
            //    1:    a
            //    2:        b
            //    3:            c
            //    4:                d
            //    5: (2*a  +b        )/3
            //    6: (  a+2*b        )/3
            //    7: (    2*b+  c    )/3
            //    8: (      b+2*c    )/3
            //    9: (  a    +2*c    )/3
            //   10: (2*a    +  c    )/3
            //   11: (2*a        +  d)/3
            //   12: (  a        +2*d)/3
            //   13: (     b     +2*d)/3
            //   14: (   2*b     +  d)/3
            //   15: (       +  c+2*d)/3
            //   16: (       +2*c+  d)/3
            //   17: (  a+  b+  c    )/3
            //   18: (  a+  b    +  d)/3
            //   19: (      b+  c+  d)/3
            //   20: (  a+      c+  d)/3
            //
            //    Leo Rebholz used the following ordering:
            //
            //    1:    a
            //    2:        b
            //    3:            c
            //    4:                d
            //    5: (2*a  +b        )/3
            //    6: (2*a    +  c    )/3
            //    7: (  a+2*b        )/3
            //    8: (  a    +2*c    )/3
            //    9: (  a+  b+  c    )/3
            //   10: (    2*b+  c    )/3
            //   11: (      b+2*c    )/3
            //   12: (2*a        +  d)/3
            //   13: (   2*b     +  d)/3
            //   14: (       +2*c+  d)/3
            //   15: (  a+  b    +  d)/3
            //   16: (      b+  c+  d)/3
            //   17: (  a+      c+  d)/3
            //   18: (  a        +2*d)/3
            //   19: (     b     +2*d)/3
            //   20: (       +  c+2*d)/3
            //
            //    Since the only 20 node data we have is from Leo, we will assume that
            //    all 20 node input data is in Leo's format, and needs to be converted
            //    to the Gmsh convention.
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
            int i2;
            int[] leo_to_gmsh =
            {
                0, 1, 2, 3, 4,
                6, 9, 10, 7, 5,
                11, 17, 18, 12, 19,
                13, 8, 14, 15, 16
            };
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
                gmsh.Add(node + 1
                              + "  " + node_x[0 + node * m]
                              + "  " + node_x[1 + node * m]
                              + "  " + node_x[2 + node * m] + "");
            }

            gmsh.Add("$EndNodes");

            if (element_order == 4)
            {
                element_type = 4;
            }
            else if (element_order == 10)
            {
                element_type = 11;
            }
            else if (element_order == 20)
            {
                element_type = 29;
            }

            tag_num = 2;
            tag1 = 0;
            gmsh.Add("$Elements");
            gmsh.Add(element_num + "");
            for (element = 0; element < element_num; element++)
            {
                string tmp = element + 1
                                     + "  " + element_type
                                     + "  " + tag_num
                                     + "  " + tag1
                                     + "  " + element + 1;
                if (element_order == 20)
                {
                    for (i = 0; i < element_order; i++)
                    {
                        i2 = leo_to_gmsh[i];
                        tmp += "  " + element_node[i2 + element * element_order];
                    }

                    gmsh.Add(tmp);
                }
                else
                {
                    for (i = 0; i < element_order; i++)
                    {
                        tmp += "  " + element_node[i + element * element_order];
                    }

                    gmsh.Add(tmp);
                }
            }

            gmsh.Add("$EndElements");

            File.WriteAllLines(gmsh_filename, gmsh);
        }
    }
}