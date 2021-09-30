using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.MeshNS;

namespace Burkardt.GMesh
{
    public class Mesh1D
    {
        public static void gmsh_mesh1d_write(string gmsh_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GMSH_MESH1D_WRITE writes 1d mesh data as a Gmsh mesh file.
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
            int element_type;
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
                    + "  0.0  0.0");
            }

            gmsh.Add("$EndNodes");

            element_type = 1;

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