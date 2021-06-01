using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.FEM
{
    public static class GMesh
    {
        public static void gmsh_mesh1d_write(string gmsh_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node)
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
            //
            //  Detect and correct 0-based node indexing.
            //
            mesh_base_one(node_num, element_order, element_num, ref element_node);

            List<string> lines = new List<string>();
            
            //
            //  Write the data.
            //
            lines.Add("$MeshFormat");
            lines.Add("2.2 0 8");
            lines.Add("$EndMeshFormat");

            lines.Add("$Nodes");
            lines.Add(node_num + "");
            for (int node = 0; node < node_num; node++)
            {
                lines.Add(node + 1
                               + "  " + node_x[0 + node * m]
                               + "  0.0  0.0");
            }

            lines.Add("$EndNodes");

            int element_type = 1;

            int tag_num = 2;
            int tag1 = 0;
            lines.Add("$Elements");
            lines.Add(element_num + "");
            for (int element = 0; element < element_num; element++)
            {
                string line = element + 1
                     + "  " + element_type
                     + "  " + tag_num
                     + "  " + tag1
                     + "  " + element + 1;
                for (int i = 0; i < element_order; i++)
                {
                    line += "  " + element_node[i + element * element_order];
                }

                lines.Add(line);
            }

            lines.Add("$EndElements");

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
            int element_type = 0;
            //
            //  Detect and correct 0-based node indexing.
            //
            mesh_base_one(node_num, element_order, element_num, ref element_node);

            List<string> lines = new List<string>();
            //
            //  Write the data.
            //
            lines.Add("$MeshFormat");
            lines.Add("2.2 0 8");
            lines.Add("$EndMeshFormat");

            lines.Add("$Nodes");
            lines.Add(node_num + "");
            for (int node = 0; node < node_num; node++)
            {
                lines.Add(node + 1
                               + "  " + node_x[0 + node * m]
                               + "  " + node_x[1 + node * m]
                               + "  0.0");
            }

            lines.Add("$EndNodes");

            if (element_order == 3)
            {
                element_type = 2;
            }
            else if (element_order == 6)
            {
                element_type = 9;
            }

            int tag_num = 2;
            int tag1 = 0;
            lines.Add("$Elements");
            lines.Add(element_num + "");
            for (int element = 0; element < element_num; element++)
            {
                string line = element + 1
                     + "  " + element_type
                     + "  " + tag_num
                     + "  " + tag1
                     + "  " + element + 1;
                for (int i = 0; i < element_order; i++)
                {
                    line += "  " + element_node[i + element * element_order];
                }

                lines.Add(line);
            }

            lines.Add("$EndElements");

            File.WriteAllLines(gmsh_filename, lines);
        }

        public static void gmsh_mesh3d_write(string gmsh_filename, int m, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )
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
            int element_type = 0;
            int[] leo_to_gmsh =  {
                0, 1, 2, 3, 4,
                6, 9, 10, 7, 5,
                11, 17, 18, 12, 19,
                13, 8, 14, 15, 16
            }
            ;
            //
            //  Detect and correct 0-based node indexing.
            //
            mesh_base_one(node_num, element_order, element_num, ref element_node);

            //
            //  Write the data.
            //
            List<string> lines = new List<string>();
            lines.Add("$MeshFormat");
            lines.Add("2.2 0 8");
            lines.Add("$EndMeshFormat");

            lines.Add("$Nodes");
            lines.Add(node_num + "");
            for (int node = 0; node < node_num; node++)
            {
                lines.Add(node + 1
                               + "  " + node_x[0 + node * m]
                               + "  " + node_x[1 + node * m]
                               + "  " + node_x[2 + node * m] + "");
            }

            lines.Add("$EndNodes");

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

            int tag_num = 2;
            int tag1 = 0;
            lines.Add("$Elements");
            lines.Add(element_num + "");
            for (int element = 0; element < element_num; element++)
            {
                string line =  element + 1
                     + "  " + element_type
                     + "  " + tag_num
                     + "  " + tag1
                     + "  " + element + 1;
                if (element_order == 20)
                {
                    for (int i = 0; i < element_order; i++)
                    {
                        int i2 = leo_to_gmsh[i];
                        line += "  " + element_node[i2 + element * element_order];
                    }

                }
                else
                {
                    for (int i = 0; i < element_order; i++)
                    {
                        line += "  " + element_node[i + element * element_order];
                    }

                }
                lines.Add(line + "");
            }

            lines.Add("$EndElements");

            File.WriteAllLines(gmsh_filename, lines);
        }

        static void mesh_base_one(int node_num, int element_order, int element_num,
            ref int[] element_node)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MESH_BASE_ONE ensures that the element definition is one-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 October 2014
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
            for (int element = 0; element < element_num; element++)
            {
                for (int order = 0; order < element_order; order++)
                {
                    int node = element_node[order + element * element_order];
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
                for (int element = 0; element < element_num; element++)
                {
                    for (int order = 0; order < element_order; order++)
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

            return;
        }
    }
}