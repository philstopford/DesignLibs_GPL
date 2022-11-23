using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.GMesh;

public static class IO
{
    public static void gmsh_data_read(string gmsh_filename, int node_dim, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GMSH_DATA_READ reads data from a GMSH file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, character *GMSH_FILENAME, the GMSH filename.
        //
        //    Input, int NODE_DIM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_X[NODE_DIM*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    the nodes that make up each element.
        //
    {
        int i;
        bool ierror = false;
        string[] input;
        int j = 0;
        int length = 0;

        try
        {
            input = File.ReadAllLines(gmsh_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("GMSH_DATA_READ - Fatal error!");
            Console.WriteLine("  Could not open input file \"" + gmsh_filename + "\"");
            throw;
        }

        int level = 0;

        int index = 0;

        string text = "";
        while (true)
        {
            try
            {
                text = input[index];
                index++;
            }
            catch
            {
                break;
            }

            if (level == 0)
            {
                if (text.StartsWith("$Nodes"))
                {
                    level = 1;
                }
            }
            else if (level == 1)
            {
                typeMethods.s_to_i4(text, ref length, ref ierror);
                level = 2;
                j = 0;
            }
            else if (level == 2)
            {
                if (text.StartsWith("$EndNodes"))
                {
                    break;
                }

                typeMethods.s_to_i4(text, ref length, ref ierror);
                text = text.Substring(length);
                for (i = 0; i < node_dim; i++)
                {
                    double x = typeMethods.s_to_r8(text, ref length, ref ierror);
                    text = text.Substring(length);
                    node_x[i + j * node_dim] = x;
                }

                j += 1;
            }
        }

        //
        //  Now read element information.
        //
        level = 0;
        while (true)
        {
            try
            {
                text = input[index];
                index++;
            }
            catch
            {
                break;
            }

            if (level == 0)
            {
                if (text.StartsWith("$Elements"))
                {
                    level = 1;
                }
            }
            else if (level == 1)
            {
                typeMethods.s_to_i4(text, ref length, ref ierror);
                level = 2;
                j = 0;
            }
            else if (level == 2)
            {
                if (text.StartsWith("$EndElements"))
                {
                    break;
                }

                int k;
                for (k = 1; k <= 5; k++)
                {
                    typeMethods.s_to_i4(text, ref length, ref ierror);
                    text = text.Substring(length);
                }

                for (i = 0; i < element_order; i++)
                {
                    k = typeMethods.s_to_i4(text, ref length, ref ierror);
                    text = text.Substring(length);
                    element_node[i + j * element_order] = k;
                }

                j += 1;
            }
        }
    }

    public static void gmsh_size_read(string gmsh_filename, ref int node_num, ref int node_dim,
            ref int element_num, ref int element_order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GMSH_SIZE_READ reads sizes from a GMSH file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string GMSH_FILENAME, the GMSH filename.
        //
        //    Output, int &NODE_NUM, the number of nodes.
        //
        //    Output, int &NODE_DIM, the spatial dimension.
        //
        //    Output, int &ELEMENT_NUM, the number of elements.
        //
        //    Output, int &ELEMENT_ORDER, the order of the elements.
        //
    {
        bool ierror = false;
        string[] input;
        int length = 0;

        node_num = 0;
        node_dim = 0;

        double x_max = -typeMethods.r8_big();
        double x_min = +typeMethods.r8_big();
        double y_max = -typeMethods.r8_big();
        double y_min = +typeMethods.r8_big();
        double z_max = -typeMethods.r8_big();
        double z_min = +typeMethods.r8_big();

        try
        {
            input = File.ReadAllLines(gmsh_filename);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("GMSH_SIZE_READ - Fatal error!");
            Console.WriteLine("  Could not open input file \"" + gmsh_filename + "\"");
            throw;
        }

        int level = 0;
        int index = 0;
        string text = "";

        while (true)
        {
            try
            {
                text = input[index];
                index++;
            }
            catch
            {
                break;
            }

            if (level == 0)
            {
                if (text.StartsWith("$Nodes"))
                {
                    level = 1;
                }
            }
            else if (level == 1)
            {
                node_num = typeMethods.s_to_i4(text, ref length, ref ierror);
                level = 2;
            }
            else if (level == 2)
            {
                if (text.StartsWith("$EndNodes"))
                {
                    break;
                }

                typeMethods.s_to_i4(text, ref length, ref ierror);
                text = text.Substring(length);

                double x = typeMethods.s_to_r8(text, ref length, ref ierror);
                x_min = Math.Min(x_min, x);
                x_max = Math.Max(x_max, x);
                text = text.Substring(length);

                double y = typeMethods.s_to_r8(text, ref length, ref ierror);
                y_min = Math.Min(y_min, y);
                y_max = Math.Max(y_max, y);
                text = text.Substring(length);

                double z = typeMethods.s_to_r8(text, ref length, ref ierror);
                z_min = Math.Min(z_min, z);
                z_max = Math.Max(z_max, z);
                text = text.Substring(length);
            }
        }

        //
        //  Make a very simple guess as to the dimensionality of the data.
        //
        node_dim = 3;
        if (Math.Abs(z_max - z_min) <= typeMethods.r8_epsilon())
        {
            node_dim = 2;
            if (Math.Abs(y_max - y_min) <= typeMethods.r8_epsilon())
            {
                node_dim = 1;
            }
        }

        //
        //  Now read element information.
        //
        level = 0;

        while (true)
        {
            try
            {
                text = input[index];
                index++;
            }
            catch
            {
                break;
            }

            if (level == 0)
            {
                if (text.StartsWith("$Elements"))
                {
                    level = 1;
                }
            }
            else if (level == 1)
            {
                element_num = typeMethods.s_to_i4(text, ref length, ref ierror);
                level = 2;
            }
            else if (level == 2)
            {
                if (text.StartsWith("$EndElements"))
                {
                    break;
                }

                /*
                for (;;)
                {
                    typeMethods.s_to_i4(text, ref length, ref ierror);
                    text = text.Substring(length);
                    if (ierror)
                    {
                        break;
                    }

                    k = k + 1;
                }

                element_order = k - 5;
                */

                string[] tokens = text.Split(' ');
                element_order += tokens.Count(t => t != "");

                element_order -= 5;
                break;
            }
        }
    }

    public static void gmsh_write(string gmsh_filename, int dim_num, int node_num,
            double[] node_xyz, int element_order, int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GMSH_WRITE writes the tetrahedral mesh data as a Gmsh mesh file.
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
        //    03 June 2013
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
        //    Input, string GMSH_FILENAME, the name of the Gmsh file 
        //    to create.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, inte NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
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
        List<string> gmsh = new();
        int[] leo_to_gmsh =  {
                0, 1, 2, 3, 4,
                6, 9, 10, 7, 5,
                11, 17, 18, 12, 19,
                13, 8, 14, 15, 16
            }
            ;
        int node;

        gmsh.Add("$MeshFormat");
        gmsh.Add("2.2 0 8");
        gmsh.Add("$EndMeshFormat");

        gmsh.Add("$Nodes");
        gmsh.Add(node_num + "");
        for (node = 0; node < node_num; node++)
        {
            gmsh.Add(node + 1
                          + "  " + node_xyz[0 + node * 3]
                          + "  " + node_xyz[1 + node * 3]
                          + "  " + node_xyz[2 + node * 3] + "");
        }

        gmsh.Add("$EndNodes");

        element_type = element_order switch
        {
            4 => 4,
            10 => 11,
            20 => 29,
            _ => element_type
        };

        const int tag_num = 2;
        const int tag1 = 0;
        gmsh.Add("$Elements");
        gmsh.Add(element_num + "");
        for (element = 0; element < element_num; element++)
        {
            string tmp = element + 1
                                 + "  " + element_type
                                 + "  " + tag_num
                                 + "  " + tag1
                                 + "  " + (element + 1);
            int i;
            switch (element_order)
            {
                case 20:
                {
                    for (i = 0; i < element_order; i++)
                    {
                        int 
                            i2 = leo_to_gmsh[i];
                        tmp += "  " + element_node[i2 + element * element_order];
                    }

                    gmsh.Add(tmp);
                    break;
                }
                default:
                {
                    for (i = 0; i < element_order; i++)
                    {
                        tmp += "  " + element_node[i + element * element_order];
                    }

                    gmsh.Add(tmp);
                    break;
                }
            }
        }

        gmsh.Add("$EndElements");

        File.WriteAllLines(gmsh_filename, gmsh);
    }
}