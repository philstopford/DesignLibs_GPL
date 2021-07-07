using System;
using System.IO;
using Burkardt.Types;

namespace Burkardt.GMesh
{
    public static class IO
    {
        public static void gmsh_data_read(string gmsh_filename, int node_dim, int node_num,
            double[] node_x, int element_order, int element_num, int[] element_node )

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
            int k;
            int length = 0;
            int level;
            double x;

            try
            {
                input = File.ReadAllLines(gmsh_filename);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("GMSH_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open input file \"" + gmsh_filename + "\"");
                throw;
            }
            
            level = 0;

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
                        x = typeMethods.s_to_r8(text, ref length, ref ierror);
                        text = text.Substring(length);
                        node_x[i + j * node_dim] = x;
                    }

                    j = j + 1;
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

                    j = j + 1;
                }
            }
        }

        public static void gmsh_size_read(string gmsh_filename, ref int node_num, ref int node_dim,
        ref int element_num, ref int element_order )

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
            int k;
            int length = 0;
            int level;
            const double r8_big = 1.0E+30;
            double x;
            double x_max;
            double x_min;
            double y;
            double y_max;
            double y_min;
            double z;
            double z_max;
            double z_min;

            node_num = 0;
            node_dim = 0;

            x_max = -r8_big;
            x_min = +r8_big;
            y_max = -r8_big;
            y_min = +r8_big;
            z_max = -r8_big;
            z_min = +r8_big;

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

            level = 0;
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

                    x = typeMethods.s_to_r8(text, ref length, ref ierror);
                    x_min = Math.Min(x_min, x);
                    x_max = Math.Max(x_max, x);
                    text = text.Substring(length);

                    y = typeMethods.s_to_r8(text, ref length, ref ierror);
                    y_min = Math.Min(y_min, y);
                    y_max = Math.Max(y_max, y);
                    text = text.Substring(length);

                    z = typeMethods.s_to_r8(text, ref length, ref ierror);
                    z_min = Math.Min(z_min, z);
                    z_max = Math.Max(z_max, z);
                    text = text.Substring(length);
                }
            }

            //
            //  Make a very simple guess as to the dimensionality of the data.
            //
            node_dim = 3;
            if (z_max == z_min)
            {
                node_dim = 2;
                if (y_max == y_min)
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
                    k = 0;
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
                    break;
                }
            }
        }
    }
}