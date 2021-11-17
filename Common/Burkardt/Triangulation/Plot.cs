using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Plot
{
    public static void triangulation_plot(string filename, int node_num, double[] node_xy,
            int element_order, int element_num, int[] element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_PLOT plots a triangulation in SVG format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string FILENAME, the name of the output file.
        //
        //    Input, int NODE_NUM, the number of points.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the nodes.
        //
        //    Input, int ELEMENT_ORDER, the order of the element.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
        //    lists, for each element, the indices of the points that form the vertices 
        //    of the element.
        //
    {
        int i;
        int i4;
        int i4_max;
        int i4_min;
        int ii;
        int j;
        int j4;
        int j4_max;
        int j4_min;
        int node;
        int[] order6 = {0, 3, 1, 4, 2, 5};
        List<string> output = new();
        double x_max;
        double x_min;
        double x_scale;
        double y_max;
        double y_min;
        double y_scale;
        //
        //  Determine SCALE, the maximum data range.
        //
        x_max = node_xy[0 + 0 * 2];
        x_min = node_xy[0 + 0 * 2];
        for (j = 0; j < node_num; j++)
        {
            x_max = Math.Max(x_max, node_xy[0 + j * 2]);
            x_min = Math.Min(x_min, node_xy[0 + j * 2]);
        }

        x_scale = x_max - x_min;
        x_max += 0.05 * x_scale;
        x_min -= 0.05 * x_scale;
        x_scale = x_max - x_min;

        y_max = node_xy[1 + 0 * 2];
        y_min = node_xy[1 + 0 * 2];
        for (j = 0; j < node_num; j++)
        {
            y_max = Math.Max(y_max, node_xy[1 + j * 2]);
            y_min = Math.Min(y_min, node_xy[1 + j * 2]);
        }

        y_scale = y_max - y_min;
        y_max += 0.05 * y_scale;
        y_min -= 0.05 * y_scale;
        y_scale = y_max - y_min;

        i4_min = 1;
        j4_min = 1;
        if (x_scale < y_scale)
        {
            i4_max = (int) (0.5 + 500.0 * x_scale / y_scale);
            j4_max = 500;
        }
        else
        {
            i4_max = 500;
            j4_max = (int) (0.5 + 500.0 * y_scale / x_scale);
        }

        //
        //  Write that junk.
        //
        output.Add("<?xml version = \"1.0\" standalone=\"no\"?>");
        output.Add("");
        output.Add("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"");
        output.Add("  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">");
        output.Add("");
        output.Add("<svg");
        output.Add("  width=\"" + i4_max + "\"");
        output.Add("  height=\"" + j4_max + "\"");
        output.Add("  viewbox=\"" + i4_min
                                  + "," + j4_min
                                  + "," + i4_max
                                  + "," + j4_max + "\"");
        output.Add("  xmlns=\"http://www.w3.org/2000/svg\"");
        output.Add("  version=\"1.1\">");
        output.Add("  <desc>");
        output.Add("    Triangulation created by triangulation_svg.cpp");
        output.Add("  </desc>");

        for (j = 0; j < element_num; j++)
        {
            output.Add("  <polygon");
            output.Add("    fill=\"red\"");
            output.Add("    stroke=\"black\"");
            output.Add("    stroke-width=\"2\"");
            output.Add("    points=\"");

            switch (element_order)
            {
                case 3:
                {
                    for (i = 0; i < 3; i++)
                    {
                        node = element_node[i + j * element_order];
                        i4 = typeMethods.r8_to_i4(x_min, x_max, node_xy[0 + node * 2], i4_min, i4_max);
                        j4 = typeMethods.r8_to_i4(y_max, y_min, node_xy[1 + node * 2], j4_min, j4_max);
                        output.Add("      " + i4 + "," + j4 + "");
                    }

                    break;
                }
                case 4:
                {
                    for (i = 0; i < 3; i++)
                    {
                        node = element_node[i + j * element_order];
                        i4 = typeMethods.r8_to_i4(x_min, x_max, node_xy[0 + node * 2], i4_min, i4_max);
                        j4 = typeMethods.r8_to_i4(y_max, y_min, node_xy[1 + node * 2], j4_min, j4_max);
                        output.Add("      " + i4 + "," + j4 + "");
                    }

                    break;
                }
                case 6:
                {
                    for (i = 0; i < 6; i++)
                    {
                        ii = order6[i];
                        node = element_node[ii + j * element_order];
                        i4 = typeMethods.r8_to_i4(x_min, x_max, node_xy[0 + node * 2], i4_min, i4_max);
                        j4 = typeMethods.r8_to_i4(y_max, y_min, node_xy[1 + node * 2], j4_min, j4_max);
                        output.Add("      " + i4 + "," + j4 + "");
                    }

                    break;
                }
            }

            output.Add("  \" />");
        }

        output.Add("</svg>");

        File.WriteAllLines(filename, output);

        Console.WriteLine("");
        Console.WriteLine("  Graphics data written to file \"" + filename + "\"");

    }
}