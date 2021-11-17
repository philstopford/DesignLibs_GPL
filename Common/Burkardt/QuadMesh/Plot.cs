using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.QuadMesh;

public static class Plot
{
    public static void plot_q4_mesh(int node_num, int element_num, double[] node_xy,
            int[] element_node, int node_show, int element_show, string output_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLOT_Q4_MESH plots a Q4 mesh.
        //
        //  Discussion:
        //
        //    The triangulation is most usually a Delaunay triangulation,
        //    but this is not necessary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes that form the elements.
        //
        //    Input, int NODE_SHOW:
        //    0, do not show nodes;
        //    1, show nodes;
        //    2, show nodes and label them.
        //
        //    Input, int ELEMENT_SHOW:
        //    0, do not show elements;
        //    1, show elements;
        //    2, show elements and label them.
        //
        //    Input, string OUTPUT_FILENAME, the name of the output file.
        //
    {
        double ave_x;
        double ave_y;
        int circle_size;
        int delta;
        int e;
        int element;
        int element_order = 4;
        int i;
        int node;
        List<string> output_unit = new();
        double x_max;
        double x_min;
        int x_ps;
        int x_ps_max = 576;
        int x_ps_max_clip = 594;
        int x_ps_min = 36;
        int x_ps_min_clip = 18;
        double x_scale;
        double y_max;
        double y_min;
        int y_ps;
        int y_ps_max = 666;
        int y_ps_max_clip = 684;
        int y_ps_min = 126;
        int y_ps_min_clip = 108;
        double y_scale;
        //
        //  We need to do some figuring here, so that we can determine
        //  the range of the data, and hence the height and width
        //  of the piece of paper.
        //
        x_max = -typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (x_max < node_xy[0 + node * 2])
            {
                x_max = node_xy[0 + node * 2];
            }
        }

        x_min = typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[0 + node * 2] < x_min)
            {
                x_min = node_xy[0 + node * 2];
            }
        }

        x_scale = x_max - x_min;

        x_max += 0.05 * x_scale;
        x_min -= 0.05 * x_scale;
        x_scale = x_max - x_min;

        y_max = -typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (y_max < node_xy[1 + node * 2])
            {
                y_max = node_xy[1 + node * 2];
            }
        }

        y_min = typeMethods.r8_huge();
        for (node = 0; node < node_num; node++)
        {
            if (node_xy[1 + node * 2] < y_min)
            {
                y_min = node_xy[1 + node * 2];
            }
        }

        y_scale = y_max - y_min;

        y_max += 0.05 * y_scale;
        y_min -= 0.05 * y_scale;
        y_scale = y_max - y_min;

        if (x_scale < y_scale)
        {
            delta = (int)typeMethods.r8_nint((x_ps_max - x_ps_min)
                * (y_scale - x_scale) / (2.0 * y_scale));

            x_ps_max -= delta;
            x_ps_min += delta;

            x_ps_max_clip -= delta;
            x_ps_min_clip += delta;

            x_scale = y_scale;
        }
        else if (y_scale < x_scale)
        {
            delta = (int)typeMethods.r8_nint((y_ps_max - y_ps_min)
                * (x_scale - y_scale) / (2.0 * x_scale));

            y_ps_max -= delta;
            y_ps_min += delta;

            y_ps_max_clip -= delta;
            y_ps_min_clip += delta;

            y_scale = x_scale;
        }

        output_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
        output_unit.Add("%%Creator: plot_q4_mesh.C");
        output_unit.Add("%%Title: " + output_filename + "");

        output_unit.Add("%%Pages: 1");
        output_unit.Add("%%BoundingBox:  "
                        + x_ps_min + "  "
                        + y_ps_min + "  "
                        + x_ps_max + "  "
                        + y_ps_max + "");
        output_unit.Add("%%Document-Fonts: Times-Roman");
        output_unit.Add("%%LanguageLevel: 1");
        output_unit.Add("%%EndComments");
        output_unit.Add("%%BeginProlog");
        output_unit.Add("/inch {72 mul} def");
        output_unit.Add("%%EndProlog");
        output_unit.Add("%%Page:      1     1");
        output_unit.Add("save");
        output_unit.Add("%");
        output_unit.Add("% Set the RGB line color to very light gray.");
        output_unit.Add("%");
        output_unit.Add(" 0.9000 0.9000 0.9000 setrgbcolor");
        output_unit.Add("%");
        output_unit.Add("% Draw a gray border around the page.");
        output_unit.Add("%");
        output_unit.Add("newpath");
        output_unit.Add(x_ps_min + "  "
                                 + y_ps_min + "  moveto");
        output_unit.Add(x_ps_max + "  "
                                 + y_ps_min + "  lineto");
        output_unit.Add(x_ps_max + "  "
                                 + y_ps_max + "  lineto");
        output_unit.Add(x_ps_min + "  "
                                 + y_ps_max + "  lineto");
        output_unit.Add(x_ps_min + "  "
                                 + y_ps_min + "  lineto");
        output_unit.Add("stroke");
        output_unit.Add("%");
        output_unit.Add("% Set RGB line color to black.");
        output_unit.Add("%");
        output_unit.Add(" 0.0000 0.0000 0.0000 setrgbcolor");
        output_unit.Add("%");
        output_unit.Add("%  Set the font and its size:");
        output_unit.Add("%");
        output_unit.Add("/Times-Roman findfont");
        output_unit.Add("0.50 inch scalefont");
        output_unit.Add("setfont");
        output_unit.Add("%");
        output_unit.Add("%  Print a title:");
        output_unit.Add("%");
        output_unit.Add("%  210  702 moveto");
        output_unit.Add("%(Pointset) show");
        output_unit.Add("%");
        output_unit.Add("% Define a clipping polygon");
        output_unit.Add("%");
        output_unit.Add("newpath");
        output_unit.Add(x_ps_min_clip + "  "
                                      + y_ps_min_clip + "  moveto");
        output_unit.Add(x_ps_max_clip + "  "
                                      + y_ps_min_clip + "  lineto");
        output_unit.Add(x_ps_max_clip + "  "
                                      + y_ps_max_clip + "  lineto");
        output_unit.Add(x_ps_min_clip + "  "
                                      + y_ps_max_clip + "  lineto");
        output_unit.Add(x_ps_min_clip + "  "
                                      + y_ps_min_clip + "  lineto");
        output_unit.Add("clip newpath");
        circle_size = node_num switch
        {
            //
            //  Draw the nodes.
            //
            <= 200 => 5,
            <= 500 => 4,
            <= 1000 => 3,
            <= 5000 => 2,
            _ => 1
        };

        switch (node_show)
        {
            case >= 1:
            {
                output_unit.Add("%");
                output_unit.Add("%  Draw filled dots at each node:");
                output_unit.Add("%");
                output_unit.Add("%  Set the color to blue:");
                output_unit.Add("%");
                output_unit.Add("0.000  0.150  0.750  setrgbcolor");
                output_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int)(
                        ((x_max - node_xy[0 + node * 2]) * x_ps_min
                         + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int)(
                        ((y_max - node_xy[1 + node * 2]) * y_ps_min
                         + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                        / (y_max - y_min));

                    output_unit.Add("newpath  "
                                    + x_ps + "  "
                                    + y_ps + "  "
                                    + circle_size + " 0 360 arc closepath fill");
                }

                break;
            }
        }

        switch (node_show)
        {
            //
            //  Label the nodes.
            //
            case >= 2:
            {
                output_unit.Add("%");
                output_unit.Add("%  Label the nodes:");
                output_unit.Add("%");
                output_unit.Add("%  Set the color to darker blue:");
                output_unit.Add("%");
                output_unit.Add("0.000  0.250  0.850  setrgbcolor");
                output_unit.Add("/Times-Roman findfont");
                output_unit.Add("0.20 inch scalefont");
                output_unit.Add("setfont");

                output_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int)(
                        ((x_max - node_xy[0 + node * 2]) * x_ps_min
                         + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int)(
                        ((y_max - node_xy[1 + node * 2]) * y_ps_min
                         + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                        / (y_max - y_min));

                    output_unit.Add("newpath  "
                                    + x_ps + "  "
                                    + (y_ps + 5) + "  moveto ("
                                    + node + ") show");
                }

                break;
            }
        }

        switch (element_show)
        {
            //
            //  Draw the elements.
            //
            case >= 1:
            {
                output_unit.Add("%");
                output_unit.Add("%  Set the RGB color to red.");
                output_unit.Add("%");
                output_unit.Add("0.900  0.200  0.100 setrgbcolor");
                output_unit.Add("%");
                output_unit.Add("%  Draw the elements.");
                output_unit.Add("%");

                for (element = 0; element < element_num; element++)
                {
                    output_unit.Add("newpath");

                    for (i = 0; i <= element_order; i++)
                    {
                        e = typeMethods.i4_wrap(i, 0, element_order - 1);

                        node = element_node[e + element * element_order];

                        x_ps = (int)(
                            ((x_max - node_xy[0 + node * 2]) * x_ps_min
                             + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                            / (x_max - x_min));

                        y_ps = (int)(
                            ((y_max - node_xy[1 + node * 2]) * y_ps_min
                             + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                            / (y_max - y_min));

                        switch (i)
                        {
                            case 0:
                                output_unit.Add(x_ps + "  " + y_ps + "  moveto");
                                break;
                            default:
                                output_unit.Add(x_ps + "  " + y_ps + "  lineto");
                                break;
                        }
                    }

                    output_unit.Add("stroke");
                }

                break;
            }
        }

        switch (element_show)
        {
            //
            //  Label the elements.
            //
            case >= 2:
            {
                output_unit.Add("%");
                output_unit.Add("%  Label the elements.");
                output_unit.Add("%");
                output_unit.Add("%  Set the RGB color to darker red.");
                output_unit.Add("%");
                output_unit.Add("0.950  0.250  0.150 setrgbcolor");
                output_unit.Add("/Times-Roman findfont");
                output_unit.Add("0.20 inch scalefont");
                output_unit.Add("setfont");
                output_unit.Add("%");

                for (element = 0; element < element_num; element++)
                {
                    ave_x = 0.0;
                    ave_y = 0.0;

                    for (i = 0; i < element_order; i++)
                    {
                        node = element_node[i + element * element_order];
                        ave_x += node_xy[0 + node * 2];
                        ave_y += node_xy[1 + node * 2];
                    }

                    ave_x /= element_order;
                    ave_y /= element_order;

                    x_ps = (int)(
                        ((x_max - ave_x) * x_ps_min
                         + (+ave_x - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int)(
                        ((y_max - ave_y) * y_ps_min
                         + (ave_y - y_min) * y_ps_max)
                        / (y_max - y_min));

                    output_unit.Add(x_ps + "  "
                                         + y_ps + "  moveto ("
                                         + element + ") show");
                }

                break;
            }
        }

        output_unit.Add("%");
        output_unit.Add("restore  showpage");
        output_unit.Add("%");
        output_unit.Add("%  End of page.");
        output_unit.Add("%");
        output_unit.Add("%%Trailer");
        output_unit.Add("%%EOF");

        try
        {
            File.WriteAllLines(output_filename, output_unit);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("PLOT_Q4_MESH - Fatal error!");
            Console.WriteLine("  Could not open the output EPS file.");
        }
    }
}