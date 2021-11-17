using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Plot
{
    public static void triangulation_order4_plot(string plot_filename, int node_num,
            double[] node_xy, int triangle_num, int[] triangle_node, int node_show,
            int triangle_show )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER4_PLOT plots a 4-node triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string PLOT_FILENAME, the name of the output file.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[4*TRIANGLE_NUM], lists, for each triangle,
        //    the indices of the nodes that form the vertices of the triangle,
        //    and the centroid.
        //
        //    Input, int NODE_SHOW:
        //    0, do not show nodes;
        //    1, show nodes;
        //    2, show nodes and label them.
        //
        //    Input, int TRIANGLE_SHOW:
        //    0, do not show triangles;
        //    1, show triangles;
        //    2, show triangles and label them.
        //
    {
        double ave_x;
        double ave_y;
        int circle_size;
        int delta;
        int e;
        List<string> plot_unit = new();
        int i;
        int node;
        int triangle;
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
            delta = (int)((x_ps_max - x_ps_min)
                * (y_scale - x_scale) / (2.0 * y_scale));

            x_ps_max -= delta;
            x_ps_min += delta;

            x_ps_max_clip -= delta;
            x_ps_min_clip += delta;

            x_scale = y_scale;
        }
        else if (y_scale < x_scale)
        {
            delta = (int)((y_ps_max - y_ps_min)
                * (x_scale - y_scale) / (2.0 * x_scale));

            y_ps_max -= delta;
            y_ps_min += delta;

            y_ps_max_clip -= delta;
            y_ps_min_clip += delta;

            y_scale = x_scale;
        }

        plot_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
        plot_unit.Add("%%Creator: triangulation_order4_plot.C");
        plot_unit.Add("%%Title: " + plot_filename + "");

        plot_unit.Add("%%Pages: 1");
        plot_unit.Add("%%BoundingBox:  "
                      + x_ps_min + "  "
                      + y_ps_min + "  "
                      + x_ps_max + "  "
                      + y_ps_max + "");
        plot_unit.Add("%%Document-Fonts: Times-Roman");
        plot_unit.Add("%%LanguageLevel: 1");
        plot_unit.Add("%%EndComments");
        plot_unit.Add("%%BeginProlog");
        plot_unit.Add("/inch {72 mul} def");
        plot_unit.Add("%%EndProlog");
        plot_unit.Add("%%Page:      1     1");
        plot_unit.Add("save");
        plot_unit.Add("%");
        plot_unit.Add("%  Increase line width from default 0.");
        plot_unit.Add("%");
        plot_unit.Add("2 setlinewidth");
        plot_unit.Add("%");
        plot_unit.Add("% Set the RGB line color to very light gray.");
        plot_unit.Add("%");
        plot_unit.Add(" 0.9000 0.9000 0.9000 setrgbcolor");
        plot_unit.Add("%");
        plot_unit.Add("% Draw a gray border around the page.");
        plot_unit.Add("%");
        plot_unit.Add("newpath");
        plot_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  moveto");
        plot_unit.Add(x_ps_max + "  "
                               + y_ps_min + "  lineto");
        plot_unit.Add(x_ps_max + "  "
                               + y_ps_max + "  lineto");
        plot_unit.Add(x_ps_min + "  "
                               + y_ps_max + "  lineto");
        plot_unit.Add(x_ps_min + "  "
                               + y_ps_min + "  lineto");
        plot_unit.Add("stroke");
        plot_unit.Add("%");
        plot_unit.Add("% Set RGB line color to black.");
        plot_unit.Add("%");
        plot_unit.Add(" 0.0000 0.0000 0.0000 setrgbcolor");
        plot_unit.Add("%");
        plot_unit.Add("%  Set the font and its size:");
        plot_unit.Add("%");
        plot_unit.Add("/Times-Roman findfont");
        plot_unit.Add("0.50 inch scalefont");
        plot_unit.Add("setfont");
        plot_unit.Add("%");
        plot_unit.Add("%  Print a title:");
        plot_unit.Add("%");
        plot_unit.Add("%  210  702 moveto");
        plot_unit.Add("%(Pointset) show");
        plot_unit.Add("%");
        plot_unit.Add("% Define a clipping polygon");
        plot_unit.Add("%");
        plot_unit.Add("newpath");
        plot_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  moveto");
        plot_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        plot_unit.Add(x_ps_max_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        plot_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_max_clip + "  lineto");
        plot_unit.Add(x_ps_min_clip + "  "
                                    + y_ps_min_clip + "  lineto");
        plot_unit.Add("clip newpath");
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
                plot_unit.Add("%");
                plot_unit.Add("%  Draw filled dots at each node:");
                plot_unit.Add("%");
                plot_unit.Add("%  Set the color to blue:");
                plot_unit.Add("%");
                plot_unit.Add("0.000  0.150  0.750  setrgbcolor");
                plot_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int) (
                        ((x_max - node_xy[0 + node * 2]) * x_ps_min
                         + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - node_xy[1 + node * 2]) * y_ps_min
                         + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                        / (y_max - y_min));

                    plot_unit.Add("newpath  "
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
                plot_unit.Add("%");
                plot_unit.Add("%  Label the nodes:");
                plot_unit.Add("%");
                plot_unit.Add("%  Set the color to darker blue:");
                plot_unit.Add("%");
                plot_unit.Add("0.000  0.250  0.850  setrgbcolor");
                plot_unit.Add("/Times-Roman findfont");
                plot_unit.Add("0.20 inch scalefont");
                plot_unit.Add("setfont");

                plot_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int) (
                        ((x_max - node_xy[0 + node * 2]) * x_ps_min
                         + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - node_xy[1 + node * 2]) * y_ps_min
                         + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                        / (y_max - y_min));

                    plot_unit.Add("newpath  "
                                  + x_ps + "  "
                                  + (y_ps + 5) + "  moveto ("
                                  + (node + 1) + ") show");
                }

                break;
            }
        }

        switch (triangle_show)
        {
            //
            //  Draw the triangles.
            //
            case >= 1:
            {
                plot_unit.Add("%");
                plot_unit.Add("%  Set the RGB color to red.");
                plot_unit.Add("%");
                plot_unit.Add("0.900  0.200  0.100 setrgbcolor");
                plot_unit.Add("%");
                plot_unit.Add("%  Draw the triangles.");
                plot_unit.Add("%");

                for (triangle = 0; triangle < triangle_num; triangle++)
                {
                    plot_unit.Add("newpath");

                    for (i = 1; i <= 4; i++)
                    {
                        e = typeMethods.i4_wrap(i, 1, 3);

                        node = triangle_node[e - 1 + triangle * 4] - 1;

                        x_ps = (int) (
                            ((x_max - node_xy[0 + node * 2]) * x_ps_min
                             + (+node_xy[0 + node * 2] - x_min) * x_ps_max)
                            / (x_max - x_min));

                        y_ps = (int) (
                            ((y_max - node_xy[1 + node * 2]) * y_ps_min
                             + (node_xy[1 + node * 2] - y_min) * y_ps_max)
                            / (y_max - y_min));

                        switch (i)
                        {
                            case 1:
                                plot_unit.Add(x_ps + "  " + y_ps + "  moveto");
                                break;
                            default:
                                plot_unit.Add(x_ps + "  " + y_ps + "  lineto");
                                break;
                        }
                    }

                    plot_unit.Add("stroke");
                }

                break;
            }
        }

        switch (triangle_show)
        {
            //
            //  Label the triangles.
            //
            case >= 2:
            {
                plot_unit.Add("%");
                plot_unit.Add("%  Label the triangles.");
                plot_unit.Add("%");
                plot_unit.Add("%  Set the RGB color to darker red.");
                plot_unit.Add("%");
                plot_unit.Add("0.950  0.250  0.150 setrgbcolor");
                plot_unit.Add("/Times-Roman findfont");
                plot_unit.Add("0.20 inch scalefont");
                plot_unit.Add("setfont");
                plot_unit.Add("%");

                for (triangle = 0; triangle < triangle_num; triangle++)
                {
                    ave_x = 0.0;
                    ave_y = 0.0;

                    for (i = 1; i <= 3; i++)
                    {
                        node = triangle_node[i - 1 + triangle * 4] - 1;
                        ave_x += node_xy[0 + node * 2];
                        ave_y += node_xy[1 + node * 2];
                    }

                    ave_x /= 3.0;
                    ave_y /= 3.0;

                    x_ps = (int) (
                        ((x_max - ave_x) * x_ps_min
                         + (+ave_x - x_min) * x_ps_max)
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - ave_y) * y_ps_min
                         + (ave_y - y_min) * y_ps_max)
                        / (y_max - y_min));

                    plot_unit.Add(x_ps + "  "
                                       + y_ps + "  moveto ("
                                       + (triangle + 1) + ") show");
                }

                break;
            }
        }

        plot_unit.Add("%");
        plot_unit.Add("restore  showpage");
        plot_unit.Add("%");
        plot_unit.Add("%  End of page.");
        plot_unit.Add("%");
        plot_unit.Add("%%Trailer");
        plot_unit.Add("%%EOF");


        try
        {
            File.WriteAllLines(plot_filename, plot_unit);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_ORDER4_PLOT - Fatal error!");
            Console.WriteLine("  Could not open the output EPS file.");
        }
    }
}