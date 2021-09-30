using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Plot
    {
        public static void triangulation_order6_plot(string file_name, int node_num,
            double[] node_xy, int triangle_num, int[] triangle_node, int node_show,
        int triangle_show )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
        //
        //  Discussion:
        //
        //    The triangulation is most usually a Delaunay triangulation,
        //    but this is not necessary.
        //
        //    This routine has been specialized to deal correctly ONLY with
        //    a mesh of 6 node elements, with the property that starting
        //    from local node 1 and traversing the edges of the element will
        //    result in encountering local nodes 1, 4, 2, 5, 3, 6 in that order.
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
        //    Input, string FILE_NAME, the name of the file to create.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists, for each triangle,
        //    the indices of the nodes that form the vertices and midsides
        //    of the triangle.
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
            List<string> file_unit = new List<string>();
            int i;
            int ip1;
            int node;
            int[] order =  {
                1, 4, 2, 5, 3, 6
            }
            ;
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

            x_max = x_max + 0.05 * x_scale;
            x_min = x_min - 0.05 * x_scale;
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

            y_max = y_max + 0.05 * y_scale;
            y_min = y_min - 0.05 * y_scale;
            y_scale = y_max - y_min;

            if (x_scale < y_scale)
            {
                delta = (int)((double) (x_ps_max - x_ps_min)
                    * (y_scale - x_scale) / (2.0 * y_scale));

                x_ps_max = x_ps_max - delta;
                x_ps_min = x_ps_min + delta;

                x_ps_max_clip = x_ps_max_clip - delta;
                x_ps_min_clip = x_ps_min_clip + delta;

                x_scale = y_scale;
            }
            else if (y_scale < x_scale)
            {
                delta = (int)((double) (y_ps_max - y_ps_min)
                    * (x_scale - y_scale) / (2.0 * x_scale));

                y_ps_max = y_ps_max - delta;
                y_ps_min = y_ps_min + delta;

                y_ps_max_clip = y_ps_max_clip - delta;
                y_ps_min_clip = y_ps_min_clip + delta;

                y_scale = x_scale;
            }

            file_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
            file_unit.Add("%%Creator: triangulation_order6_plot.C");
            file_unit.Add("%%Title: " + file_name + "");

            file_unit.Add("%%Pages: 1");
            file_unit.Add("%%BoundingBox:  "
                      + x_ps_min + "  "
                      + y_ps_min + "  "
                      + x_ps_max + "  "
                      + y_ps_max + "");
            file_unit.Add("%%Document-Fonts: Times-Roman");
            file_unit.Add("%%LanguageLevel: 1");
            file_unit.Add("%%EndComments");
            file_unit.Add("%%BeginProlog");
            file_unit.Add("/inch {72 mul} def");
            file_unit.Add("%%EndProlog");
            file_unit.Add("%%Page:      1     1");
            file_unit.Add("save");
            file_unit.Add("%");
            file_unit.Add("%  Increase line width from default 0.");
            file_unit.Add("%");
            file_unit.Add("2 setlinewidth");
            file_unit.Add("%");
            file_unit.Add("% Set the RGB line color to very light gray.");
            file_unit.Add("%");
            file_unit.Add(" 0.9000 0.9000 0.9000 setrgbcolor");
            file_unit.Add("%");
            file_unit.Add("% Draw a gray border around the page.");
            file_unit.Add("%");
            file_unit.Add("newpath");
            file_unit.Add(x_ps_min + "  "
                + y_ps_min + "  moveto");
            file_unit.Add(x_ps_max + "  "
                + y_ps_min + "  lineto");
            file_unit.Add(x_ps_max + "  "
                + y_ps_max + "  lineto");
            file_unit.Add(x_ps_min + "  "
                + y_ps_max + "  lineto");
            file_unit.Add(x_ps_min + "  "
                + y_ps_min + "  lineto");
            file_unit.Add("stroke");
            file_unit.Add("%");
            file_unit.Add("% Set RGB line color to black.");
            file_unit.Add("%");
            file_unit.Add(" 0.0000 0.0000 0.0000 setrgbcolor");
            file_unit.Add("%");
            file_unit.Add("%  Set the font and its size:");
            file_unit.Add("%");
            file_unit.Add("/Times-Roman findfont");
            file_unit.Add("0.50 inch scalefont");
            file_unit.Add("setfont");
            file_unit.Add("%");
            file_unit.Add("%  Print a title:");
            file_unit.Add("%");
            file_unit.Add("%  210  702 moveto");
            file_unit.Add("%(Pointset) show");
            file_unit.Add("%");
            file_unit.Add("% Define a clipping polygon");
            file_unit.Add("%");
            file_unit.Add("newpath");
            file_unit.Add(x_ps_min_clip + "  "
                + y_ps_min_clip + "  moveto");
            file_unit.Add(x_ps_max_clip + "  "
                + y_ps_min_clip + "  lineto");
            file_unit.Add(x_ps_max_clip + "  "
                + y_ps_max_clip + "  lineto");
            file_unit.Add(x_ps_min_clip + "  "
                + y_ps_max_clip + "  lineto");
            file_unit.Add(x_ps_min_clip + "  "
                + y_ps_min_clip + "  lineto");
            file_unit.Add("clip newpath");
            //
            //  Draw the nodes.
            //
            if (node_num <= 200)
            {
                circle_size = 5;
            }
            else if (node_num <= 500)
            {
                circle_size = 4;
            }
            else if (node_num <= 1000)
            {
                circle_size = 3;
            }
            else if (node_num <= 5000)
            {
                circle_size = 2;
            }
            else
            {
                circle_size = 1;
            }

            if (1 <= node_show)
            {
                file_unit.Add("%");
                file_unit.Add("%  Draw filled dots at each node:");
                file_unit.Add("%");
                file_unit.Add("%  Set the color to blue:");
                file_unit.Add("%");
                file_unit.Add("0.000  0.150  0.750  setrgbcolor");
                file_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int) (
                        ((x_max - node_xy[0 + node * 2]) * (double) (x_ps_min)
                         + (+node_xy[0 + node * 2] - x_min) * (double) (x_ps_max))
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - node_xy[1 + node * 2]) * (double) (y_ps_min)
                         + (node_xy[1 + node * 2] - y_min) * (double) (y_ps_max))
                        / (y_max - y_min));

                    file_unit.Add("newpath  "
                              + x_ps + "  "
                              + y_ps + "  "
                              + circle_size + " 0 360 arc closepath fill");
                }
            }

            //
            //  Label the nodes.
            //
            if (2 <= node_show)
            {
                file_unit.Add("%");
                file_unit.Add("%  Label the nodes:");
                file_unit.Add("%");
                file_unit.Add("%  Set the color to darker blue:");
                file_unit.Add("%");
                file_unit.Add("0.000  0.250  0.850  setrgbcolor");
                file_unit.Add("/Times-Roman findfont");
                file_unit.Add("0.20 inch scalefont");
                file_unit.Add("setfont");

                file_unit.Add("%");

                for (node = 0; node < node_num; node++)
                {
                    x_ps = (int) (
                        ((x_max - node_xy[0 + node * 2]) * (double) (x_ps_min)
                         + (+node_xy[0 + node * 2] - x_min) * (double) (x_ps_max))
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - node_xy[1 + node * 2]) * (double) (y_ps_min)
                         + (node_xy[1 + node * 2] - y_min) * (double) (y_ps_max))
                        / (y_max - y_min));

                    file_unit.Add("newpath  "
                              + x_ps + "  "
                              + (y_ps + 5) + "  moveto ("
                              + (node + 1) + ") show");
                }
            }

            //
            //  Draw the triangles.
            //
            if (1 <= triangle_show)
            {
                file_unit.Add("%");
                file_unit.Add("%  Set the RGB color to red.");
                file_unit.Add("%");
                file_unit.Add("0.900  0.200  0.100 setrgbcolor");
                file_unit.Add("%");
                file_unit.Add("%  Draw the triangles.");
                file_unit.Add("%");

                for (triangle = 0; triangle < triangle_num; triangle++)
                {
                    node = triangle_node[order[0] - 1 + triangle * 6] - 1;

                    x_ps = (int) (
                        ((x_max - node_xy[0 + node * 2]) * (double) (x_ps_min)
                         + (+node_xy[0 + node * 2] - x_min) * (double) (x_ps_max))
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - node_xy[1 + node * 2]) * (double) (y_ps_min)
                         + (node_xy[1 + node * 2] - y_min) * (double) (y_ps_max))
                        / (y_max - y_min));

                    file_unit.Add("newpath  " + x_ps + "  " + y_ps + "  moveto");

                    for (i = 1; i <= 6; i++)
                    {
                        ip1 = (i % 6) + 1;
                        node = triangle_node[order[ip1 - 1] - 1 + triangle * 6] - 1;

                        x_ps = (int) (
                            ((x_max - node_xy[0 + node * 2]) * (double) (x_ps_min)
                             + (+node_xy[0 + node * 2] - x_min) * (double) (x_ps_max))
                            / (x_max - x_min));

                        y_ps = (int) (
                            ((y_max - node_xy[1 + node * 2]) * (double) (y_ps_min)
                             + (node_xy[1 + node * 2] - y_min) * (double) (y_ps_max))
                            / (y_max - y_min));

                        file_unit.Add(x_ps + "  " + y_ps + "  lineto");
                    }

                    file_unit.Add("stroke");
                }
            }

            //
            //  Label the triangles.
            //
            if (2 <= triangle_show)
            {
                file_unit.Add("%");
                file_unit.Add("%  Label the triangles.");
                file_unit.Add("%");
                file_unit.Add("%  Set the RGB color to darker red.");
                file_unit.Add("%");
                file_unit.Add("0.950  0.250  0.150 setrgbcolor");
                file_unit.Add("/Times-Roman findfont");
                file_unit.Add("0.20 inch scalefont");
                file_unit.Add("setfont");
                file_unit.Add("%");

                for (triangle = 0; triangle < triangle_num; triangle++)
                {
                    ave_x = 0.0;
                    ave_y = 0.0;

                    for (i = 0; i < 6; i++)
                    {
                        node = triangle_node[i + triangle * 6] - 1;
                        ave_x = ave_x + node_xy[0 + node * 2];
                        ave_y = ave_y + node_xy[1 + node * 2];
                    }

                    ave_x = ave_x / 6.0;
                    ave_y = ave_y / 6.0;

                    x_ps = (int) (
                        ((x_max - ave_x) * (double) (x_ps_min)
                         + (+ave_x - x_min) * (double) (x_ps_max))
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - ave_y) * (double) (y_ps_min)
                         + (ave_y - y_min) * (double) (y_ps_max))
                        / (y_max - y_min));

                    file_unit.Add(x_ps.ToString().PadLeft(4) + "  "
                        + y_ps.ToString().PadLeft(4) + "  "
                        + "moveto (" + (triangle + 1) + ") show");
                }
            }

            file_unit.Add("%");
            file_unit.Add("restore showpage");
            file_unit.Add("%");
            file_unit.Add("% End of page");
            file_unit.Add("%");
            file_unit.Add("%%Trailer");
            file_unit.Add("%%EOF");

            try
            {
                File.WriteAllLines(file_name, file_unit);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_ORDER6_PLOT - Fatal error!");
                Console.WriteLine("  Could not open the output EPS file.");
            }
        }
    }
}