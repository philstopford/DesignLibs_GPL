using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.Types
{
    public class points
    {
        public static void points_plot(string file_name, int node_num, double[] node_xy,
                bool node_label)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_PLOT plots a pointset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 March 2005
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
            //    Input, double NODE_XY[2*NODE_NUM], the nodes.
            //
            //    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
            //
            //  Local parameters:
            //
            //    int CIRCLE_SIZE, controls the size of the circles depicting
            //    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
            //    barely visible.
            //
        {
            int circle_size = 3;
            int delta;
            List<string> file_unit = new List<string>();
            int node;
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
                delta = (int) ((double) (x_ps_max - x_ps_min)
                    * (y_scale - x_scale) / (2.0 * y_scale));

                x_ps_max = x_ps_max - delta;
                x_ps_min = x_ps_min + delta;

                x_ps_max_clip = x_ps_max_clip - delta;
                x_ps_min_clip = x_ps_min_clip + delta;

                x_scale = y_scale;
            }
            else if (y_scale < x_scale)
            {
                delta = (int) ((double) (y_ps_max - y_ps_min)
                    * (x_scale - y_scale) / (2.0 * x_scale));

                y_ps_max = y_ps_max - delta;
                y_ps_min = y_ps_min + delta;

                y_ps_max_clip = y_ps_max_clip - delta;
                y_ps_min_clip = y_ps_min_clip + delta;

                y_scale = x_scale;
            }

            file_unit.Add("%!PS-Adobe-3.0 EPSF-3.0");
            file_unit.Add("%%Creator: points_plot.C");
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

            //
            //  Label the nodes.
            //
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
                              + y_ps + 5 + "  moveto ("
                              + node + 1 + ") show");
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
                Console.WriteLine("POINTS_PLOT - Fatal error!");
                Console.WriteLine("  Could not open the output EPS file.");
            }
        }
    }
}