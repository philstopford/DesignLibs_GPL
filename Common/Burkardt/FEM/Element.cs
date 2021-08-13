using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.OrderNS;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class Element
    {
        public static string element_code(int i)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELEMENT_CODE returns the code for each element.
            //
            //  List:
            //
            //    I  ELEMENT_CODE   Definition
            //    -  ------------   ----------
            //    1  Q4             4 node linear Lagrange/serendipity quadrilateral;
            //    2  Q8             8 node quadratic serendipity quadrilateral;
            //    3  Q9             9 node quadratic Lagrange quadrilateral;
            //    4  Q12            12 node cubic serendipity quadrilateral;
            //    5  Q16            16 node cubic Lagrange quadrilateral;
            //    6  QL             6 node linear/quadratic quadrilateral;
            //    7  T3             3 node linear triangle;
            //    8  T4             4 node cubic bubble triangle
            //    9  T6             6 node quadratic triangle;
            //   10  T10            10 node cubic triangle.
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, the number of the element.  
            //
            //    Output, string ELEMENT_CODE, the code for the element.
            //
        {
            string value;

            if (i == 1)
            {
                value = "Q4";
            }
            else if (i == 2)
            {
                value = "Q8";
            }
            else if (i == 3)
            {
                value = "Q9";
            }
            else if (i == 4)
            {
                value = "Q12";
            }
            else if (i == 5)
            {
                value = "Q16";
            }
            else if (i == 6)
            {
                value = "QL";
            }
            else if (i == 7)
            {
                value = "T3";
            }
            else if (i == 8)
            {
                value = "T4";
            }
            else if (i == 9)
            {
                value = "T6";
            }
            else if (i == 10)
            {
                value = "T10";
            }
            else
            {
                value = "????";
            }

            return value;
        }

        public static void elements_eps(string file_name, int node_num, double[] node_xy, string code,
                int element_num, bool[] element_mask, int[] element_node, int node_show,
                int element_show)

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    ELEMENTS_EPS creates an EPS file image of the elements of a grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 April 2005
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
            //    Input, string CODE, the code for the element.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, bool ELEMENT_MASK[ELEMENT_NUM], a mask for the elements.
            //    Only elements with a TRUE mask will be shown.
            //
            //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes making up
            //    each element.
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
            int circle_size = 3;
            int delta;
            int element;
            int element_order;
            List<string> file_unit = new List<string>();
            int local;
            int node;
            bool[] node_mask;
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

            element_order = Order.order_code(code);
            //
            //  Determine which nodes are visible, controlled by which elements are visible.
            //
            node_mask = new bool[node_num];
            for (node = 0; node < node_num; node++)
            {
                node_mask[node] = false;
            }

            for (element = 0; element < element_num; element++)
            {
                if (element_mask[element])
                {
                    for (local = 0; local < element_order; local++)
                    {
                        node = element_node[local + element * element_order] - 1;
                        node_mask[node] = true;
                    }
                }
            }

            //
            //  We need to do some figuring here, so that we can determine
            //  the range of the data, and hence the height and width
            //  of the piece of paper.
            //
            x_max = -typeMethods.r8_huge();
            for (node = 0; node < node_num; node++)
            {
                if (node_mask[node])
                {
                    if (x_max < node_xy[0 + node * 2])
                    {
                        x_max = node_xy[0 + node * 2];
                    }
                }
            }

            x_min = typeMethods.r8_huge();
            for (node = 0; node < node_num; node++)
            {
                if (node_mask[node])
                {
                    if (node_xy[0 + node * 2] < x_min)
                    {
                        x_min = node_xy[0 + node * 2];
                    }
                }
            }

            x_scale = x_max - x_min;

            x_max = x_max + 0.05 * x_scale;
            x_min = x_min - 0.05 * x_scale;
            x_scale = x_max - x_min;

            y_max = -typeMethods.r8_huge();
            for (node = 0; node < node_num; node++)
            {
                if (node_mask[node])
                {
                    if (y_max < node_xy[1 + node * 2])
                    {
                        y_max = node_xy[1 + node * 2];
                    }
                }
            }

            y_min = typeMethods.r8_huge();
            for (node = 0; node < node_num; node++)
            {
                if (node_mask[node])
                {
                    if (node_xy[1 + node * 2] < y_min)
                    {
                        y_min = node_xy[1 + node * 2];
                    }
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
            file_unit.Add("%%Creator: elements_eps.C");
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
                    if (node_mask[node])
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
                    if (node_mask[node])
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
                }
            }

            //
            //  Draw the elements.
            //
            file_unit.Add("%");
            file_unit.Add("%  Draw the element sides:");
            file_unit.Add("%");
            file_unit.Add(" 9.0000 0.0000 0.0000 setrgbcolor");

            for (element = 0; element < element_num; element++)
            {
                if (element_mask[element])
                {
                    local = 1;
                    node = element_node[local - 1 + element * element_order] - 1;

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
                                  + y_ps + "  moveto");

                    for (;;)
                    {
                        local = NextBoundaryNode.next_boundary_node(local, code);
                        node = element_node[local - 1 + element * element_order] - 1;

                        x_ps = (int) (
                            ((x_max - node_xy[0 + node * 2]) * (double) (x_ps_min)
                             + (+node_xy[0 + node * 2] - x_min) * (double) (x_ps_max))
                            / (x_max - x_min));

                        y_ps = (int) (
                            ((y_max - node_xy[1 + node * 2]) * (double) (y_ps_min)
                             + (node_xy[1 + node * 2] - y_min) * (double) (y_ps_max))
                            / (y_max - y_min));

                        file_unit.Add("  "
                                      + x_ps + "  "
                                      + y_ps + "  lineto");

                        if (local == 1)
                        {
                            break;
                        }
                    }

                    file_unit.Add("stroke");
                }
            }

            //
            //  Label the elements.
            //
            file_unit.Add("%");
            file_unit.Add("%  Label the elements:");
            file_unit.Add("%");
            file_unit.Add(" 1.0000 0.0000 0.0000 setrgbcolor");
            file_unit.Add("/Times-Roman findfont");
            file_unit.Add("0.30 inch scalefont setfont");

            for (element = 0; element < element_num; element++)
            {
                if (element_mask[element])
                {
                    ave_x = 0.0;
                    ave_y = 0.0;

                    for (local = 0; local < element_order; local++)
                    {
                        node = element_node[local + element_order * element] - 1;

                        ave_x = ave_x + node_xy[0 + node * 2];
                        ave_y = ave_y + node_xy[1 + node * 2];
                    }

                    ave_x = ave_x / (double) (element_order);
                    ave_y = ave_y / (double) (element_order);

                    x_ps = (int) (
                        ((x_max - ave_x) * (double) (x_ps_min)
                         + (+ave_x - x_min) * (double) (x_ps_max))
                        / (x_max - x_min));

                    y_ps = (int) (
                        ((y_max - ave_y) * (double) (y_ps_min)
                         + (ave_y - y_min) * (double) (y_ps_max))
                        / (y_max - y_min));

                    file_unit.Add("newpath  "
                                  + x_ps + "  "
                                  + y_ps + 5 + "  moveto ("
                                  + element + 1 + ") show");
                }
            }

            //
            //  Finish up.
            //
            file_unit.Add("%");
            file_unit.Add("restore  showpage");
            file_unit.Add("%");
            file_unit.Add("%  End of page.");
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
                Console.WriteLine("ELEMENTS_EPS - Fatal error!");
                Console.WriteLine("  Could not open the output EPS file.");
            }

        }

        public static void element_data_read(string element_file, int element_num, int element_order,
                int element_att_num, ref int[] element_node, ref double[] element_att)

            //*****************************************************************************80
            //
            //  Purpose:
            //
            //    ELEMENT_DATA_READ reads the header information from an element file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string ELEMENT_FILE, the name of the file to be read.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Input, int ELEMENT_ATT_NUM, number of element attributes listed on each 
            //    node record.
            //
            //    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
            //    nodes that make up each element.
            //
            //    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
            //    of each element.
            //
        {
            int element;
            int i;
            int i1;
            int i2;
            int i3;
            string[] inputlines;
            int ival;
            double value;

            element = -1;

            try
            {
                inputlines = File.ReadAllLines(element_file);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("ELEMENT_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open file.");
                return;
            }

            foreach (string input in inputlines)
            {
                /*
                Read, but ignore, dimension line.
                */
                string[] tokens = input.Split(new[] {' '});
                if (element == -1)
                {
                    i1 = Convert.ToInt32(tokens[0]);
                    i2 = Convert.ToInt32(tokens[1]);
                    i3 = Convert.ToInt32(tokens[2]);
                }
                else
                {
                    int index = 0;
                    ival = Convert.ToInt32(tokens[index]);
                    index++;

                    for (i = 0; i < element_order; i++)
                    {
                        ival = Convert.ToInt32(tokens[index]);
                        index++;
                        element_node[i + element * element_order] = ival;
                    }

                    for (i = 0; i < element_att_num; i++)
                    {
                        value = Convert.ToInt32(tokens[index]);
                        index++;
                        element_att[i + element * element_att_num] = value;
                    }
                }

                element = element + 1;

                if (element_num <= element)
                {
                    break;
                }
            }

        }

        public static void element_size_read(string element_file, ref int element_num,
                ref int element_order, ref int element_att_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELEMENT_SIZE_READ reads the header information from an element file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string ELEMENT_FILE, the name of the file to be read.
            //
            //    Output, int *ELEMENT_NUM, the number of elements.
            //
            //    Output, int *ELEMENT_ORDER, the order of the elements.
            //
            //    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
            //
        {
            string input;

            try
            {
                input = File.ReadAllLines(element_file)[0];
                string[] tokens = input.Split(new[] {' '});
                element_num = Convert.ToInt32(tokens[0]);
                element_order = Convert.ToInt32(tokens[1]);
                element_att_num = Convert.ToInt32(tokens[2]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("ELEMENT_SIZE_READ - Fatal error!");
                Console.WriteLine("  Could not open file.");
            }


        }

    }
}