using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.Pyramid
{
    public static class Grid
    {
        public static int pyramid_grid_size(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_GRID_SIZE sizes a pyramid grid.
            //
            //  Discussion:
            //
            //    0:  x
            //
            //    1:  x  x
            //        x  x
            //
            //    2:  x  x  x
            //        x  x  x
            //        x  x  x
            //
            //    3:  x  x  x  x
            //        x  x  x  x
            //        x  x  x  x
            //        x  x  x  x
            //
            //    N  Size
            //
            //    0     1
            //    1     5 = 1 + 4
            //    2    14 = 1 + 4 + 9
            //    3    30 = 1 + 4 + 9 + 16
            //    4    55 = 1 + 4 + 9 + 16 + 25
            //    5    91 = 1 + 4 + 9 + 16 + 25 + 36
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of subintervals.
            //
            //    Output, int PYRAMID_GRID_SIZE, the number of
            //    nodes in the grid of size N.
            //
        {
            int np1;
            int value;

            np1 = n + 1;

            value = (np1 * (np1 + 1) * (2 * np1 + 1)) / 6;

            return value;
        }

        public static double[] pyramid_unit_grid(int n, int ng)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_UNIT_GRID computes grid points in the unit pyramid.
            //
            //  Discussion:
            //
            //    The unit pyramid has base (-1,-1,0), (+1,1,0), (+1,+1,0), (-1,+1,0)
            //    and vertex (0,0,1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of subintervals.
            //
            //    Input, int NG, the number of nodes to generate,
            //    as determined by pyramid_grid_size().
            //
            //    Output, double PYRAMID_UNIT_GRID[3*NG], the grid point coordinates.
            //
        {
            int g;
            int hi;
            int i;
            int j;
            int k;
            int lo;
            double[] pg;

            pg = new double[3 * ng];

            g = 0;

            for (k = n; 0 <= k; k--)
            {
                hi = n - k;
                lo = -hi;
                for (j = lo; j <= hi; j = j + 2)
                {
                    for (i = lo; i <= hi; i = i + 2)
                    {
                        pg[0 + g * 3] = (double)(i) / (double)(n);
                        pg[1 + g * 3] = (double)(j) / (double)(n);
                        pg[2 + g * 3] = (double)(k) / (double)(n);
                        g = g + 1;
                    }
                }
            }

            return pg;
        }

        public static void pyramid_unit_grid_plot(int n, int ng, double[] pg, string header)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_UNIT_GRID_PLOT sets up a GNUPLOT plot of a unit pyramid grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of subintervals.
            //
            //    Input, int NG, the number of nodes to generate,
            //    as determined by pyramid_grid_size().
            //
            //    Input, double PG[3*NG], the grid point coordinates.
            //
            //    Input, string HEADER, the header for the files.
            //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            int j;
            string node_filename;
            List<string> node_unit = new List<string>();
            string plot_filename;
            double[] v1 = new double[3];
            double[] v2 = new double[3];
            double[] v3 = new double[3];
            double[] v4 = new double[3];
            double[] v5 = new double[3];
            string vertex_filename;
            List<string> vertex_unit = new List<string>();
            //
            //  Create the vertex file.
            //
            pyramid_unit_vertices(ref v1, ref v2, ref v3, ref v4, ref v5);

            vertex_filename = header + "_vertices.txt";

            vertex_unit.Add(v2[0] + "  "
                                  + v2[1] + "  "
                                  + v2[2] + "");
            vertex_unit.Add(v3[0] + "  "
                                  + v3[1] + "  "
                                  + v3[2] + "");
            vertex_unit.Add(v4[0] + "  "
                                  + v4[1] + "  "
                                  + v4[2] + "");
            vertex_unit.Add(v5[0] + "  "
                                  + v5[1] + "  "
                                  + v5[2] + "");
            vertex_unit.Add(v2[0] + "  "
                                  + v2[1] + "  "
                                  + v2[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v1[0] + "  "
                                  + v1[1] + "  "
                                  + v1[2] + "");
            vertex_unit.Add(v2[0] + "  "
                                  + v2[1] + "  "
                                  + v2[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v1[0] + "  "
                                  + v1[1] + "  "
                                  + v1[2] + "");
            vertex_unit.Add(v3[0] + "  "
                                  + v3[1] + "  "
                                  + v3[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v1[0] + "  "
                                  + v1[1] + "  "
                                  + v1[2] + "");
            vertex_unit.Add(v4[0] + "  "
                                  + v4[1] + "  "
                                  + v4[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v1[0] + "  "
                                  + v1[1] + "  "
                                  + v1[2] + "");
            vertex_unit.Add(v5[0] + "  "
                                  + v5[1] + "  "
                                  + v5[2] + "");
            vertex_unit.Add("");

            File.WriteAllLines(vertex_filename, vertex_unit);

            Console.WriteLine("");
            Console.WriteLine("  Created vertex file '" + vertex_filename + "'");
            //
            //  Create the node file.
            //
            node_filename = header + "_nodes.txt";

            for (j = 0; j < ng; j++)
            {
                node_unit.Add(pg[0 + j * 3] + "  "
                                            + pg[1 + j * 3] + "  "
                                            + pg[2 + j * 3] + "");
            }

            File.WriteAllLines(node_filename, node_unit);

            Console.WriteLine(" Created node file '" + node_filename + "'");
            //
            //  Create the command file.
            //
            command_filename = header + "_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");

            plot_filename = header + ".png";

            command_unit.Add("set output '" + plot_filename + "'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + header + "'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("set view equal xyz");
            command_unit.Add("set view 80, 40");
            command_unit.Add("set style data lines");
            command_unit.Add("set timestamp");
            command_unit.Add("splot '" + vertex_filename + "' with lines lw 3, \\");
            command_unit.Add("      '" + node_filename + "' with points pt 7 lt 0");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        public static void pyramid_unit_vertices(ref double[] v1, ref double[] v2, ref double[] v3,
                ref double[] v4, ref double[] v5)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_UNIT_VERTICES returns the vertices of the unit pyramid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], the vertices.
            //
        {
            v1[0] = 0.0;
            v1[1] = 0.0;
            v1[2] = +1.0;

            v2[0] = -1.0;
            v2[1] = -1.0;
            v2[2] = 0.0;

            v3[0] = +1.0;
            v3[1] = -1.0;
            v3[2] = 0.0;

            v4[0] = +1.0;
            v4[1] = +1.0;
            v4[2] = 0.0;

            v5[0] = -1.0;
            v5[1] = +1.0;
            v5[2] = 0.0;
        }
    }
}