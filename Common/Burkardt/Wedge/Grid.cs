using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.Wedge
{
    public static class Grid
    {
        public static double[] wedge_grid(int n, int ng)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEDGE_GRID computes grid points in the unit wedge in 3D.
            //
            //  Discussion:
            //
            //    The interior of the unit wedge in 3D is defined by the constraints:
            //      0 <= X
            //      0 <= Y
            //           X + Y <= 1
            //     -1 <= Z <= +1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of subintervals.
            //    0 <= N.
            //
            //    Input, int NG, the number of grid points.
            //    This can be computed by WEDGE_GRID_SIZE, or else determined by
            //    NG =(N+1)*((N+1)*(N+2))/2.
            //
            //    Output, double WEDGE+GRID[3*NG], the coordinates of the grid points.
            //
        {
            int i;
            double ir;
            int j;
            double jr;
            int k;
            double kr;
            double nr;
            int p;
            double[] g;

            g = new double[3 * ng];

            if (n == 0)
            {
                g[0 + 0 * 3] = 0.5;
                g[1 + 0 * 3] = 0.5;
                g[2 + 0 * 3] = 0.0;
                return g;
            }

            p = 0;
            nr = (double) (n);

            for (k = 0; k <= n; k++)
            {
                kr = (double) (2 * k - n) / nr;
                for (j = 0; j <= n; j++)
                {
                    jr = (double) (j) / nr;
                    for (i = 0; i <= n - j; i++)
                    {
                        ir = (double) (i) / nr;
                        g[0 + p * 3] = ir;
                        g[1 + p * 3] = jr;
                        g[2 + p * 3] = kr;
                        p = p + 1;
                    }
                }
            }

            return g;
        }

        public static int wedge_grid_size(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEDGE_GRID_SIZE counts the points in a grid of the unit wedge in 3D.
            //
            //  Discussion:
            //
            //    The interior of the unit wedge in 3D is defined by the constraints:
            //      0 <= X
            //      0 <= Y
            //           X + Y <= 1
            //     -1 <= Z <= +1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of subintervals.
            //    0 <= N.
            //
            //    Output, int WEDGE_GRID_SIZE, the number of grid points.
            //
        {
            int ng;

            ng = (n + 1) * ((n + 1) * (n + 2)) / 2;

            return ng;
        }

        public static void wedge_grid_plot(int n, int ng, double[] g, string header )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEDGE_GRID_PLOT sets up a GNUPLOT plot of a unit wedge grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, int NG, the number of nodes.
        //
        //    Input, double G[3*NG], the grid point coordinates.
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
            double[] v6 = new double[3];
            string vertex_filename;
            List<string> vertex_unit = new List<string>();
            //
            //  Create the vertex file.
            //
            wedge_vertices(ref v1, ref v2, ref v3, ref v4, ref v5, ref v6);

            vertex_filename = header + "_vertices.txt";

            vertex_unit.Add(v1[0] + "  "
                + v1[1] + "  "
                + v1[2] + "");
            vertex_unit.Add(v2[0] + "  "
                + v2[1] + "  "
                + v2[2] + "");
            vertex_unit.Add(v3[0] + "  "
                + v3[1] + "  "
                + v3[2] + "");
            vertex_unit.Add(v1[0] + "  "
                + v1[1] + "  "
                + v1[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v4[0] + "  "
                + v4[1] + "  "
                + v4[2] + "");
            vertex_unit.Add(v5[0] + "  "
                + v5[1] + "  "
                + v5[2] + "");
            vertex_unit.Add(v6[0] + "  "
                + v6[1] + "  "
                + v6[2] + "");
            vertex_unit.Add(v4[0] + "  "
                + v4[1] + "  "
                + v4[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v1[0] + "  "
                + v1[1] + "  "
                + v1[2] + "");
            vertex_unit.Add(v4[0] + "  "
                + v4[1] + "  "
                + v4[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v2[0] + "  "
                + v2[1] + "  "
                + v2[2] + "");
            vertex_unit.Add(v5[0] + "  "
                + v5[1] + "  "
                + v5[2] + "");
            vertex_unit.Add("");

            vertex_unit.Add(v3[0] + "  "
                + v3[1] + "  "
                + v3[2] + "");
            vertex_unit.Add(v6[0] + "  "
                + v6[1] + "  "
                + v6[2] + "");
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
                node_unit.Add(g[0 + j * 3] + "  "
                    + g[1 + j * 3] + "  "
                    + g[2 + j * 3] + "");
            }

            File.WriteAllLines(node_filename, node_unit);
            Console.WriteLine("  Created node file '" + node_filename + "'");
            //
            //  Create the command file.
            //
            command_filename = header + "_commands.txt";

            command_unit.Clear();
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");

            plot_filename = header + ".png";

            command_unit.Add("set output '" + plot_filename + "");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set zlabel '<--- Z --->'");
            command_unit.Add("set title '" + header + "'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("#set view equal xyz");
            command_unit.Add("set view 80, 85");
            command_unit.Add("set style data lines");
            command_unit.Add("set timestamp");
            command_unit.Add("splot '" + vertex_filename + "' with lines lt 3, \\");
            command_unit.Add("      '" + node_filename + "' with points pt 7 lt 0");

            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file '" + command_filename + "'");
        }

        public static void wedge_vertices(ref double[] v1, ref double[] v2, ref double[] v3, ref double[] v4,
        ref double[] v5, ref double[] v6 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEDGE_VERTICES returns the vertices of the unit wege.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double V1[3], V2[3], V3[3], V4[3], V5[3], V6[3],
        //    the vertices.
        //
        {
            double[] v1_save = {
                0.0, 0.0, -1.0
            }
            ;
            double[] v2_save = {
                1.0, 0.0, -1.0
            }
            ;
            double[] v3_save = {
                0.0, 1.0, -1.0
            }
            ;
            double[] v4_save = {
                0.0, 0.0, +1.0
            }
            ;
            double[] v5_save = {
                1.0, 0.0, +1.0
            }
            ;
            double[] v6_save = {
                0.0, 1.0, +1.0
            }
            ;

            typeMethods.r8vec_copy(3, v1_save, ref v1);
            typeMethods.r8vec_copy(3, v2_save, ref v2);
            typeMethods.r8vec_copy(3, v3_save, ref v3);
            typeMethods.r8vec_copy(3, v4_save, ref v4);
            typeMethods.r8vec_copy(3, v5_save, ref v5);
            typeMethods.r8vec_copy(3, v6_save, ref v6);
        }
    }
}