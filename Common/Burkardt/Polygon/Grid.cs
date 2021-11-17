using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.Polygon;

public static class Grid
{
    public static int polygon_grid_count(int n, int nv)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_GRID_COUNT counts the grid points inside a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals on a side.
        //
        //    Input, int NV, the number of vertices.
        //    3 <= NV.
        //
        //    Output, int POLYGON_GRID_COUNT, the number of grid points.
        //
    {
        int ng;

        ng = 1 + nv * n * (n + 1) / 2;

        return ng;
    }

    public static void polygon_grid_display(int n, int nv, double[] v, int ng, double[] xg,
            string prefix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_GRID_DISPLAY displays grid points inside a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, int NV, the number of vertices in the polygon.
        //
        //    Input, double V[2*NV], the coordinates of the vertices.
        //
        //    Input, int NG, the number of grid points.
        //
        //    Input, double XG[2*NG], the grid points.
        //
        //    Input, string PREFIX, a string used to name the files.
        //
    {
        List<string> command_unit = new();
        string command_filename;
        List<string> grid_unit = new();
        string grid_filename;
        int j;
        string plot_filename;
        double[] vc = new double[2];
        List<string> vertex_unit = new();
        string vertex_filename;
        //
        //  Determine the centroid.
        //
        vc[0] = 0.0;
        vc[1] = 0.0;
        for (j = 0; j < nv; j++)
        {
            vc[0] += v[0 + j * 2];
            vc[1] += v[1 + j * 2];
        }

        vc[0] /= nv;
        vc[1] /= nv;
        //
        //  Write the vertex file.
        //
        vertex_filename = prefix + "_vertex.txt";

        for (j = 0; j < nv; j++)
        {
            vertex_unit.Add("  " + v[0 + 2 * j]
                                 + "  " + v[1 + 2 * j] + "");
        }

        vertex_unit.Add("  " + v[0 + 0 * j]
                             + "  " + v[1 + 0 * j] + "");
        for (j = 0; j < nv; j++)
        {
            vertex_unit.Add("");
            vertex_unit.Add("  " + v[0 + j * 2]
                                 + "  " + v[1 + j * 2] + "");
            vertex_unit.Add("  " + vc[0]
                                 + "  " + vc[1] + "");
        }

        File.WriteAllLines(vertex_filename, vertex_unit);
        Console.WriteLine("");
        Console.WriteLine("  Created vertex file '" + vertex_filename + "'");
        //
        //  Write the gridpoint file.
        //
        grid_filename = prefix + "_grid.txt";
        for (j = 0; j < ng; j++)
        {
            grid_unit.Add("  " + xg[0 + j * 2]
                               + "  " + xg[1 + j * 2] + "");
        }

        File.WriteAllLines(grid_filename, grid_unit);
        Console.WriteLine("");
        Console.WriteLine("  Created grid file '" + grid_filename + "'");
        //
        //  Write the command file.
        //
        plot_filename = prefix + ".png";

        command_filename = prefix + "_commands.txt";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set xlabel '<--- X --->'");
        command_unit.Add("set ylabel '<--- Y --->'");
        command_unit.Add("set title '" + prefix + "'");
        command_unit.Add("set grid");
        command_unit.Add("set key off");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set style data lines");

        command_unit.Add(
            "plot '" + grid_filename + "' using 1:2 with points lt 3 pt 3,\\");
        command_unit.Add(
            "     '" + vertex_filename + "' using 1:2 lw 3 linecolor rgb 'black'");
        command_unit.Add("quit");
        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Created command file '" + command_filename + "'");
    }

    public static double[] polygon_grid_points(int n, int nv, double[] v, int ng)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_GRID_POINTS computes points on a polygonal grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, int NV, the number of vertices in the polygon.
        //
        //    Input, double V[2*NV], the coordinates of the vertices.
        //
        //    Input, int NG, the number of grid points.
        //
        //    Output, double POLYGON_GRID_POINTS[2*NG], the coordinates of the 
        //    grid points.
        //
    {
        int i;
        int j;
        int k;
        int l;
        int lp1;
        int p;
        double[] vc = new double[2];
        double[] xg;

        xg = new double[2 * ng];
        p = 0;
        //
        //  Determine the centroid.
        //
        vc[0] = 0.0;
        vc[1] = 0.0;
        for (j = 0; j < nv; j++)
        {
            vc[0] += v[0 + j * 2];
            vc[1] += v[1 + j * 2];
        }

        vc[0] /= nv;
        vc[1] /= nv;
        //
        //  The centroid is the first point.
        //
        xg[0 + p * 2] = vc[0];
        xg[1 + p * 2] = vc[1];
        p += 1;
        //
        //  Consider each triangle formed by two consecutive vertices and the centroid,
        //  but skip the first line of points.
        //
        for (l = 0; l < nv; l++)
        {
            lp1 = (l + 1) % nv;
            for (i = 1; i <= n; i++)
            {
                for (j = 0; j <= n - i; j++)
                {
                    k = n - i - j;
                    xg[0 + p * 2] = (i * v[0 + l * 2]
                                     + j * v[0 + lp1 * 2]
                                     + k * vc[0])
                                    / n;
                    xg[1 + p * 2] = (i * v[1 + l * 2]
                                     + j * v[1 + lp1 * 2]
                                     + k * vc[1])
                                    / n;
                    p += 1;
                }
            }
        }

        return xg;
    }
}