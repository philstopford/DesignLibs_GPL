﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Polygon;
using Burkardt.Table;
using Burkardt.Types;

namespace PolygonTriangulateTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POLYGON_TRIANGULATE_TEST.
        //
        //  Discussion:
        //
        //    POLYGON_TRIANGULATE_TEST tests the POLYGON_TRIANGULATE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POLYGON_TRIANGULATE_TEST");
        Console.WriteLine("  Test the POLYGON_TRIANGULATE library.");

        test01();

        test02("comb");
        test02("hand");
        test02("i18");

        test03("comb");
        test03("hand");
        test03("i18");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("POLYGON_TRIANGULATE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests the "comb_10" polygon.
        //
        //  Discussion:
        //
        //    There are N-3 triangles in the triangulation.
        //
        //    For the first N-2 triangles, the first edge listed is always an
        //    internal diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2014
        //
    {
        const int n = 10;
        double[] x =
        {
            8.0, 7.0, 6.0, 5.0, 4.0,
            3.0, 2.0, 1.0, 0.0, 4.0
        };
        double[] y =
        {
            0.0, 10.0, 0.0, 10.0, 0.0,
            10.0, 0.0, 10.0, 0.0, -2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Triangulate the comb_10 polygon.");

        int[] triangles = Triangulate.polygon_triangulate(n, x, y);

        typeMethods.i4mat_transpose_print(3, n - 2, triangles, "  Triangles");
    }

    private static void test02(string prefix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 triangulates a polygon described in a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2014
        //
        //  Author:
        //
        //    John Burkardt.
        //
    {
        int i;
        //
        //  Create filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Read polygon coordinates in \"" + node_filename + "\"");
        //
        //  Read the node coordinates.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);

        int n = h.n;

        double[] xy = typeMethods.r8mat_data_read(node_filename, 2, n);
        //
        //  Get the triangulation.
        //
        double[] x = new double[n];
        double[] y = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = xy[0 + i * 2];
            y[i] = xy[1 + i * 2];
        }

        int[] triangles = Triangulate.polygon_triangulate(n, x, y);
        //
        //  Write the triangulation to a file.
        //
        int triangle_num = n - 2;
        typeMethods.i4mat_write(element_filename, 3, triangle_num, triangles);

        Console.WriteLine("  Write triangulation to \"" + element_filename + "\"");
    }

    private static void test03(string prefix)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 plots a triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2014
        //
        //  Author:
        //
        //    John Burkardt.
        //
    {
        List<string> command_unit = new();
        List<string> diagonal_unit = new();
        List<string> edge_unit = new();
        int i;
        int j;
        //
        //  Create filenames.
        //
        string node_filename = prefix + "_nodes.txt";
        string edge_filename = prefix + "_edges.txt";
        string diagonal_filename = prefix + "_diagonals.txt";
        string command_filename = prefix + "_commands.txt";
        string plot_filename = prefix + ".png";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Read node coordinates in \"" + node_filename + "\"");
        //
        //  Read the node coordinates.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int n = h.n;

        double[] xy = typeMethods.r8mat_data_read(node_filename, 2, n);
        //
        //  Get the triangulation.
        //
        double[] x = new double[n];
        double[] y = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = xy[0 + i * 2];
            y[i] = xy[1 + i * 2];
        }

        int[] triangles = Triangulate.polygon_triangulate(n, x, y);
        //
        //  Plot the edges.
        //
        for (j = 0; j < n + 1; j++)
        {
            int j2 = j % n;
            edge_unit.Add(xy[0 + j2 * 2] + "  "
                                         + xy[1 + j2 * 2] + "");
        }

        File.WriteAllLines(edge_filename, edge_unit);

        //
        //  Plot the diagonals.
        //

        for (j = 0; j < n - 3; j++)
        {
            for (i = 0; i < 2; i++)
            {
                int node = triangles[i + j * 3];
                diagonal_unit.Add(xy[0 + node * 2] + "  "
                                                   + xy[1 + node * 2] + "");
            }

            diagonal_unit.Add("");
        }

        File.WriteAllLines(diagonal_filename, diagonal_unit);
        //
        //  Write the GNUPLOT command file.
        //

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output \"" + plot_filename + "\"");
        command_unit.Add("set nokey");
        command_unit.Add("set size ratio 1");
        command_unit.Add("set timestamp");
        command_unit.Add("set xlabel \"<---X--->\"");
        command_unit.Add("set ylabel \"<---Y--->");
        command_unit.Add("set title \"Edges (green) and Diagonals (red)\"");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot \"" + edge_filename
                                   + "\" using 1:2 lw 3 linecolor rgb \"green\",\\");
        command_unit.Add("     \"" + diagonal_filename
                                   + "\" using 1:2 lw 3 linecolor rgb \"red\",\\");
        command_unit.Add("     \"" + node_filename
                                   + "\" using 1:2 with points pt 7 ps 2 lc rgb \"black\"");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("");
        Console.WriteLine("  Write edges to \"" + edge_filename + "\"");
        Console.WriteLine("  Write diagonals to \"" + diagonal_filename + "\"");
        Console.WriteLine("  Write gnuplot commands to \"" + command_filename + "\"");
    }
}