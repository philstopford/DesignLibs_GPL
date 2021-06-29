using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;
using Grid = Burkardt.Disk.Grid;

namespace DiskGridTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for DISK_GRID_TEST.
            //
            //  Discussion:
            //
            //    DISK_GRID_TEST tests the DISK_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("DISK_GRID_TEST:");
            Console.WriteLine("  Test the DISK_GRID library.");

            disk_grid_test01();
            disk_grid_test02();

            Console.WriteLine("");
            Console.WriteLine("DISK_GRID_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void disk_grid_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DISK_GRID_TEST01 tests DISK_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string boundary_filename = "disk_grid_test01_boundary.txt";
            List<string> boundary_unit = new List<string>();
            double[] c = new double[2];
            double[] cg;
            string command_filename = "disk_grid_test01_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "disk_grid_test01_data.txt";
            List<string> data_unit = new List<string>();
            string filename = "disk_grid_test01.xy";
            int i;
            int n;
            int ng;
            const double pi = 3.141592653589793;
            double r;
            double t;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  DISK_GRID can define a grid of points");
            Console.WriteLine("  with N+1 points on a horizontal or vertical radius,");
            Console.WriteLine("  based on any disk.");

            n = 20;
            r = 2.0;
            c[0] = 1.0;
            c[1] = 5.0;

            Console.WriteLine("");
            Console.WriteLine("  We use N = " + n + "");
            Console.WriteLine("  Radius R = " + r + "");
            ;
            Console.WriteLine("  Center C = (" + c[0] + "," + c[1] + ")");

            ng = Grid.disk_grid_count(n, r, c);

            Console.WriteLine("");
            Console.WriteLine("  Number of grid points will be " + ng + "");

            cg = Grid.disk_grid(n, r, c, ng);

            typeMethods.r82vec_print_part(ng, cg, 20, "  Part of the grid point array:");
            //
            //  Write the coordinate data to a file.
            //
            typeMethods.r8mat_write(filename, 2, ng, cg);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file \"" + filename + "\"");
            //
            //  Create graphics data files.
            //
            for (i = 0; i <= 50; i++)
            {
                t = 2.0 * pi * (double) (i) / 50.0;
                boundary_unit.Add("  " + c[0] + r * Math.Cos(t)
                    + "  " + c[1] + r * Math.Sin(t) + "");
            }
            
            File.WriteAllLines(boundary_filename, boundary_unit);

            Console.WriteLine("");
            Console.WriteLine("  Created boundary file \"" + boundary_filename + "\".");

            for (i = 0; i < ng; i++)
            {
                data_unit.Add("  " + cg[0 + i * 2]
                    + "  " + cg[1 + i * 2] + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created data file \"" + data_filename + "\"");
            //
            //  Create graphics command file.
            //
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'disk_grid_test01.png'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set title 'Disk Grid'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("set size ratio -1");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename
                + "' using 1:2 with points lt 3 pt 3,\\");
            command_unit.Add("    '" + boundary_filename
                + "' using 1:2 lw 3 linecolor rgb 'black'");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created command file \"" + command_filename + "\"");
        }

        static void disk_grid_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DISK_GRID_TEST02 tests DISK_GRID_FIBONACCI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string boundary_filename = "disk_grid_test02_boundary.txt";
            List<string> boundary_unit = new List<string>();
            double[] c = new double[2];
            string command_filename = "disk_grid_test02_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "disk_grid_test02_data.txt";
            List<string> data_unit = new List<string>();
            string filename = "disk_grid_test02.xy";
            double[] g;
            int i;
            int n;
            const double pi = 3.141592653589793;
            double r;
            double t;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  DISK_GRID_FIBONACCI can define a grid of N points");
            Console.WriteLine("  based on a Fibonacci spiral inside a disk.");

            n = 1000;
            r = 2.0;
            c[0] = 1.0;
            c[1] = 5.0;

            Console.WriteLine("");
            Console.WriteLine("  We use N = " + n + "");
            Console.WriteLine("  Radius R = " + r + "");
            ;
            Console.WriteLine("  Center C = (" + c[0] + "," + c[1] + ")");

            g = Grid.disk_grid_fibonacci(n, r, c);

            typeMethods.r82vec_print_part(n, g, 20, "  Part of the grid point array:");
            //
            //  Write the coordinate data to a file.
            //
            typeMethods.r8mat_write(filename, 2, n, g);

            Console.WriteLine("");
            Console.WriteLine("  Data written to the file \"" + filename + "\"");
            //
            //  Create graphics data files.
            //
            for (i = 0; i <= 50; i++)
            {
                t = 2.0 * pi * (double) (i) / 50.0;
                boundary_unit.Add("  " + c[0] + r * Math.Cos(t)
                    + "  " + c[1] + r * Math.Sin(t) + "");
            }

            File.WriteAllLines(boundary_filename, boundary_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created boundary file \"" + boundary_filename + "\".");

            for (i = 0; i < n; i++)
            {
                data_unit.Add("  " + g[0 + i * 2]
                    + "  " + g[1 + i * 2] + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created data file \"" + data_filename + "\"");
            //
            //  Create graphics command file.
            //
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'disk_grid_test02.png'");
            command_unit.Add("set xlabel '<--- X --->'");
            command_unit.Add("set ylabel '<--- Y --->'");
            command_unit.Add("set title 'Fibonacci Disk Grid'");
            command_unit.Add("set grid");
            command_unit.Add("set key off");
            command_unit.Add("set size ratio -1");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename
                + "' using 1:2 with points lt 3 pt 3,\\");
            command_unit.Add("    '" + boundary_filename
                + "' using 1:2 lw 3 linecolor rgb 'black'");
            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created command file \"" + command_filename + "\"");
            
            return;
        }
    }
}