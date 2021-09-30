using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.PathsNS
{
    public static class Paths
    {
        public static void paths_plot(int n, int n2, double[] rho, double[] x, string header,
        string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PATHS_PLOT plots a sequence of paths or simulations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points in each path.
        //
        //    Input, int N2, the number of paths.
        //
        //    Input, double RHO[N], the independent value.
        //
        //    Input, double X[N*N2], the path values.
        //
        //    Input, string HEADER, an identifier for the files.
        //
        //    Input, string TITLE, a title for the plot.
        //
        {
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            int i;
            int j;
            //double rho0;

            data_filename = header + "_path_data.txt";

            for (i = 0; i < n; i++)
            {
                string line = "  " + rho[i];
                for (j = 0; j < n2; j++)
                {
                    line += "  " + x[i + j * n];
                }

                data_unit.Add(line);;
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file \"" + data_filename + "\".");

            command_filename = header + "_path_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output \"" + header + "_paths.png\"");
            command_unit.Add("set xlabel 'Rho'");
            command_unit.Add("set ylabel 'X(Rho)'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("set key off");
            if (n2 == 1)
            {
                command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3");
            }
            else
            {
                command_unit.Add("plot '" + data_filename + "' using 1:2, \\");
                for (i = 2; i < n2; i++)
                {
                    command_unit.Add("     '" + data_filename + "' using 1:" + (i + 1) + ", \\");
                }

                command_unit.Add("     '" + data_filename + "' using 1:" + (n2 + 1) + "");
            }

            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file \"" + command_filename + "\".");
        }
    }
}