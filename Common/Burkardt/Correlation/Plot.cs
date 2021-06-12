using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.CorrelationNS
{
    public static class Plot
    {
        public static void correlation_plot(int n, double[] rho, double[] c, string header,
                string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_PLOT makes a plot of a correlation function.
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
            //    Input, int N, the number of arguments.
            //
            //    Input, double RHO[N], the arguments.
            //
            //    Input, double C[N], the correlations.
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

            data_filename = header + "_data.txt";

            for (i = 0; i < n; i++)
            {
                data_unit.Add("  " + rho[i] + "  " + c[i] + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file \"" + data_filename + "\".");
            
            command_filename = header + "_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output \"" + header + "_plot.png\"");
            command_unit.Add("set xlabel 'Distance Rho'");
            command_unit.Add("set ylabel 'Correlation C(Rho)'");
            command_unit.Add("set title '" + title + "'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3 linecolor rgb 'blue'");
            command_unit.Add("quit");

            File.WriteAllLines(command_filename, command_unit);
            
            Console.WriteLine("  Created command file \"" + command_filename + "\".");
        }

        public static void correlation_plots(int n, int n2, double[] rho, double[] rho0, double[] c,
                string header, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CORRELATION_PLOTS plots correlations for a range of correlation lengths.
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
            //    Input, int N, the number of values of RHO.
            //
            //    Input, int N2, the number of values of RHO0.
            //
            //    Input, double RHO[N], the independent value.
            //
            //    Input, double RHO0[N2], the correlation lengths.
            //
            //    Input, double C[N*N2], the correlations.
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

            data_filename = header + "_plots_data.txt";

            for (i = 0; i < n; i++)
            {
                string line = "  " + rho[i];
                for (j = 0; j < n2; j++)
                {
                    line += "  " + c[i + j * n];
                }

                data_unit.Add(line);
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("  Created data file \"" + data_filename + "\".");

            command_filename = header + "_plots_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output \"" + header + "_plots.png\"");
            command_unit.Add("set xlabel 'Rho'");
            command_unit.Add("set ylabel 'Correlation(Rho)'");
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
                command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3, \\");
                for (i = 2; i < n2; i++)
                {
                    command_unit.Add("     '" + data_filename + "' using 1:" + i + 1 + " lw 3, \\");
                }

                command_unit.Add("     '" + data_filename + "' using 1:" + n2 + 1 + " lw 3");
            }

            command_unit.Add("quit");
            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created command file \"" + command_filename + "\".");

            return;
        }
    }
}