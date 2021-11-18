using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.PlotNS;

public static class Plot
{
    public static void plot_file(int m, int n, int[] c1, string title, string plot_filename,
            string png_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLOT_FILE writes the current configuration to a GNUPLOT plot file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int C1[M*N], the current state of the system.
        //
        //    Input, string TITLE, a title for the plot.
        //
        //    Input, string PLOT_FILENAME, a name for the GNUPLOT
        //    command file to be created.
        //
        //    Input, string PNG_FILENAME, the name of the PNG graphics
        //    file to be created.
        //
    {
        int j;
        List<string> plot_unit = new();

        double ratio = n / (double) m;

        plot_unit.Add("set term png");
        plot_unit.Add("set output \"" + png_filename + "\"");
        plot_unit.Add("set xrange [ 0 : " + m + " ]");
        plot_unit.Add("set yrange [ 0 : " + n + " ]");
        plot_unit.Add("set nokey");
        plot_unit.Add("set title \"" + title + "\"");
        plot_unit.Add("unset tics");

        plot_unit.Add("set size ratio " + ratio + "");
        for (j = 0; j < n; j++)
        {
            int y1 = j;
            int y2 = j + 1;
            int i;
            for (i = 0; i < m; i++)
            {
                int x1 = m - i - 1;
                int x2 = m - i;
                switch (c1[i + j * m])
                {
                    case < 0:
                        plot_unit.Add("set object rectangle from " + x1 + "," + y1 + " to "
                                      + x2 + "," + y2 + " fc rgb 'blue'");
                        break;
                    default:
                        plot_unit.Add("set object rectangle from " + x1 + "," + y1 + " to "
                                      + x2 + "," + y2 + " fc rgb 'red'");
                        break;
                }
            }
        }

        plot_unit.Add("plot 1");
        plot_unit.Add("quit");
        
        File.WriteAllLines(plot_filename, plot_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created the gnuplot graphics file \"" + plot_filename + "\"");
    }

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
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;

        string data_filename = header + "_data.txt";

        for (i = 0; i < n; i++)
        {
            data_unit.Add("  " + rho[i] + "  " + c[i] + "");
        }

        File.WriteAllLines(data_filename, data_unit);
        Console.WriteLine("  Created data file \"" + data_filename + "\".");

        string command_filename = header + "_commands.txt";

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
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;

        string data_filename = header + "_plots_data.txt";

        for (i = 0; i < n; i++)
        {
            string line = "  " + rho[i];
            int j;
            for (j = 0; j < n2; j++)
            {
                line += "  " + c[i + j * n];
            }

            data_unit.Add(line);
        }

        File.WriteAllLines(data_filename, data_unit);
        Console.WriteLine("  Created data file \"" + data_filename + "\".");

        string command_filename = header + "_plots_commands.txt";

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
        switch (n2)
        {
            case 1:
                command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3");
                break;
            default:
            {
                command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3, \\");
                for (i = 2; i < n2; i++)
                {
                    command_unit.Add("     '" + data_filename + "' using 1:" + (i + 1) + " lw 3, \\");
                }

                command_unit.Add("     '" + data_filename + "' using 1:" + (n2 + 1) + " lw 3");
                break;
            }
        }

        command_unit.Add("quit");
        File.WriteAllLines(command_filename, command_unit);
        Console.WriteLine("  Created command file \"" + command_filename + "\".");
    }

    public static void energy_plot(int it_num, double[] e_plot, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENERGY_PLOT plots the energy as a function of the iterations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int IT_NUM, the number of iterations to take.
        //
        //    Input, double E_PLOT[IT_NUM+1], the energy per iteration.
        //
        //    Input, string HEADER, an identifying string.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int it;
        //
        //  Write data file.
        //
        string data_filename = header + "_energy_data.txt";
        for (it = 0; it <= it_num; it++)
        {
            switch (e_plot[it])
            {
                case > 0.0:
                    data_unit.Add(it + "  "
                                     + Math.Log(e_plot[it]) + "");
                    break;
            }
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Gnuplot data written to file '" + data_filename + "'");
        //
        //  Write command file.
        //
        string command_filename = header + "_energy_commands.txt";

        string plot_filename = header + "_energy.png";

        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("set timestamp");
        command_unit.Add("unset key");
        command_unit.Add("set xlabel '<---Iteration--->'");
        command_unit.Add("set ylabel '<---Log(Energy)--->'");
        command_unit.Add("set title 'Energy Decrease with Iteration'");
        command_unit.Add("plot '" + data_filename + "' using 1:2 with points pt 7 ps 1");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Gnuplot commands written to '" + command_filename + "'");
    }

    public static void evolution_plot(int n, int it_num, double[] x_plot, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EVOLUTION_PLOT plots all points as a function of the iterations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, int IT_NUM, the number of iterations to take.
        //
        //    Input, double X_PLOT[N*IT_NUM], the point locations over time.
        //
        //    Input, string HEADER, an identifying string.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int it;
        //
        //  Write data file.
        //
        string data_filename = header + "_evolution_data.txt";

        for (it = 0; it <= it_num; it++)
        {
            string tmp = it + "  ";
            int i;
            for (i = 0; i < n; i++)
            {
                tmp += x_plot[i + it * n] + "  ";
            }

            data_unit.Add(tmp);
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Gnuplot data written to file '" + data_filename + "'");
        //
        //  Write command file.
        //
        string command_filename = header + "_evolution_commands.txt";

        string plot_filename = header + "_evolution.png";

        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("set timestamp");
        command_unit.Add("unset key");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---Iteration--->'");
        command_unit.Add("set title 'Point Motion with Iteration'");
        command_unit.Add("plot for [i=2:" + (n + 1) + "] '"
                         + data_filename + "' using i:1 with points pt 7 ps 1");

        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Gnuplot commands written to '" + command_filename + "'");
    }

    public static void motion_plot(int it_num, double[] xm_plot, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOTION_PLOT plots the motion as a function of the iterations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int IT_NUM, the number of iterations to take.
        //
        //    Input, double XM_PLOT[IT_NUM], the average motion per iteration.
        //
        //    Input, string HEADER, an identifying string.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int it;
        //
        //  Write data file.
        //
        string data_filename = header + "_motion_data.txt";

        for (it = 1; it <= it_num; it++)
        {
            switch (xm_plot[it - 1])
            {
                case > 0.0:
                    data_unit.Add(it + "  "
                                     + Math.Log(xm_plot[it - 1]) + "");
                    break;
            }
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Gnuplot data written to file '" + data_filename + "'");
        //
        //  Write command file.
        //
        string command_filename = header + "_motion_commands.txt";

        string plot_filename = header + "_motion.png";

        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("set timestamp");
        command_unit.Add("unset key");
        command_unit.Add("set xlabel '<---Iteration--->'");
        command_unit.Add("set ylabel '<---Average Motion--->'");
        command_unit.Add("set title 'Generator Motion with Iteration'");
        command_unit.Add("plot '" + data_filename
                                  + "' using 1:2 with points pt 7 ps 1");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Gnuplot commands written to '" + command_filename + "'");
    }
}