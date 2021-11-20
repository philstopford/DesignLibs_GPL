using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.NavierStokesNS;

public static class NavierStokes2DExact
{
    public static void grid_2d(int x_num, double x_lo, double x_hi, int y_num, double y_lo,
            double y_hi, ref double[] x, ref double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRID_2D returns a regular 2D grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int X_NUM, the number of X values to use.
        //
        //    double X_LO, X_HI, the range of X values.
        //
        //    int Y_NUM, the number of Y values to use.
        //
        //    double Y_LO, Y_HI, the range of Y values.
        //
        //  Output:
        //
        //    double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM], 
        //    the coordinates of the grid.
        //
    {
        int i;
        int j;

        switch (x_num)
        {
            case 1:
            {
                for (j = 0; j < y_num; j++)
                {
                    for (i = 0; i < x_num; i++)
                    {
                        x[i + j * x_num] = (x_lo + x_hi) / 2.0;
                    }
                }

                break;
            }
            default:
            {
                for (i = 0; i < x_num; i++)
                {
                    double xi = ((x_num - i - 1) * x_lo
                                 + i * x_hi)
                                / (x_num - 1);
                    for (j = 0; j < y_num; j++)
                    {
                        x[i + j * x_num] = xi;
                    }
                }

                break;
            }
        }

        switch (y_num)
        {
            case 1:
            {
                for (j = 0; j < y_num; j++)
                {
                    for (i = 0; i < x_num; i++)
                    {
                        y[i + j * x_num] = (y_lo + y_hi) / 2.0;
                    }
                }

                break;
            }
            default:
            {
                for (j = 0; j < y_num; j++)
                {
                    double yj = ((y_num - j - 1) * y_lo
                                 + j * y_hi)
                                / (y_num - 1);
                    for (i = 0; i < x_num; i++)
                    {
                        y[i + j * x_num] = yj;
                    }
                }

                break;
            }
        }

    }

    public static void ns2de_gnuplot(string header, int n, double[] x, double[] y, double[] u,
            double[] v, double[] p, double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NS2DE_GNUPLOT writes the Navier-Stokes solution to files for GNUPLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    string HEADER, a header to be used to name the files.
        //
        //    int N, the number of evaluation points.
        //
        //    double X[N], Y[N], the coordinates of the evaluation points.
        //
        //    double U[N], V[N], P[N], the solution samples.
        //
        //    double S, a scale factor for the velocity vectors.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        //
        //  Write the data file.
        //
        string data_filename = header + "_data.txt";


        for (i = 0; i < n; i++)
        {
            data_unit.Add("  " + x[i]
                               + "  " + y[i]
                               + "  " + u[i]
                               + "  " + v[i]
                               + "  " + s * u[i]
                               + "  " + s * v[i]
                               + "  " + p[i] + "");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Data written to '" + data_filename + "'");
        //
        //  Write the command file.
        //
        string command_filename = header + "_commands.txt";
        string plot_filename = header + ".png";

        command_unit.Add("#  " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("#");
        command_unit.Add("#  Add titles and labels.");
        command_unit.Add("#");
        command_unit.Add("set xlabel '<--- X --->'");
        command_unit.Add("set ylabel '<--- Y --->'");
        command_unit.Add("set title 'Navier-Stokes velocity field'");
        command_unit.Add("unset key");
        command_unit.Add("#");
        command_unit.Add("#  Add grid lines.");
        command_unit.Add("#");
        command_unit.Add("set grid");
        command_unit.Add("set size ratio -1");
        command_unit.Add("#");
        command_unit.Add("#  Timestamp the plot.");
        command_unit.Add("#");
        command_unit.Add("set timestamp");
        command_unit.Add("plot '" + data_filename
                                  + "' using 1:2:5:6 with vectors \\");
        command_unit.Add("  head filled lt 2 linecolor rgb 'blue'");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Commands written to '" + command_filename + "'");
    }
}