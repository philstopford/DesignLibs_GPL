﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.CorrelationNS;

public static partial class Correlation
{
    public static double[] correlation_brownian(int m, int n, double[] s, double[] t, double rho0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_BROWNIAN computes the Brownian correlation function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of arguments.
        //
        //    Input, double S[M], T[N], two samples.
        //    0 <= S(*), T(*).
        //
        //    Input, double RHO0, the correlation length.
        //
        //    Output, double C[M*N], the correlations.
        //
    {
        int j;

        double[] c = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                c[i + j * m] = Math.Max(s[i], t[j]) switch
                {
                    > 0.0 => Math.Sqrt(Math.Min(s[i], t[j]) / Math.Max(s[i], t[j])),
                    _ => 1.0
                };
            }
        }

        return c;
    }

    public static void correlation_brownian_display()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CORRELATION_BROWNIAN_DISPLAY displays 4 slices of the Brownian Correlation.
        //
        //  Discussion:
        //
        //    The correlation function is C(S,T) = sqrt ( min ( s, t ) / max ( s, t ) ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string command_filename = "brownian_plots_commands.txt";
        List<string> command_unit = new();
        const string data_filename = "brownian_plots_data.txt";
        List<string> data_unit = new();
        const int n = 101;
        const int n2 = 4;
        double[] t = {0.25, 1.50, 2.50, 3.75};

        double[] s = typeMethods.r8vec_linspace_new(n, 0.0, 5.0);

        double[] c = new double[n * n2];

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n2; j++)
            {
                c[i + j * n] = Math.Sqrt(Math.Min(s[i], t[j]) / Math.Max(s[i], t[j]));
            }
        }

        for (int i = 0; i < n; i++)
        {
            string line = "  " + s[i];
            for (int j = 0; j < n2; j++)
            {
                line += "  " + c[i + j * n];
            }

            data_unit.Add(line);
        }

        File.WriteAllLines(data_filename, data_unit);
        Console.WriteLine("  Created data file \"" + data_filename + "\".");

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set key off");
        command_unit.Add("set output \"brownian_plots.png\"");
        command_unit.Add("set title 'Brownian correlation C(S,T), S = 0.25, 1.5, 2.5, 3.75'");
        command_unit.Add("set xlabel 'S'");
        command_unit.Add("set ylabel 'C(s,t)'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot \"" + data_filename + "\" using 1:2 lw 3 linecolor rgb 'blue',\\");
        command_unit.Add("     \"" + data_filename + "\" using 1:3 lw 3 linecolor rgb 'blue',\\");
        command_unit.Add("     \"" + data_filename + "\" using 1:4 lw 3 linecolor rgb 'blue',\\");
        command_unit.Add("     \"" + data_filename + "\" using 1:5 lw 3 linecolor rgb 'blue'");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);
        Console.WriteLine("  Created command file \"" + command_filename + "\".");

    }
}