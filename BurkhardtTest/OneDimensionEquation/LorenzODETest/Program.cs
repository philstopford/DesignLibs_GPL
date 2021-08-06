using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.Types;

namespace LorenzODETest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LORENZ_ODE.
            //
            //  Discussion:
            //
            //    Thanks to Ben Whitney for pointing out an error in specifying a loop,
            //    24 May 2016.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string command_filename = "lorenz_ode_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "lorenz_ode_data.txt";
            List<string> data_unit = new List<string>();
            double dt;
            int i;
            int j;
            int m = 3;
            int n = 200000;
            double[] t;
            double t_final;
            double[] x;
            double[] xnew;

            Console.WriteLine("");
            Console.WriteLine("LORENZ_ODE");
            Console.WriteLine("  Compute solutions of the Lorenz system.");
            Console.WriteLine("  Write data to a file for use by gnuplot.");
            //
            //  Data
            //
            t_final = 40.0;
            dt = t_final / (double) (n);
            //
            //  Store the initial conditions in entry 0.
            //
            t = typeMethods.r8vec_linspace_new(n + 1, 0.0, t_final);
            x = new double[m * (n + 1)];
            x[0 + 0 * m] = 8.0;
            x[0 + 1 * m] = 1.0;
            x[0 + 2 * m] = 1.0;
            //
            //  Compute the approximate solution at equally spaced times.
            //
            for (j = 0; j < n; j++)
            {
                xnew = RungeKutta.rk4vec(t[j], m, x, dt, Lorenz.lorenz_rhs, index: +j * m);
                for (i = 0; i < m; i++)
                {
                    x[i + (j + 1) * m] = xnew[i];
                }
            }

            //
            //  Create the plot data file.
            //
            for (j = 0; j <= n; j = j + 50)
            {
                data_unit.Add("  " + t[j].ToString().PadLeft(14)
                                   + "  " + x[0 + j * m].ToString().PadLeft(14)
                                   + "  " + x[1 + j * m].ToString().PadLeft(14)
                                   + "  " + x[2 + j * m].ToString().PadLeft(14) + "");
            }

            File.WriteAllLines(data_filename, data_unit);

            Console.WriteLine("  Created data file \"" + data_filename + "\".");
            /*
            Create the plot command file.
            */
            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'xyz_time.png'");
            command_unit.Add("set xlabel '<--- T --->'");
            command_unit.Add("set ylabel '<--- X(T), Y(T), Z(T) --->'");
            command_unit.Add("set title 'X, Y and Z versus Time'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename
                                      + "' using 1:2 lw 3 linecolor rgb 'blue',"
                                      + "'' using 1:3 lw 3 linecolor rgb 'red',"
                                      + "'' using 1:4 lw 3 linecolor rgb 'green'");
            command_unit.Add("set output 'xyz_3d.png'");
            command_unit.Add("set xlabel '<--- X(T) --->'");
            command_unit.Add("set ylabel '<--- Y(T) --->'");
            command_unit.Add("set zlabel '<--- Z(T) --->'");
            command_unit.Add("set title '(X(T),Y(T),Z(T)) trajectory'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("splot '" + data_filename
                                       + "' using 2:3:4 lw 1 linecolor rgb 'blue'");
            command_unit.Add("quit");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created command file '" + command_filename + "'");

            Console.WriteLine("");
            Console.WriteLine("LORENZ_ODE:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}