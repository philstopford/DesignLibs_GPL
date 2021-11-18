using System;
using System.Collections.Generic;
using System.IO;

namespace Spring2ODE;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPRING_ODE2.
        //
        //  Discussion:
        //
        //    This is a revision of the SPRING_ODE code.
        //
        //    In this revision of the program, we want to use vectors (C arrays) to 
        //    store the data, and we want to write the data out to a file in a form 
        //    that Gnuplot (or other plotting programs) can use.
        //
        //    Hooke's law for a spring observes that the restoring force is
        //    proportional to the displacement: F = - k x
        //
        //    Newton's law relates the force to acceleration: F = m a
        //
        //    Putting these together, we have
        //
        //      m * d^2 x/dt^2 = - k * x
        //
        //    We can add a damping force with coefficient c:
        //
        //      m * d^2 x/dt^2 = - k * x - c * dx/dt
        //
        //    If we write this as a pair of first order equations for (x,v), we have
        //
        //          dx/dt = v
        //      m * dv/dt = - k * x - c * v
        //
        //    and now we can approximate these values for small time steps.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    None
        //
    {
        double c;
        string command_filename = "spring_ode2_commands.txt";
        List<string> command_unit = new();
        string data_filename = "spring_ode2_data.txt";
        List<string> data_unit = new();
        double dt;
        int i;
        double k;
        double m;
        int n = 101;
        double[] t = new double[101];
        double t_final;
        double[] v = new double[101];
        double[] x = new double[101];

        Console.WriteLine("");
        Console.WriteLine("SPRING_ODE2");
        Console.WriteLine("  Approximate the solution of a spring equation.");
        Console.WriteLine("  Write data to a file for use by gnuplot.");
        //
        //  Data
        //
        m = 1.0;
        k = 1.0;
        c = 0.3;
        t_final = 20.0f;
        dt = t_final / (n - 1);
        //
        //  Store the initial conditions in entry 0.
        //
        t[0] = 0.0;
        x[0] = 1.0;
        v[0] = 0.0;
        //
        //  Compute the approximate solution at equally spaced times 
        //  in entries 1 through N-1.
        //
        for (i = 1; i < n; i++)
        {
            t[i] = i * t_final / (n - 1);
            x[i] = x[i - 1] + dt * v[i - 1];
            v[i] = v[i - 1] + dt / m * (-k * x[i - 1] - c * v[i - 1]);
        }

        //
        //  Create the plot data file.
        //
        for (i = 0; i < n; i++)
        {
            data_unit.Add("  " + t[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + v[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("  Created data file \"" + data_filename + "\".");
        //
        //  Create the plot command file.
        //
        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output 'xv_time.png'");
        command_unit.Add("set xlabel '<--- T --->'");
        command_unit.Add("set ylabel '<--- X(T), V(T) --->'");
        command_unit.Add("set title 'Position and Velocity versus Time'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename
                                  + "' using 1:2 lw 3 linecolor rgb 'blue',"
                                  + " '' using 1:3 lw 3 linecolor rgb 'red'");
        command_unit.Add("set output 'xv_phase.png'");
        command_unit.Add("set xlabel '<--- X(T) --->'");
        command_unit.Add("set ylabel '<--- V(T) --->'");
        command_unit.Add("set title 'Position versus Velocity'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename
                                  + "' using 2:3 lw 3 linecolor rgb 'green'");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);
        Console.WriteLine("  Created command file '" + command_filename + "'");

        Console.WriteLine("");
        Console.WriteLine("SPRING_ODE2:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}