using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.FDM;
using Burkardt.Types;

namespace FD1DAdvectionLAXWendroffTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_ADVECTION_LAX_WENDROFF: advection equation using Lax-Wendroff method.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 March 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string command_filename = "advection_commands.txt";
        List<string> command_unit = new();
        const string data_filename = "advection_data.txt";
        List<string> data_unit = new();
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("FD1D_ADVECTION_LAX_WENDROFF:");
        Console.WriteLine("");
        Console.WriteLine("  Solve the constant-velocity advection equation in 1D,");
        Console.WriteLine("    du/dt = - c du/dx");
        Console.WriteLine("  over the interval:");
        Console.WriteLine("    0.0 <= x <= 1.0");
        Console.WriteLine("  with periodic boundary conditions, and");
        Console.WriteLine("  with a given initial condition");
        Console.WriteLine("    u(0,x) = (10x-4)^2 (6-10x)^2 for 0.4 <= x <= 0.6");
        Console.WriteLine("           = 0 elsewhere.");
        Console.WriteLine("");
        Console.WriteLine("  We modify the FTCS method using the Lax-Wendroff method:");

        const int nx = 101;
        const double dx = 1.0 / (nx - 1);
        const double a = 0.0;
        const double b = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(nx, a, b);
        const int nt = 1000;
        const double dt = 1.0 / nt;
        const double c = 1.0;
        const double c1 = 0.5 * c * dt / dx;
        double c2 = 0.5 * Math.Pow(c * dt / dx, 2);

        double[] u = InitialCondition.initial_condition(nx, x);
        //
        //  Open data file, and write solutions as they are computed.
        //

        double t = 0.0;
        data_unit.Add("  " + x[0]
                           + "  " + t
                           + "  " + u[0] + "");
        for (j = 0; j < nx; j++)
        {
            data_unit.Add("  " + x[j]
                               + "  " + t
                               + "  " + u[j] + "");
        }

        data_unit.Add("");

        int nt_step = 100;

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes NX = " + nx + "");
        Console.WriteLine("  Number of time steps NT = " + nt + "");
        Console.WriteLine("  Constant velocity C = " + c + "");
        Console.WriteLine("  CFL condition: dt (" + dt + ") <= dx / c (" + dx / c + ")");

        double[] unew = new double[nx];

        for (i = 0; i < nt; i++)
        {
            for (j = 0; j < nx; j++)
            {
                int jm1 = typeMethods.i4_wrap(j - 1, 0, nx - 1);
                int jp1 = typeMethods.i4_wrap(j + 1, 0, nx - 1);
                unew[j] = u[j] - c1 * (u[jp1] - u[jm1]) + c2 * (u[jp1] - 2.0 * u[j] + u[jm1]);
            }

            for (j = 0; j < nx; j++)
            {
                u[j] = unew[j];
            }

            if (i != nt_step - 1)
            {
                continue;
            }

            t = i * dt;
            for (j = 0; j < nx; j++)
            {
                data_unit.Add("  " + x[j]
                                   + "  " + t
                                   + "  " + u[j] + "");
            }

            data_unit.Add("");
            nt_step += 100;
        }

        //
        //  Close the data file once the computation is done.
        //
        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Plot data written to the file \"" + data_filename + "\"");
        //
        //  Write gnuplot command file.
        //

        command_unit.Add("set term png");
        command_unit.Add("set output 'advection_lax_wendroff.png'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("unset key");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---Time--->'");
        command_unit.Add("splot '" + data_filename + "' using 1:2:3 with lines");
        command_unit.Add("quit");

        File.WriteAllLines(command_filename, command_unit);
            
        Console.WriteLine("  Gnuplot command data written to the file \"" + command_filename + "\"");

        Console.WriteLine("");
        Console.WriteLine("FD1D_ADVECTION_LAX_WENDROFF:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}