using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.FDM;
using Burkardt.Types;

namespace FD1DAdvectionFTCSTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_ADVECTION_FTCS solves the advection equation using the FTCS method.
        //
        //  Discussion:
        //
        //    The FTCS method is unstable for the advection problem.
        //
        //    Given a smooth initial condition, successive FTCS approximations will
        //    exhibit erroneous oscillations of increasing magnitude.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 December 2012
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
        Console.WriteLine("FD1D_ADVECTION_FTCS:");
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
        Console.WriteLine("  We use a method known as FTCS:");
        Console.WriteLine("   FT: Forward Time  : du/dt = (u(t+dt,x)-u(t,x))/dt");
        Console.WriteLine("   CS: Centered Space: du/dx = (u(t,x+dx)-u(t,x-dx))/2/dx");

        const int nx = 101;
        const double dx = 1.0 / (nx - 1);
        const double a = 0.0;
        const double b = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(nx, a, b);
        const int nt = 1000;
        const double dt = 1.0 / nt;
        const double c = 1.0;

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

        double[] unew = new double[nx];

        for (i = 0; i < nt; i++)
        {
            for (j = 0; j < nx; j++)
            {
                int jm1 = typeMethods.i4_wrap(j - 1, 0, nx - 1);
                int jp1 = typeMethods.i4_wrap(j + 1, 0, nx - 1);
                unew[j] = u[j] - c * dt / dx / 2.0 * (u[jp1] - u[jm1]);
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
        command_unit.Add("set output 'advection.png'");
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
        Console.WriteLine("FD1D_ADVECTION_FTCS");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
        

}