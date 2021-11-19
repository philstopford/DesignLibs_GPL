using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.FDM;
using Burkardt.Types;

namespace FD1DAdvectionLAXTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_ADVECTION_LAX solves the advection equation using the Lax method.
        //
        //  Discussion:
        //
        //    The Lax method is stable for the advection problem, if the time step
        //    satisifies the Courant-Friedrichs-Levy (CFL) condition:
        //
        //      dt <= dx / c
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double c;
        string command_filename = "advection_commands.txt";
        List<string> command_unit = new();
        string data_filename = "advection_data.txt";
        List<string> data_unit = new();
        double dt;
        double dx;
        int i;
        int j;
        int jm1;
        int jp1;
        int nx;
        int nt;
        int nt_step;
        double t;
        double[] u;
        double[] unew;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("FD1D_ADVECTION_LAX:");
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
        Console.WriteLine("  We modify the FTCS method using the Lax method:");
        Console.WriteLine("    du/dt = (u(t+dt,x)-0.5*u(t,x-dx)-0.5*u(t,x+dx))/dt");
        Console.WriteLine("    du/dx = (u(t,x+dx)-u(t,x-dx))/2/dx");

        nx = 101;
        dx = 1.0 / (nx - 1);
        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(nx, a, b);
        nt = 1000;
        dt = 1.0 / nt;
        c = 1.0;

        u = InitialCondition.initial_condition(nx, x);

        t = 0.0;
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

        nt_step = 100;

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes NX = " + nx + "");
        Console.WriteLine("  Number of time steps NT = " + nt + "");
        Console.WriteLine("  Constant velocity C = " + c + "");
        Console.WriteLine("  CFL condition: dt (" + dt + ") <= dx / c (" + dx / c + ")");

        unew = new double[nx];

        for (i = 0; i < nt; i++)
        {
            for (j = 0; j < nx; j++)
            {
                jm1 = typeMethods.i4_wrap(j - 1, 0, nx - 1);
                jp1 = typeMethods.i4_wrap(j + 1, 0, nx - 1);
                unew[j] = 0.5 * u[jp1] + 0.5 * u[jm1]
                          - c * dt / dx / 2.0 * (u[jp1] - u[jm1]);
            }

            for (j = 0; j < nx; j++)
            {
                u[j] = unew[j];
            }

            if (i == nt_step - 1)
            {
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
        command_unit.Add("set output 'advection_lax.png'");
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
        Console.WriteLine("FD1D_ADVECTION_LAX");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}