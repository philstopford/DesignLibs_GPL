using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace FD1DAdvectionFTCSTest
{
    class Program
    {
        static void Main(string[] args)
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
            double a;
            double b;
            double c;
            string command_filename = "advection_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "advection_data.txt";
            List<string> data_unit = new List<string>();
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

            nx = 101;
            dx = 1.0 / (double) (nx - 1);
            a = 0.0;
            b = 1.0;
            x = typeMethods.r8vec_linspace_new(nx, a, b);
            nt = 1000;
            dt = 1.0 / (double) (nt);
            c = 1.0;

            u = initial_condition(nx, x);
            //
            //  Open data file, and write solutions as they are computed.
            //

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

            unew = new double[nx];

            for (i = 0; i < nt; i++)
            {
                for (j = 0; j < nx; j++)
                {
                    jm1 = typeMethods.i4_wrap(j - 1, 0, nx - 1);
                    jp1 = typeMethods.i4_wrap(j + 1, 0, nx - 1);
                    unew[j] = u[j] - c * dt / dx / 2.0 * (u[jp1] - u[jm1]);
                }

                for (j = 0; j < nx; j++)
                {
                    u[j] = unew[j];
                }

                if (i == nt_step - 1)
                {
                    t = (double) (i) * dt;
                    for (j = 0; j < nx; j++)
                    {
                        data_unit.Add("  " + x[j]
                            + "  " + t
                            + "  " + u[j] + "");
                    }

                    data_unit.Add("");
                    nt_step = nt_step + 100;
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
        
        static double[] initial_condition ( int nx, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INITIAL_CONDITION sets the initial condition.
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
            //  Parameters:
            //
            //    Input, int NX, the number of nodes.
            //
            //    Input, double X[NX], the coordinates of the nodes.
            //
            //    Output, double INITIAL_CONDITION[NX], the value of the initial condition.
            //
        {
            int i;
            double[] u;

            u = new double[nx];

            for ( i = 0; i < nx; i++ )
            {
                if  ( 0.4 <= x[i] && x[i] <= 0.6 )
                {
                    u[i] = Math.Pow ( 10.0 * x[i] - 4.0, 2 )
                           * Math.Pow ( 6.0 - 10.0 * x[i], 2 );
                }
                else
                {
                    u[i] = 0.0;
                }
            }
            return u;
        }
    }
}