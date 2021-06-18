using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.Types;

namespace FD1DAdvectionDiffusionSteadyTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FD1D_ADVECTION_DIFFUSION_STEADY solves steady advection diffusion equation.
            //
            //  Discussion:
            //
            //    The steady advection diffusion equation has the form:
            //
            //      v ux - k * uxx = 0
            //
            //    where V (the advection velocity) and K (the diffusivity) are positive 
            //    constants, posed in the region
            //
            //      a = 0 < x < 1 = b
            //
            //    with boundary conditions
            //
            //      u(0) = 0, u(1) = 1.
            //
            //    The discrete solution is unreliable when dx > 2 * k / v / ( b - a ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double[] a3;
            double b;
            string command_filename = "fd1d_advection_diffusion_steady_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "fd1d_advection_diffusion_steady_data.txt";
            List<string> data_unit = new List<string>();
            double dx;
            double[] f;
            int i;
            int j;
            double k;
            int nx;
            double r;
            double[] u;
            double v;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("FD1D_ADVECTION_DIFFUSION_STEADY:");
            Console.WriteLine("");
            Console.WriteLine("  Solve the 1D steady advection diffusion equation:,");
            Console.WriteLine("    v du/dx - k d2u/dx2 = 0");
            Console.WriteLine("  with constant, positive velocity V and diffusivity K");
            Console.WriteLine("  over the interval:");
            Console.WriteLine("    0.0 <= x <= 1.0");
            Console.WriteLine("  with boundary conditions:");
            Console.WriteLine("    u(0) = 0, u(1) = 1.");
            Console.WriteLine("");
            Console.WriteLine("  Use finite differences");
            Console.WriteLine("   d u/dx  = (u(t,x+dx)-u(t,x-dx))/2/dx");
            Console.WriteLine("   d2u/dx2 = (u(x+dx)-2u(x)+u(x-dx))/dx^2");
            //
            //  Physical constants.
            //
            v = 1.0;
            k = 0.05;
            Console.WriteLine("");
            Console.WriteLine("  Diffusivity K = " + k + "");
            Console.WriteLine("  Velocity V    = " + v + "");
            //
            //  Spatial discretization.
            //
            nx = 101;
            a = 0.0;
            b = 1.0;
            dx = (b - a) / (double) (nx - 1);
            x = typeMethods.r8vec_linspace_new(nx, a, b);

            Console.WriteLine("  Number of nodes NX = " + nx + "");
            Console.WriteLine("  DX = " + dx + "");
            Console.WriteLine("  Maximum safe DX is " + 2.0 * k / v / (b - a) + "");
            //
            //  Set up the tridiagonal linear system corresponding to the boundary 
            //  conditions and advection-diffusion equation.
            //
            a3 = new double[nx * 3];
            f = new double[nx];

            a3[0 + 1 * nx] = 1.0;
            f[0] = 0.0;

            for (i = 1; i < nx - 1; i++)
            {
                a3[i + 0 * nx] = -v / dx / 2.0 - k / dx / dx;
                a3[i + 1 * nx] = +2.0 * k / dx / dx;
                a3[i + 2 * nx] = +v / dx / 2.0 - k / dx / dx;
                f[i] = 0.0;
            }

            a3[nx - 1 + 1 * nx] = 1.0;
            f[nx - 1] = 1.0;

            u = Trisolve.trisolve(nx, a3, f);
            //
            //  The exact solution to the differential equation is known.
            //
            r = v * (b - a) / k;

            w = new double[nx];

            for (i = 0; i < nx; i++)
            {
                w[i] = (1.0 - Math.Exp(r * x[i])) / (1.0 - Math.Exp(r));
            }

            //
            //  Write data file.
            //
            for (j = 0; j < nx; j++)
            {
                data_unit.Add(x[j] + "  "
                    + u[j] + "  "
                    + w[j] + "");
            }

            File.WriteAllLines(data_filename, data_unit);

            Console.WriteLine("");
            Console.WriteLine("  Gnuplot data written to file '" + data_filename + "'.");
            //
            //  Write command file.
            //

            command_unit.Add("set term png");
            command_unit.Add("set output 'fd1d_advection_diffusion_steady.png'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("unset key");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---U(X)--->'");
            command_unit.Add("set title 'Exact: green line, Approx: red dots'");
            command_unit.Add("plot '" + data_filename
                + "' using 1:2 with points pt 7 ps 2,\\");
            command_unit.Add("'' using 1:3 with lines lw 3");
            command_unit.Add("quit");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Gnuplot commands written to '" + command_filename + "'");
            Console.WriteLine("");
            Console.WriteLine("FD1D_ADVECTION_DIFFUSION_STEADY");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}