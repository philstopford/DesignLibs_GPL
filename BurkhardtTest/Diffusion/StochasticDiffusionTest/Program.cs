using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Diffusion;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.StochasticDiffusionTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for stochastic_diffusion_test.
//
//  Discussion:
//
//    stochastic_diffusion_test tests STOCHASTIC_DIFFUSION.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2019
//
//  Author:
//
//    John Burkardt
//
        {
            Console.WriteLine("");
            Console.WriteLine("stochastic_diffusion_test");
            Console.WriteLine("  Test stochastic_diffusion.");

            diffusivity_1d_pwc_test();
            diffusivity_1d_xk_test();
            diffusivity_2d_bnt_contour();
            diffusivity_2d_elman_contour();
            diffusivity_2d_ntw_contour();
            diffusivity_2d_pwc_test();

            Console.WriteLine("");
            Console.WriteLine("stochastic_diffusion_test");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void diffusivity_1d_pwc_test()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_1d_pwc_test tests diffusivity_1d_pwc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2019
//
//  Author:
//
//    John Burkardt
//
        {
            string command_filename = "diffusivity_1d_pwc_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_1d_pwc_data.txt";
            List<string> data_unit = new List<string>();
            double[] vc =  {
                1.0, 1.5, 3.0, 1.2, 1.0, 0.8, 0.2, 0.4, 0.8, 1.0
            }
            ;
            double[] xc =  {
                -0.9, -0.5, -0.45, -0.1, 0.2, 0.3, 0.32, 0.7, 0.85
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("diffusivity_1d_pwc_test");
            Console.WriteLine("  Test diffusivity_1d_pwc.");
//
//  Set up the spatial grid.
//
            int nc = 10;
//
//  Sample the diffusivity.
//
            int np = 100;
            double x_min = -1.0;
            double x_max = +1.0;
            double[] xp = typeMethods.r8vec_linspace_new(np, x_min, x_max);
//
//  Compute the diffusivity field.
//
            double[] vp = Stochastic.diffusivity_1d_pwc(nc, xc, vc, np, xp);
//
//  Create data file.
//
            for (int j = 0; j < np; j++)
            {
                data_unit.Add("  " + xp[j].ToString().PadLeft(14)
                    + "  " + vp[j].ToString().PadLeft(14) + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
//
//  Create the command file.
//
            double vp_max = typeMethods.r8vec_max(np, vp);

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_1d_pwc.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Rho(X)--->'");
            command_unit.Add("set yrange [0.0:" + vp_max + "]");
            command_unit.Add("set title 'PWC 1D Stochastic diffusivity function'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3 linecolor rgb 'red'");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        }

        static void diffusivity_1d_xk_test()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_1d_xk_test tests diffusivity_1d_xk.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu, George Karniadakis,
//    Modeling uncertainty in steady state diffusion problems via
//    generalized polynomial chaos,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 191, 2002, pages 4927-4948.
//
        {
            string command_filename = "diffusivity_1d_xk_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_1d_xk_data.txt";
            List<string> data_unit = new List<string>();

            Console.WriteLine("");
            Console.WriteLine("diffusivity_1d_xk_test");
            Console.WriteLine("  Plot the stochastic diffusivity function");
            Console.WriteLine("  defined by diffusivity_1d_xk.");
//
//  Set up the spatial grid.
//
            int n = 51;
            double x_min = -1.0;
            double x_max = +1.0;
            double[] x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
//
//  Sample the OMEGA values.
//
            int m = 5;
            int seed = 123456789;
            double[] omega = typeMethods.r8vec_normal_01_new(m, ref seed);
//
//  Compute the diffusivity field.
//
            double dc0 = 10.0;
            double[] dc = Stochastic.diffusivity_1d_xk(dc0, m, omega, n, x);
//
//  Create data file.
//
            for (int j = 0; j < n; j++)
            {
                data_unit.Add("  " + x[j].ToString().PadLeft(14)
                    + "  " + dc[j].ToString().PadLeft(14) + "");
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
//
//  Create the command file.
//
            double dc_max = typeMethods.r8vec_max(n, dc);

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_1d_xk.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---DC(X)--->'");
            command_unit.Add("set yrange [0.0:" + dc_max + "]");
            command_unit.Add("set title 'XK Stochastic diffusivity function'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            command_unit.Add("plot '" + data_filename + "' using 1:2 lw 3 linecolor rgb 'red'");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        }

        static void diffusivity_2d_bnt_contour()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_2d_bnt_contour displays contour plots of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is computed by diffusivity_2d_bnt.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ivo Babuska, Fabio Nobile, Raul Tempone,
//    A stochastic collocation method for elliptic partial differential equations
//    with random input data,
//    SIAM Journal on Numerical Analysis,
//    Volume 45, Number 3, 2007, pages 1005-1034.
//
        {
            string command_filename = "diffusivity_2d_bnt_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_2d_bnt_data.txt";
            List<string> data_unit = new List<string>();
            int m = 4;
            int nx = 41;
            int ny = 31;

            Console.WriteLine("");
            Console.WriteLine("diffusivity_2d_bnt_contour");
            Console.WriteLine("  Display contour or surface plots of the stochastic");
            Console.WriteLine("  diffusivity function defined by diffusivity_2d_bnt.");
            Console.WriteLine("");
            Console.WriteLine("  The first plot uses uniform random values for OMEGA.");
            Console.WriteLine("  The second uses Gaussian (normal) random values.");
//
//  Set the spatial grid.
//
            double[] xvec = typeMethods.r8vec_linspace_new(nx, -1.5, 0.0);
            double[] yvec = typeMethods.r8vec_linspace_new(ny, -0.4, 0.8);

            double[] xmat = new double[nx * ny];
            double[] ymat = new double[nx * ny];
            typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, xmat, ymat);
//
//  Sample OMEGA.
//
            int seed = 123456789;
            double[] omega = UniformRNG.r8vec_uniform_01_new(m, ref seed);
//
//  Compute the diffusivity field.
//
            double dc0 = 10.0;
            int n = nx * ny;
            double[] dc = Stochastic.diffusivity_2d_bnt(dc0, omega, n, xmat, ymat);
//
//  Create a data file.
//
            for (int j = 0; j < ny; j++)
            {
                string dutemp = "";
                for (int i = 0; i < nx; i++)
                {
                    dutemp
                        += "  " + xmat[i + j * nx].ToString().PadLeft(14)
                        + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                        + "  " + dc[i + j * nx].ToString().PadLeft(14) + "";
                }

                data_unit.Add(dutemp);
            }

            File.WriteAllLines(data_filename, data_unit);
            
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'.");
//
//  Create the command file.
//

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_2d_bnt.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y--->'");
            command_unit.Add("set zlabel '<---DC(X,Y)--->'");
            command_unit.Add("set title 'BNT Stochastic diffusivity function'");
            command_unit.Add("set contour");
            command_unit.Add("set timestamp");
            command_unit.Add("set cntrparam levels 10");
            command_unit.Add("#set view map");
            command_unit.Add("set view 75, 75");
            command_unit.Add("unset key");
            command_unit.Add("splot '" + data_filename + "'");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        }

        static void diffusivity_2d_elman_contour()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_2d_elman_contour displays a contour plot of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is computed by DIFFUSIVITY_2D_ELMAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Darran Furnaval,
//    Solving the stochastic steady-state diffusion problem using multigrid,
//    IMA Journal on Numerical Analysis,
//    Volume 27, Number 4, 2007, pages 675-688.
//
        {
            string command_filename = "diffusivity_2d_elman_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_2d_elman_data.txt";
            List<string> data_unit = new List<string>();
            int m_1d = 5;
            int nx = 51;
            int ny = 51;

            Console.WriteLine("");
            Console.WriteLine("diffusivity_2d_elman_contour");
            Console.WriteLine("  Display contour or surface plots of the stochastic");
            Console.WriteLine("  diffusivity function defined by diffusivity_2d_elman.");
//
//  Set the spatial grid.
//
            double a = 1.0;
            double[] xvec = typeMethods.r8vec_linspace_new(nx, -a, a);
            double[] yvec = typeMethods.r8vec_linspace_new(ny, -a, a);

            double[] xmat = new double[nx * ny];
            double[] ymat = new double[nx * ny];
            typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, xmat, ymat);
//
//  Sample OMEGA.
//
            int seed = 123456789;
            double[] omega = typeMethods.r8vec_normal_01_new(m_1d * m_1d, ref seed);
//
//  Compute the diffusivity field.
//
            double cl = 0.1;
            double dc0 = 10.0;
            double[] dc = Stochastic.diffusivity_2d_elman(a, cl, dc0, m_1d, omega, nx, nx, xmat, ymat);
//
//  Create a data file.
//
            for (int j = 0; j < ny; j++)
            {
                string dutemp = "";
                for (int i = 0; i < nx; i++)
                {
                    dutemp +=
                        "  " + xmat[i + j * nx].ToString().PadLeft(14)
                        + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                        + "  " + dc[i + j * nx].ToString().PadLeft(14);
                }

                data_unit.Add(dutemp);
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
//
//  Create the command file.
//

            command_unit.Add( "# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_2d_elman.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y--->'");
            command_unit.Add("set zlabel '<---DC(X,Y)--->'");
            command_unit.Add("set title 'Elman Stochastic diffusivity function'");
            command_unit.Add("set contour");
            command_unit.Add("set timestamp");
            command_unit.Add("set cntrparam levels 10");
            command_unit.Add("#set view map");
            command_unit.Add("set view 75, 75");
            command_unit.Add("unset key");
            command_unit.Add("splot '" + data_filename + "'");

            File.WriteAllLines(command_filename, command_unit);
            
            Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        }

        static void diffusivity_2d_ntw_contour()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_2d_ntw_contour displays a contour plot of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is computed by diffusivity_2d_ntw.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
        {
            string command_filename = "diffusivity_2d_ntw_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_2d_ntw_data.txt";
            List<string> data_unit = new List<string>();
            int m = 21;
            int nx = 101;
            int ny = 101;

            Console.WriteLine("");
            Console.WriteLine("diffusivity_2d_ntw_contour");
            Console.WriteLine("  Display contour or surface plots of the stochastic");
            Console.WriteLine("  diffusivity function defined by diffusivity_2d_ntw.");
//
//  Set the spatial grid.
//
            double d = 1.0;
            double[] xvec = typeMethods.r8vec_linspace_new(nx, 0.0, d);
            double[] yvec = typeMethods.r8vec_linspace_new(ny, 0.0, d);

            double[] xmat = new double[nx * ny];
            double[] ymat = new double[nx * ny];
            typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, xmat, ymat);
//
//  Sample OMEGA.
//  We rescale to  [-sqrt(3),sqrt(3)].
//
            int seed = 123456789;
            double[] omega = UniformRNG.r8vec_uniform_01_new(m, ref seed);
            for (int i = 0; i < m; i++)
            {
                omega[i] = (1.0 - omega[i]) * (-Math.Sqrt(3.0))
                           + omega[i] * Math.Sqrt(3.0);
            }

//
//  Evaluate the diffusivity field.
//
            double cl = 0.1;
            double dc0 = 0.5;
            double[] dc = Stochastic.diffusivity_2d_ntw(cl, dc0, m, omega, nx * ny, xmat, ymat);
//
//  Create a data file.
//
            for (int j = 0; j < ny; j++)
            {
                string dutemp = "";
                for (int i = 0; i < nx; i++)
                {
                    dutemp
                        += "  " + xmat[i + j * nx].ToString().PadLeft(14)
                        + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                        + "  " + dc[i + j * nx].ToString().PadLeft(14) + "";
                }

                data_unit.Add(dutemp);
            }

            File.WriteAllLines(data_filename, data_unit);

            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
//
//  Create the command file.
//

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_2d_ntw.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y--->'");
            command_unit.Add("set zlabel '<---DC(X,Y)--->'");
            command_unit.Add("set title 'NTW Stochastic diffusivity function'");
            command_unit.Add("set contour");
            command_unit.Add("set timestamp");
            command_unit.Add("set cntrparam levels 15");
            command_unit.Add("#set view map");
            command_unit.Add("set view 65, 65");
            command_unit.Add("set key");
            command_unit.Add("splot '" + data_filename + "'");

            File.WriteAllLines(command_filename, command_unit);

            Console.WriteLine("  Created graphics command file '" + command_filename + "'.");

        }

        static void diffusivity_2d_pwc_test()
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_2d_pwc_test tests diffusivity_2d_pwc.
//
//  Discussion:
//
//    This function calls diffusivity_2d_pwc_contour, which evaluates
//    diffusivity_2d_pwc() over a given [a,b]x[c,d] grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2019
//
//  Author:
//
//    John Burkardt
//
        {
            double a;
            double b;
            double c;
            double d;
            int h;
            int w;

            h = 5;
            w = 4;
            a = -1.0;
            b = +1.0;
            c = -1.0;
            d = +1.0;
            diffusivity_2d_pwc_contour(h, w, a, b, c, d);
        }

        static void diffusivity_2d_pwc_contour(int h, int w, double a, double b, double c,
            double d)
//****************************************************************************80
//
//  Purpose:
//
//    diffusivity_2d_pwc_contour displays a contour plot of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is computed by diffusivity_2d_pwc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 March 2019
//
//  Author:
//
//    John Burkardt
//
        {
            string command_filename = "diffusivity_2d_pwc_commands.txt";
            List<string> command_unit = new List<string>();
            string data_filename = "diffusivity_2d_pwc_data.txt";
            List<string> data_unit = new List<string>();
            int n;
            int nx;
            int ny;
            double[] omega;
            double[] rho;
            int seed;
            double[] xmat;
            double[] xvec;
            double[] ymat;
            double[] yvec;

            Console.WriteLine("");
            Console.WriteLine("diffusivity_2d_pwc_contour");
            Console.WriteLine("  Display contour or surface plots of the stochastic");
            Console.WriteLine("  diffusivity function defined by diffusivity_2d_pwc.");
            Console.WriteLine("");
            Console.WriteLine("  Underlying grid is " + w + " elements wide (X) and " + h + " high (Y)");

            nx = 101;
            ny = 101;
//
//  Set the spatial grid.
//
            xvec = typeMethods.r8vec_linspace_new(nx, a, b);
            yvec = typeMethods.r8vec_linspace_new(ny, c, d);

            xmat = new double[nx * ny];
            ymat = new double[nx * ny];
            typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, xmat, ymat);
//
//  Sample OMEGA.
//
            seed = 123456789;
            omega = UniformRNG.r8vec_uniform_ab_new(h * w, 0.5, 1.5, ref seed);

            n = nx * ny;
//
//  Compute the diffusivity field.
//
            rho = Stochastic.diffusivity_2d_pwc(h, w, a, b, c, d, omega, n, xmat, ymat);
//
//  Create a data file.
//
            for (int j = 0; j < ny; j++)
            {
                string dutemp = "";
                for (int i = 0; i < nx; i++)
                {
                    dutemp
                        += "  " + xmat[i + j * nx].ToString().PadLeft(14)
                        + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                        + "  " + rho[i + j * nx].ToString().PadLeft(14) + "";
                }

                data_unit.Add(dutemp);
            }

            File.WriteAllLines(data_filename, data_unit);
            
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'");
//
//  Create the command file.
//

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set output 'diffusivity_2d_pwc.png'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---Y--->'");
            command_unit.Add("set zlabel '<---Rho(X,Y)--->'");
            command_unit.Add("set title 'PWC Stochastic diffusivity function'");
            command_unit.Add("set contour");
            command_unit.Add("set timestamp");
            command_unit.Add("set cntrparam levels 10");
            command_unit.Add("#set view map");
            command_unit.Add("set view 75, 75");
            command_unit.Add("unset key");
            command_unit.Add("splot '" + data_filename + "'");

            File.WriteAllLines(command_filename,command_unit);
            
            Console.WriteLine("  Created graphics command file '" + command_filename + "'");

        }
    }
}