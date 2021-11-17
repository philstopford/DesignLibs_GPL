using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.StochasticDifferentialEquations;
using Burkardt.Types;
using Burkardt.Uniform;

namespace StochasticDiffusionTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STOCHASTIC_DIFFUSION_TEST.
        //
        //  Discussion:
        //
        //    STOCHASTIC_DIFFUSION_TEST tests the STOCHASTIC_DIFFUSION library.
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
    {
        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_DIFFUSION_TEST");
        Console.WriteLine("  Test the STOCHASTIC_DIFFUSION library.");

        bnt_contour();
        elman_contour();
        ntw_contour();
        xk_contour();

        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_DIFFUSION_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void bnt_contour()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BNT_CONTOUR displays contour plots of a 2D stochastic diffusivity function.
        //
        //  Discussion:
        //
        //    The diffusivity function is compute by DIFFUSIVITY_2D_BNT.
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
        string command_filename = "bnt_commands.txt";
        List<string> command_unit = new();
        string data_filename = "bnt_data.txt";
        List<string> data_unit = new();
        double[] dc;
        double dc0;
        int i;
        int j;
        int m = 4;
        int n;
        int nx = 41;
        int ny = 31;
        double[] omega;
        int seed;
        double[] xmat;
        double[] xvec;
        double[] ymat;
        double[] yvec;

        Console.WriteLine("");
        Console.WriteLine("BNT_CONTOUR");
        Console.WriteLine("  Display contour or surface plots of the stochastic");
        Console.WriteLine("  diffusivity function defined by DIFFUSIVITY_2D_BNT.");
        Console.WriteLine("");
        Console.WriteLine("  The first plot uses uniform random values for OMEGA.");
        Console.WriteLine("  The second uses Gaussian (normal) random values.");
        //
        //  Set the spatial grid.
        //
        xvec = typeMethods.r8vec_linspace_new(nx, -1.5, 0.0);
        yvec = typeMethods.r8vec_linspace_new(ny, -0.4, 0.8);

        xmat = new double[nx * ny];
        ymat = new double[nx * ny];
        typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, ref xmat, ref ymat);
        //
        //  Sample OMEGA.
        //
        seed = 123456789;
        omega = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        //
        //  Compute the diffusivity field.
        //
        dc0 = 10.0;
        n = nx * ny;
        dc = Diffusion.diffusivity_2d_bnt(dc0, omega, n, xmat, ymat);

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                data_unit.Add("  " + xmat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + dc[i + j * nx].ToString().PadLeft(14) + "");
            }

            data_unit.Add("");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file '" + data_filename + "'.");

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output 'bnt_contour.png'");
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

    private static void elman_contour()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELMAN_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
        //
        //  Discussion:
        //
        //    The diffusivity function is compute by DIFFUSIVITY_2D_ELMAN.
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
        double a;
        double cl;
        string command_filename = "elman_commands.txt";
        List<string> command_unit = new();
        string data_filename = "elman_data.txt";
        List<string> data_unit = new();
        double[] dc;
        double dc0;
        int i;
        int j;
        int m_1d = 5;
        int nx = 51;
        int ny = 51;
        double[] omega;
        int seed;
        double[] xmat;
        double[] xvec;
        double[] ymat;
        double[] yvec;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("ELMAN_CONTOUR");
        Console.WriteLine("  Display contour or surface plots of the stochastic");
        Console.WriteLine("  diffusivity function defined by DIFFUSIVITY_2D_ELMAN.");
        //
        //  Set the spatial grid.
        //
        a = 1.0;
        xvec = typeMethods.r8vec_linspace_new(nx, -a, a);
        yvec = typeMethods.r8vec_linspace_new(ny, -a, a);

        xmat = new double[nx * ny];
        ymat = new double[nx * ny];
        typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, ref xmat, ref ymat);
        //
        //  Sample OMEGA.
        //
        seed = 123456789;
        omega = typeMethods.r8vec_normal_01_new(m_1d * m_1d, ref data, ref seed);
        //
        //  Compute the diffusivity field.
        //
        cl = 0.1;
        dc0 = 10.0;
        dc = Diffusion.diffusivity_2d_elman(a, cl, dc0, m_1d, omega, nx, nx, xmat, ymat);

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                data_unit.Add("  " + xmat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + dc[i + j * nx].ToString().PadLeft(14) + "");
            }

            data_unit.Add("");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file '" + data_filename + "'");

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output 'elman_contour.png'");
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

    private static void ntw_contour()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NTW_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
        //
        //  Discussion:
        //
        //    The diffusivity function is compute by DIFFUSIVITY_2D_NTW.
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
        double cl;
        string command_filename = "ntw_commands.txt";
        List<string> command_unit = new();
        double d;
        string data_filename = "ntw_data.txt";
        List<string> data_unit = new();
        double[] dc;
        double dc0;
        int i;
        int j;
        int m = 21;
        int nx = 101;
        int ny = 101;
        double[] omega;
        int seed;
        double[] xmat;
        double[] xvec;
        double[] ymat;
        double[] yvec;

        Console.WriteLine("");
        Console.WriteLine("NTW_CONTOUR");
        Console.WriteLine("  Display contour or surface plots of the stochastic");
        Console.WriteLine("  diffusivity function defined by DIFFUSIVITY_2D_NTW.");
        //
        //  Set the spatial grid.
        //
        d = 1.0;
        xvec = typeMethods.r8vec_linspace_new(nx, 0.0, d);
        yvec = typeMethods.r8vec_linspace_new(ny, 0.0, d);

        xmat = new double[nx * ny];
        ymat = new double[nx * ny];
        typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, ref xmat, ref ymat);
        //
        //  Sample OMEGA.
        //  We rescale to  [-sqrt(3),sqrt(3)].
        //
        seed = 123456789;
        omega = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        for (i = 0; i < m; i++)
        {
            omega[i] = (1.0 - omega[i]) * -Math.Sqrt(3.0)
                       + omega[i] * Math.Sqrt(3.0);
        }

        //
        //  Evaluate the diffusivity field.
        //
        cl = 0.1;
        dc0 = 0.5;
        dc = Diffusion.diffusivity_2d_ntw(cl, dc0, m, omega, nx * ny, xmat, ymat);

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                data_unit.Add("  " + xmat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + ymat[i + j * nx].ToString().PadLeft(14)
                                   + "  " + dc[i + j * nx].ToString().PadLeft(14) + "");
            }

            data_unit.Add("");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file '" + data_filename + "'");

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output 'ntw_contour.png'");
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

    private static void xk_contour()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XK_CONTOUR displays contour plots of a 1D stochastic diffusivity function.
        //
        //  Discussion:
        //
        //    The diffusivity function is compute by DIFFUSIVITY_1D_XK.
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
        double[] dc;
        double dc_max;
        double dc0;
        string command_filename = "xk_commands.txt";
        List<string> command_unit = new();
        string data_filename = "xk_data.txt";
        List<string> data_unit = new();
        int j;
        int m;
        int n;
        int seed;
        double[] omega;
        double[] x;
        double x_min;
        double x_max;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("XK_CONTOUR");
        Console.WriteLine("  Plot the stochastic diffusivity function");
        Console.WriteLine("  defined by DIFFUSIVITY_1D_XK.");
        //
        //  Set up the spatial grid.
        //
        n = 51;
        x_min = -1.0;
        x_max = +1.0;
        x = typeMethods.r8vec_linspace_new(n, x_min, x_max);
        //
        //  Sample the OMEGA values.
        //
        m = 5;
        seed = 123456789;
        omega = typeMethods.r8vec_normal_01_new(m, ref data, ref seed);
        //
        //  Compute the diffusivity field.
        //
        dc0 = 10.0;
        dc = Diffusion.diffusivity_1d_xk(dc0, m, omega, n, x);

        for (j = 0; j < n; j++)
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
        dc_max = typeMethods.r8vec_max(n, dc);

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output 'xk_contour.png'");
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
}