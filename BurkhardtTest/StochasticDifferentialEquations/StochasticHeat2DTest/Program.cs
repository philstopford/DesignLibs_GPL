using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.StochasticDifferentialEquations;
using Burkardt.Types;

namespace StochasticHeat2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STOCHASTIC_HEAT2D_TEST.
        //
        //  Discussion:
        //
        //    STOCHASTIC_HEAT2D_TEST tests the STOCHASTIC_HEAT2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_HEAT2D_TEST:");
        Console.WriteLine("  Test the STOCHASTIC_HEAT2D library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("STOCHASTIC_HEAT2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 plots a sample solution of a 2D stochastic diffusivity equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string command_filename = "solution_commands.txt";
        List<string> command_unit = new();
        string data_filename = "solution_data.txt";
        List<string> data_unit = new();
        int i;
        int j;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Consider the steady heat equation in the unit square,");
        Console.WriteLine("  with 0 Dirichlet boundary conditions, ");
        Console.WriteLine("  and a heat source term F that is a Gaussian centered at (0.60,0.80).");
        Console.WriteLine("");
        Console.WriteLine("  Model the diffusivity coefficient as spatially varying,");
        Console.WriteLine("  with a stochastic dependence on parameters OMEGA(1:4),");
        Console.WriteLine("  as described in Babuska, Nobile, Tempone (BNT).");
        Console.WriteLine("");
        Console.WriteLine("  Compute and display the solution U for a given choice");
        Console.WriteLine("  of the parameters OMEGA.");
        //
        //  Create the X and Y coordinate vectors.
        //
        int nx = 21;
        double xmin = 0.0;
        double xmax = 1.0;
        double[] xvec = typeMethods.r8vec_linspace_new(nx, xmin, xmax);

        int ny = 21;
        double ymin = 0.0;
        double ymax = 1.0;
        double[] yvec = typeMethods.r8vec_linspace_new(ny, ymin, ymax);
        //
        //  Create the X and Y coordinate matrices.
        //
        double[] xmat = new double[nx * ny];
        double[] ymat = new double[nx * ny];
        typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, ref xmat, ref ymat);
        //
        //  Sample OMEGA:
        //
        int seed = 123456789;
        double[] omega = typeMethods.r8vec_normal_01_new(4, ref data, ref seed);
        for (i = 0; i < 4; i++)
        {
            omega[i] = 2.0 * omega[i];
        }

        typeMethods.r8vec_print(4, omega, "  Sampled OMEGA values:");
        //
        //  Solve the finite difference approximation to the steady 2D heat equation
        //  for this set of OMEGA values.
        //
        double[] umat = Diffusion.stochastic_heat2d(omega, nx, ny, xvec, yvec, test01_f);

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                data_unit.Add("  " + xmat[i + j * nx]
                                   + "  " + ymat[i + j * nx]
                                   + "  " + umat[i + j * nx] + "");
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
        command_unit.Add("set output 'solution.png'");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---Y--->'");
        command_unit.Add("set zlabel '<---U(X,Y)--->'");
        command_unit.Add("set title 'Sample Solution'");
        command_unit.Add("set contour");
        command_unit.Add("set timestamp");
        command_unit.Add("set cntrparam levels 10");
        command_unit.Add("set view 75, 75");
        command_unit.Add("unset key");
        command_unit.Add("splot '" + data_filename + "'");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        //
        //  Report the average value of U.
        //
        double u_mean = typeMethods.r8mat_mean(nx, ny, umat);

        Console.WriteLine("");
        Console.WriteLine("  Mean value of U is " + u_mean + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 looks at mean temperature as a function of OMEGA(1) and OMEGA(2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string command_filename = "umean_commands.txt";
        List<string> command_unit = new();
        string data_filename = "umean_data.txt";
        List<string> data_unit = new();
        int i;
        int j;
        double[] omega = new double[4];

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Fix OMEGA(3) = 4, OMEGA(4) = 0, and");
        Console.WriteLine("  examine dependence of average temperature on OMEGA(1) and OMEGA(2)");
        Console.WriteLine("  over the range [-10,+10].");
        //
        //  Create the X and Y coordinate vectors.
        //
        int nx = 21;
        double xmin = 0.0;
        double xmax = 1.0;
        double[] xvec = typeMethods.r8vec_linspace_new(nx, xmin, xmax);

        int ny = 21;
        double ymin = 0.0;
        double ymax = 1.0;
        double[] yvec = typeMethods.r8vec_linspace_new(ny, ymin, ymax);
        //
        //  Create the X and Y coordinate matrices.
        //
        double[] xmat = new double[nx * ny];
        double[] ymat = new double[nx * ny];
        typeMethods.r8vec_mesh_2d(nx, ny, xvec, yvec, ref xmat, ref ymat);
        //
        //  Create OMEGA1 and OMEGA2 vectors.
        //
        int omega1_num = 21;
        double omega1_min = -10.0;
        double omega1_max = +10.0;
        double[] omega1_vec = typeMethods.r8vec_linspace_new(omega1_num, omega1_min, omega1_max);

        int omega2_num = 21;
        double omega2_min = -10.0;
        double omega2_max = +10.0;
        double[] omega2_vec = typeMethods.r8vec_linspace_new(omega2_num, omega2_min, omega2_max);
        //
        //  Create the OMEGA1 and OMEGA2 coordinate matrices.
        //
        double[] omega1_mat = new double[omega1_num * omega2_num];
        double[] omega2_mat = new double[omega1_num * omega2_num];
        typeMethods.r8vec_mesh_2d(omega1_num, omega2_num, omega1_vec, omega2_vec, ref omega1_mat, ref omega2_mat);
        //
        //  Set OMEGA(3) and OMEGA(4).
        //
        omega[2] = 4.0;
        omega[3] = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  Omega(3) fixed at " + omega[2] + "");
        Console.WriteLine("  Omega(4) fixed at " + omega[3] + "");
        //
        //  Solve the finite difference approximation to the steady 2D heat equation,
        //  and save the mean value of the solution, which is a slightly biased
        //  estimate of the heat integral over the unit square.
        //
        double[] u_mean_mat = new double[omega1_num * omega2_num];

        for (j = 0; j < omega2_num; j++)
        {
            omega[1] = omega2_vec[j];
            for (i = 0; i < omega1_num; i++)
            {
                omega[0] = omega1_vec[i];
                double[] umat = Diffusion.stochastic_heat2d(omega, nx, ny, xvec, yvec, test01_f);
                u_mean_mat[i + j * omega1_num] = typeMethods.r8mat_mean(nx, ny, umat);
            }
        }

        for (j = 0; j < ny; j++)
        {
            for (i = 0; i < nx; i++)
            {
                data_unit.Add("  " + omega1_mat[i + j * omega1_num]
                                   + "  " + omega2_mat[i + j * omega1_num]
                                   + "  " + u_mean_mat[i + j * omega1_num] + "");
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
        command_unit.Add("set output 'umean.png'");
        command_unit.Add("set xlabel '<---OMEGA1--->'");
        command_unit.Add("set ylabel '<---OMEGA2--->'");
        command_unit.Add("set zlabel '<---U_MEAN(OMEGA1,OMEGA2)--->'");
        command_unit.Add("set title 'Solution Mean as Function of Omega1, Omega2'");
        command_unit.Add("set contour");
        command_unit.Add("set timestamp");
        command_unit.Add("set cntrparam levels 10");
        command_unit.Add("set view 75, 75");
        command_unit.Add("unset key");
        command_unit.Add("splot '" + data_filename + "'");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Created graphics command file '" + command_filename + "'");
        //
        //  Print the maximum value of the mean.
        //
        double u_mean_max = typeMethods.r8mat_max(omega1_num, omega2_num, u_mean_mat);

        Console.WriteLine("");
        Console.WriteLine("  U_Mean_Max = " + u_mean_max + "");

    }


    private static double test01_f(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01_F evaluates the heat source term.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the evaluation point.
        //
        //    Output, double TEST01_F, the value of the heat source term at (X,Y).
        //
    {
        double v = 0.05;
        double arg = (Math.Pow(x - 0.60, 2) + Math.Pow(y - 0.80, 2)) / Math.Pow(v, 2);
        double value = 2000.0 * Math.Exp(-arg);

        return value;
    }
}