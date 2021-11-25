using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.PiecewiseLinear;
using Burkardt.Types;
using InterpTest;

namespace PWLInterp1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PWL_INTERP_1D_TEST.
        //
        //  Discusion:
        //
        //    PWL_INTERP_1D_TEST tests the PWL_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int prob;

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_1D_TEST:");
        Console.WriteLine("  Test the PWL_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The test needs the TEST_INTERP library.");

        pwl_basis_1d_test();

        pwl_value_1d_test();

        int prob_num = TestInterp.p00_prob_num();
        for (prob = 1; prob <= prob_num; prob++)
        {
            pwl_interp_1d_test01(prob);
        }

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void pwl_basis_1d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_BASIS_1D_TEST tests PWL_BASIS_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nd = 4;
        const int ni = 21;
        double[] xd = { 0.0, 2.0, 5.0, 10.0 };

        Console.WriteLine("");
        Console.WriteLine("PWL_BASIS_1D_TEST:");
        Console.WriteLine("  PWL_BASIS_1D evaluates the piecewise linear 1D basis");
        Console.WriteLine("  functions.");

        const double x_min = 0.0;
        const double x_max = 10.0;
        double[] xi = typeMethods.r8vec_linspace_new(ni, x_min, x_max);

        double[] lb = Interp1D.pwl_basis_1d(nd, xd, ni, xi);

        typeMethods.r8mat_print(ni, nd, lb, "  The piecewise linear basis functions:");
    }

    private static void pwl_value_1d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_VALUE_1D_TEST tests PWL_VALUE_1D.
        //
        //  Discussion:
        //
        //    f(x) = x^3 - 12 x^2 + 39 x - 28 = ( x - 1 ) * ( x - 4 ) * ( x - 7 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified: 
        //
        //    01 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nd = 4;
        const int ni = 21;
        double[] xd = { 0.0, 2.0, 5.0, 10.0 };
        double[] yd = { -28.0, +10.0, -8.0, +162.0 };

        Console.WriteLine("");
        Console.WriteLine("PWL_VALUE_1D_TEST:");
        Console.WriteLine("  PWL_VALUE_1D evaluates a piecewise linear 1D interpolant.");

        const double x_min = 0.0;
        const double x_max = 10.0;
        double[] xi = typeMethods.r8vec_linspace_new(ni, x_min, x_max);

        double[] yi = Interp1D.pwl_value_1d(nd, xd, yd, ni, xi);

        typeMethods.r8vec2_print(ni, xi, yi, "  Table of interpolant values:");
    }

    private static void pwl_interp_1d_test01(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_INTERP_1D_TEST01 tests PWL_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        List<string> interp_unit = new();
        int j;

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_1D_TEST01:");
        Console.WriteLine("  PWL_INTERP_1D evaluates the piecewise linear interpolant.");
        Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");

        int nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xy = TestInterp.p00_data(prob, 2, nd);

        typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");

        double[] xd = new double[nd];
        double[] yd = new double[nd];

        for (i = 0; i < nd; i++)
        {
            xd[i] = xy[0 + 2 * i];
            yd[i] = xy[1 + 2 * i];
        }

        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;

        double[] xi = new double[ni];
        for (i = 0; i < ni; i++)
        {
            xi[i] = xd[i];
        }

        double[] yi = Interp1D.pwl_value_1d(nd, xd, yd, ni, xi);

        double interp_error = typeMethods.r8vec_diff_norm(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + interp_error + "");
        //
        //  Create data file.
        //
        string data_filename = "data" + prob + ".txt";
        for (j = 0; j < nd; j++)
        {
            data_unit.Add("  " + xd[j]
                               + "  " + yd[j] + "");
        }

        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file \"" + data_filename + "\".");
        //
        //  Create interp file.
        //

        ni = 501;
        double xmin = typeMethods.r8vec_min(nd, xd);
        double xmax = typeMethods.r8vec_max(nd, xd);
        xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
        yi = Interp1D.pwl_value_1d(nd, xd, yd, ni, xi);

        string interp_filename = "interp" + prob + ".txt";
        for (j = 0; j < ni; j++)
        {
            interp_unit.Add("  " + xi[j]
                                 + "  " + yi[j] + "");
        }

        File.WriteAllLines(interp_filename, interp_unit);

        Console.WriteLine("  Created graphics interp file \"" + interp_filename + "\".");
        //
        //  Plot the data and the interpolant.
        //
        string command_filename = "commands" + prob + ".txt";

        string output_filename = "plot" + prob + ".png";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---Y--->'");
        command_unit.Add("set title 'Data versus Nearest Neighbor Interpolant'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("plot '" + data_filename
                                  + "' using 1:2 with points pt 7 ps 2 lc rgb 'blue',\\");
        command_unit.Add("     '" + interp_filename
                                  + "' using 1:2 lw 3 linecolor rgb 'red'");

        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Created graphics command file \"" + command_filename + "\".");
    }
}