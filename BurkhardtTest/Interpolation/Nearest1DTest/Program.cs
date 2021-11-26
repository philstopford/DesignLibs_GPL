using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace Nearest1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for NEAREST_INTERP_1D_TEST.
        //
        //  Discussion:
        //
        //    NEAREST_INTERP_1D_TEST tests the NEAREST_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 May 2013
        //
        //  Author:
        //
        //   John Burkardt
        //
    {
        int prob;

        Console.WriteLine("");
        Console.WriteLine("NEAREST_INTERP_1D_TEST:");
        Console.WriteLine("  Test the NEAREST_INTERP_1D library.");
        Console.WriteLine("  The test needs the TEST_INTERP library.");

        int prob_num = Data_1D.p00_prob_num();

        int ni = 11;
        for (prob = 1; prob <= prob_num; prob++)
        {
            nearest_interp_1d_test01(prob, ni);
        }

        for (prob = 1; prob <= prob_num; prob++)
        {
            nearest_interp_1d_test02(prob);
        }

        Console.WriteLine("");
        Console.WriteLine("NEAREST_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void nearest_interp_1d_test01(int prob, int ni)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEAREST_INTERP_1D_TEST01 tests NEAREST_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the index of the problem.
        //
        //    Input, int NI, the number of interpolation points.
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("NEAREST_INTERP_1D_TEST01");
        Console.WriteLine("  Sample the nearest neighbor interpolant for problem # " + prob + "");

        int nd = TestInterp.p00_data_num(prob);

        double[] d = TestInterp.p00_data(prob, 2, nd);

        double[] xd = new double[nd];
        double[] yd = new double[nd];

        for (j = 0; j < nd; j++)
        {
            xd[j] = d[0 + j * 2];
            yd[j] = d[1 + j * 2];
        }

        double xd_min = typeMethods.r8vec_min(nd, xd);
        double xd_max = typeMethods.r8vec_max(nd, xd);

        double[] xi = typeMethods.r8vec_linspace_new(ni, xd_min, xd_max);
        double[] yi = Nearest1D.nearest_interp_1d(nd, xd, yd, ni, xi);

        string title = "X, Y for problem " + prob;

        typeMethods.r8vec2_print(ni, xi, yi, title);
    }

    private static void nearest_interp_1d_test02(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEAREST_INTERP_1D_TEST02 tests NEAREST_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the index of the problem.
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        List<string> interp_unit = new();
        int j;

        Console.WriteLine("");
        Console.WriteLine("NEAREST_INTERP_1D_TEST02:");
        Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");

        int nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xy = TestInterp.p00_data(prob, 2, nd);

        typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");

        double[] xd = new double[nd];
        double[] yd = new double[nd];
        for (i = 0; i < nd; i++)
        {
            xd[i] = xy[0 + i * 2];
            yd[i] = xy[1 + i * 2];
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

        double[] yi = Nearest1D.nearest_interp_1d(nd, xd, yd, ni, xi);

        double interp_error = typeMethods.r8vec_diff_norm(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  Node-averaged L2 interpolation error = " + interp_error + "");
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
        yi = Nearest1D.nearest_interp_1d(nd, xd, yd, ni, xi);

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