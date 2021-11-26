﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace Newton1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEWTON_INTERP_1D_TEST tests the NEWTON_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int prob;

        Console.WriteLine("");
        Console.WriteLine("NEWTON_INTERP_1D_TEST:");
        Console.WriteLine("  Test the NEWTON_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test needs the TEST_INTERP library as well.");

        newton_coef_1d_test();

        newton_value_1d_test();

        int prob_num = Data_1D.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            newton_interp_1d_test01(prob);
        }

        Console.WriteLine("");
        Console.WriteLine("NEWTON_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void newton_coef_1d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEWTON_COEF_1D_TEST tests NEWTON_COEF_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nd = 5;
        double[] xd = {0.0, 1.0, 2.0, 3.0, 4.0};
        double[] yd = {24.0, 0.0, 0.0, 0.0, 0.0};

        Console.WriteLine("");
        Console.WriteLine("NEWTON_COEF_1D_TEST");
        Console.WriteLine("  NEWTON_COEF_1D sets the coefficients for a 1D Newton interpolation.");

        typeMethods.r8vec2_print(nd, xd, yd, "  Interpolation data:");

        double[] cd = Newton1D.newton_coef_1d(nd, xd, yd);

        typeMethods.r8vec_print(nd, cd, "  Newton interpolant coefficients:");
    }

    private static void newton_value_1d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEWTON_VALUE_1D_TEST tests NEWTON_VALUE_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] cd = {24.0, -24.0, +12.0, -4.0, 1.0};
        const int nd = 5;
        const int ni = 16;
        double[] xd = {0.0, 1.0, 2.0, 3.0, 4.0};

        Console.WriteLine("");
        Console.WriteLine("NEWTON_VALUE_1D_TEST");
        Console.WriteLine("  NEWTON_VALUE_1D evaluates a Newton 1d interpolant.");

        typeMethods.r8vec2_print(nd, xd, cd, "  The Newton interpolant data:");

        const double x_lo = 0.0;
        const double x_hi = 5.0;
        double[] xi = typeMethods.r8vec_linspace_new(ni, x_lo, x_hi);
        double[] yi = Newton1D.newton_value_1d(nd, xd, cd, ni, xi);

        typeMethods.r8vec2_print(ni, xi, yi, "  Newton interpolant values:");
    }

    private static void newton_interp_1d_test01(int prob)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEWTON_INTERP_1D_TEST01 tests NEWTON_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        List<string> interp_unit = new();
        int j;

        Console.WriteLine("");
        Console.WriteLine("NEWTON_INTERP_1D_TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");

        int nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xy = TestInterp.p00_data(prob, 2, nd);

        double[] xd = new double[nd];
        double[] yd = new double[nd];

        for (i = 0; i < nd; i++)
        {
            xd[i] = xy[0 + i * 2];
            yd[i] = xy[1 + i * 2];
        }

        typeMethods.r8vec2_print(nd, xd, yd, "  X, Y data:");
        //
        //  Get the Newton coefficients.
        //
        double[] cd = Newton1D.newton_coef_1d(nd, xd, yd);
        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = Newton1D.newton_value_1d(nd, xd, cd, ni, xi);

        double interp_error = typeMethods.r8vec_norm_affine(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + interp_error + "");

        //
        //  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
        //  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
        //  (YMAX-YMIN).
        //
        double xmin = typeMethods.r8vec_min(nd, xd);
        double xmax = typeMethods.r8vec_max(nd, xd);
        double ymin = typeMethods.r8vec_min(nd, yd);
        double ymax = typeMethods.r8vec_max(nd, yd);

        ni = 501;
        xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
        yi = Newton1D.newton_value_1d(nd, xd, cd, ni, xi);

        double ld = 0.0;
        for (i = 0; i < nd - 1; i++)
        {
            ld += Math.Sqrt(Math.Pow((xd[i + 1] - xd[i]) / (xmax - xmin), 2)
                            + Math.Pow((yd[i + 1] - yd[i]) / (ymax - ymin), 2));
        }

        double li = 0.0;
        for (i = 0; i < ni - 1; i++)
        {
            li += Math.Sqrt(Math.Pow((xi[i + 1] - xi[i]) / (xmax - xmin), 2)
                            + Math.Pow((yi[i + 1] - yi[i]) / (ymax - ymin), 2));
        }

        Console.WriteLine("");
        Console.WriteLine("  Normalized length of piecewise linear interpolant = " + ld + "");
        Console.WriteLine("  Normalized length of Newton interpolant           = " + li + "");

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
        xmin = typeMethods.r8vec_min(nd, xd);
        xmax = typeMethods.r8vec_max(nd, xd);
        xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
        yi = Newton1D.newton_value_1d(nd, xd, cd, ni, xi);

        string interp_filename = "interp" + prob + ".txt";
        if (interp_filename == null)
        {
            throw new ArgumentNullException(nameof(interp_filename));
        }

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
        command_unit.Add("set title 'Data versus Newton polynomial interpolant'");
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