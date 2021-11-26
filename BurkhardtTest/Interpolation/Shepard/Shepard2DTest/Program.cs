﻿using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace Shepard2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SHEPARD_INTERP_2D_TEST.
        //
        //  Discussion:
        //
        //    SHEPARD_INTERP_2D_TEST tests the SHEPARD_INTERP_2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] p_test =
            {
                1.0, 2.0, 4.0, 8.0
            }
            ;
        const int p_test_num = 4;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_2D_TEST:");
        Console.WriteLine("  Test the SHEPARD_INTERP_2D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();
        int g = 1;

        for (prob = 1; prob <= prob_num; prob++)
        {
            int j;
            for (j = 0; j < p_test_num; j++)
            {
                double p = p_test[j];
                test01(prob, g, p);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int g, double p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SHEPARD_INTERP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int G, the grid number.
        //
        //    Input, double P, the power used in the distance weighting.
        //
    {
        const bool debug = false;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
        Console.WriteLine("  using grid #" + g + "");
        Console.WriteLine("  using Shepard interpolation with P = " + p + "");

        int nd = Data_2D.g00_size(g);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xd = new double[nd];
        double[] yd = new double[nd];
        Data_2D.g00_xy(g, nd, ref xd, ref yd);

        double[] zd = new double[nd];
        Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

        switch (debug)
        {
            case true:
                typeMethods.r8vec3_print(nd, xd, yd, zd, "  X, Y, Z data:");
                break;
        }

        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = typeMethods.r8vec_copy_new(ni, yd);

        double[] zi = Shepard.shepard_interp_2d(nd, xd, yd, zd, p, ni, xi, yi);

        double int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

    }
}