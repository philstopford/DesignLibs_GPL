﻿using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace Shepard1DTest;

internal static class Program
{
    private static void Main()
/******************************************************************************/
//
//  Purpose:
//
//    MAIN is the main program for SHEPARD_INTERP_1D_TEST.
//
//  Discussion:
//
//    SHEPARD_INTERP_1D_TEST tests the SHEPARD_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2015
//
//  Author:
//
//    John Burkardt
//
    {
        int p_num = 5;
        double[] p_test =  {
                0.0, 1.0, 2.0, 4.0, 8.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_1D_TEST:");
        Console.WriteLine("  Test the SHEPARD_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test needs the TEST_INTERP library as well.");

        shepard_basis_1d_test();

        shepard_value_1d_test();

        int prob_num = TestInterp.p00_prob_num();

        for (int prob = 1; prob <= prob_num; prob++)
        {
            for (int j = 0; j < p_num; j++)
            {
                double p = p_test[j];
                test01(prob, p);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void shepard_basis_1d_test()
//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_BASIS_1D_TEST tests SHEPARD_BASIS_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2015
//
//  Author:
//
//    John Burkardt
//
    {
        const int nd = 4;
        const int ni = 21;
        const double p = 2.0;
        double[] xd =  {
                0.0, 2.0, 5.0, 10.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_BASIS_1D_TEST:");
        Console.WriteLine("  SHEPARD_BASIS_1D evaluates the Shepard 1D basis");
        Console.WriteLine("  functions.");
        Console.WriteLine("");
        Console.WriteLine("  Using power P = " + p + "");

        const double x_min = 0.0;
        const double x_max = 10.0;
        double[] xi = typeMethods.r8vec_linspace_new(ni, x_min, x_max);

        double[] lb = Shepard.shepard_basis_1d(nd, xd, p, ni, xi);

        typeMethods.r8mat_print(ni, nd, lb, "  The Shepard basis functions:");

    }

    private static void shepard_value_1d_test()

//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_VALUE_1D_TEST tests SHEPARD_VALUE_1D.
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
//    03 July 2015
//
//  Author:
//
//    John Burkardt
//
    {
        const int nd = 4;
        const int ni = 21;
        double[] xd =  {
                0.0, 2.0, 5.0, 10.0
            }
            ;
        double[] yd =  {
                -28.0, +10.0, -8.0, +162.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_VALUE_1D_TEST:");
        Console.WriteLine("  SHEPARD_VALUE_1D evaluates a Shepard 1D interpolant.");

        const double p = 2.0;
        Console.WriteLine("");
        Console.WriteLine("  Using power P = " + p + "");

        const double x_min = 0.0;
        const double x_max = 10.0;
        double[] xi = typeMethods.r8vec_linspace_new(ni, x_min, x_max);

        double[] yi = Shepard.shepard_value_1d(nd, xd, yd, p, ni, xi);

        typeMethods.r8vec2_print(ni, xi, yi, "  Table of interpolant values:");

    }

    private static void test01(int prob, double p)
//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SHEPARD_VALUE_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");
        Console.WriteLine("  using Shepard interpolation with P = " + p + "");

        int dim_num = TestInterp.p00_dim_num(prob);

        int nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xy = TestInterp.p00_data(prob, dim_num, nd);

        switch (p)
        {
            case 0.0:
                typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");
                break;
        }

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

        double[] yi = Shepard.shepard_value_1d(nd, xd, yd, p, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = "
                          + int_error + "");

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
        yi = Shepard.shepard_value_1d(nd, xd, yd, p, ni, xi);

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
        Console.WriteLine("  Normalized length of Shepard interpolant          = " + li + "");

    }
}