﻿using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using Burkardt.Uniform;
using InterpTest;

namespace ShepardnDTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SHEPARD_INTERP_ND_TEST.
        //
        //  Discussion:
        //
        //    SHEPARD_INTERP_ND_TEST tests the SHEPARD_INTERP_ND library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int m;
        double p;
        double[] p_test = {1.0, 2.0, 4.0, 8.0};
        const int p_test_num = 4;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_ND_TEST:");
        Console.WriteLine("  Test the SHEPARD_INTERP_ND library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test also needs the TEST_INTERP_ND library.");
        //
        //  Look at Shepard interpolant on an irregular grid.
        //
        int nd = 25;

        int prob_num = Data_nD.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (m = 2; m <= 5; m += 3)
            {
                for (j = 0; j < p_test_num; j++)
                {
                    p = p_test[j];
                    test01(prob, p, m, nd);
                }

            }
        }

        //
        //  Look at Shepard interpolant on a regular N1D^M grid.
        //
        int n1d = 5;

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (m = 2; m <= 5; m += 3)
            {
                for (j = 0; j < p_test_num; j++)
                {
                    p = p_test[j];
                    test02(prob, p, m, n1d);
                }
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SHEPARD_INTERP_ND_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01(int prob, double p, int m, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SHEPARD_INTERP on an irregular grid.
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
        //    Input, double P, the power used in the distance weighting.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int ND, the number of data points.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_ND problem #" + prob + "");
        Console.WriteLine("  using Shepard interpolation with P = " + p + "");
        Console.WriteLine("  spatial dimension M = " + m + "");
        Console.WriteLine("  and an irregular grid of ND = " + nd + " data points.");
        //
        //  Set problem parameters:
        //
        int seed = 123456789;
        double[] c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        double[] w = UniformRNG.r8vec_uniform_01_new(m, ref seed);

        double[] xd = UniformRNG.r8mat_uniform_01_new(m, nd, ref seed);

        double[] zd = Data_nD.p00_f(prob, m, c, w, nd, xd);
        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8mat_copy_new(m, ni, xd);
        double[] zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        ni = 50;
        xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
        zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);
        double[] ze = Data_nD.p00_f(prob, m, c, w, ni, xi);

        double app_error = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

        Console.WriteLine("  L2 approximation error averaged per 1000 samples =     " + app_error + "");

    }

    private static void test02(int prob, double p, int m, int n1d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SHEPARD_INTERP_ND on a regular N1D^M grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, double P, the power used in the distance weighting.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N1D, the number of points in 1D.
        //
    {
        int i;
        //
        //  Set problem parameters:
        //
        int seed = 123456789;
        double[] c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        double[] w = UniformRNG.r8vec_uniform_01_new(m, ref seed);

        int nd = (int) Math.Pow(n1d, m);

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_ND problem #" + prob + "");
        Console.WriteLine("  using Shepard interpolation with P = " + p + "");
        Console.WriteLine("  spatial dimension M = " + m + "");
        Console.WriteLine("  and a regular grid of N1D^M = " + nd + " data points.");

        const double a = 0.0;
        const double b = 1.0;

        double[] x1d = typeMethods.r8vec_linspace_new(n1d, a, b);

        double[] xd = new double[m * nd];
        for (i = 0; i < m; i++)
        {
            typeMethods.r8vecDPData data = new();
            typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
        }

        double[] zd = Data_nD.p00_f(prob, m, c, w, nd, xd);
        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8mat_copy_new(m, nd, xd);
        double[] zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
        //
        ni = 1000;
        xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);

        zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

        double[] ze = Data_nD.p00_f(prob, m, c, w, ni, xi);

        double app_error = typeMethods.r8vec_norm_affine(ni, zi, ze) / ni;

        Console.WriteLine("  L2 approximation error averaged per 1000 samples =     " + app_error + "");

    }
}