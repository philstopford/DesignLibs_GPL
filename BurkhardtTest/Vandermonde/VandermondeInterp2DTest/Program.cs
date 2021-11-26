﻿using System;
using Burkardt.PolynomialNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using Burkardt.Uniform;
using InterpTest;
using Vandermonde = Burkardt.Interpolation.Vandermonde;

namespace VandermondeInterp2D;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for VANDERMONDE_INTERP_1D_TEST.
        //
        //  Discussion:
        //
        //    VANDERMONDE_INTERP_2D_TEST tests the VANDERMONDE_INTERP_2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] m_test = {1, 2, 3, 4, 8};
        const int m_test_num = 5;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_INTERP_2D_TEST:");
        Console.WriteLine("  Test the VANDERMONDE_INTERP_2D library.");
        Console.WriteLine("  The QR_SOLVE library is needed.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();
        for (prob = 1; prob <= prob_num; prob++)
        {
            int j;
            for (j = 0; j < m_test_num; j++)
            {
                int m = m_test[j];
                test01(prob, m);
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_INTERP_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests VANDERMONDE_INTERP_2D_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int M, the degree of interpolation.
        //
    {
        const bool debug = false;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
        Console.WriteLine("  Create an interpolant of total degree " + m + "");
        int tmp1 = typeMethods.triangle_num(m + 1);
        Console.WriteLine("  Number of data values needed is " + tmp1 + "");

        int seed = 123456789;

        double[] xd = UniformRNG.r8vec_uniform_01_new(tmp1, ref seed);
        double[] yd = UniformRNG.r8vec_uniform_01_new(tmp1, ref seed);
        double[] zd = new double[tmp1];
        Data_2D.f00_f0(prob, tmp1, xd, yd, ref zd);

        switch (debug)
        {
            case true:
                typeMethods.r8vec3_print(tmp1, xd, yd, zd, "  X, Y, Z data:");
                break;
        }

        //
        //  Compute the Vandermonde matrix.
        //
        double[] a = Vandermonde.vandermonde_interp_2d_matrix(tmp1, m, xd, yd);
        //
        //  Solve linear system.
        //
        double[] c = QRSolve.qr_solve(tmp1, tmp1, a, zd);
        //
        //  #1:  Does interpolant match function at data points?
        //
        int ni = tmp1;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = typeMethods.r8vec_copy_new(ni, yd);
        double[] zi = Polynomial.r8poly_values_2d(m, c, ni, xi, yi);

        double app_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 data interpolation error = " + app_error + "");
    }
}