﻿using System;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using InterpTest;

namespace VandermondeApprox1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for VANDERMONDE_APPROX_1D_TEST.
        //
        //  Discussion:
        //
        //    VANDERMONDE_APPROX_1D_TEST tests the VANDERMONDE_APPROX_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] m_test = {0, 1, 2, 3, 4, 5, 9, 12};
        const int m_test_num = 8;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_APPROX_1D_TEST:");
        Console.WriteLine("  Test the VANDERMONDE_APPROX_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The QR_SOLVE library is needed.");
        Console.WriteLine("  The test needs the TEST_INTERP libary.");

        int prob_num = Data_1D.p00_prob_num();
        for (prob = 1; prob <= prob_num; prob++)
        {
            int j;
            for (j = 0; j < m_test_num; j++)
            {
                int m = m_test[j];
                test01(prob, m);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_APPROX_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests VANDERMONDE_APPROX_1D_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int M, the polynomial degree.
        //
    {
        const bool debug = false;
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Approximate data from TEST_INTERP problem #" + prob + "");

        int nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xy = TestInterp.p00_data(prob, 2, nd);

        switch (debug)
        {
            case true:
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
        //  Compute the Vandermonde matrix.
        //
        Console.WriteLine("  Using polynomial approximant of degree " + m + "");

        double[] a = VandermondeMatrix.vandermonde_approx_1d_matrix(nd, m, xd);
        //
        //  Solve linear system.
        //
        double[] c = QRSolve.qr_solve(nd, m + 1, a, yd);
        //
        //  #1:  Does approximant match function at data points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = Polynomial.r8poly_values(m, c, ni, xi);

        double app_error = typeMethods.r8vec_norm_affine(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 data approximation error = " + app_error + "");

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

        yi = Polynomial.r8poly_values(m, c, ni, xi);

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
        Console.WriteLine("  Normalized length of polynomial interpolant       = " + li + "");
    }
}