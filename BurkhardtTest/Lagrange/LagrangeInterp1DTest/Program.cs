using System;
using Burkardt.Lagrange;
using Burkardt.Types;
using InterpTest;

namespace LagrangeInterp1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGRANGE_INTERP_1D_TEST.
        //
        //  Discussion:
        //
        //    LAGRANGE_INTERP_1D_TEST tests the LAGRANGE_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nd_test_num = 6;

        int j;
        int nd;
        int[] nd_test =  {
                4, 8, 16, 32, 64, 256
            }
            ;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_1D_TEST:");
        Console.WriteLine("  Test the LAGRANGE_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  These tests need the TEST_INTERP_1D library.");

        int prob_num = Data_1D.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (j = 0; j < nd_test_num; j++)
            {
                nd = nd_test[j];
                test01(prob, nd);
            }
        }

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (j = 0; j < nd_test_num; j++)
            {
                nd = nd_test[j];
                test02(prob, nd);
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int nd)

        //*****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LAGRANGE_VALUE_1D with evenly spaced data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int ND, the number of data points to use.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Use even spacing for data points.");
        Console.WriteLine("  Number of data points = " + nd + "");

        const double a = 0.0;
        const double b = 1.0;

        double[] xd = typeMethods.r8vec_linspace_new(nd, a, b);

        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = Lagrange1D.lagrange_value_1d(nd, xd, yd, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
        //  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
        //  (YMAX-YMIN).
        //
        double ymin = typeMethods.r8vec_min(nd, yd);
        double ymax = typeMethods.r8vec_max(nd, yd);

        ni = 501;
        xi = typeMethods.r8vec_linspace_new(ni, a, b);
        yi = Lagrange1D.lagrange_value_1d(nd, xd, yd, ni, xi);

        double ld = 0.0;
        for (i = 0; i < nd - 1; i++)
        {
            ld += Math.Sqrt(Math.Pow((xd[i + 1] - xd[i]) / (b - a), 2)
                            + Math.Pow((yd[i + 1] - yd[i]) / (ymax - ymin), 2));
        }

        double li = 0.0;
        for (i = 0; i < ni - 1; i++)
        {
            li += Math.Sqrt(Math.Pow((xi[i + 1] - xi[i]) / (b - a), 2)
                            + Math.Pow((yi[i + 1] - yi[i]) / (ymax - ymin), 2));
        }

        Console.WriteLine("");
        Console.WriteLine("  Normalized length of piecewise linear interpolant = " + ld + "");
        Console.WriteLine("  Normalized length of polynomial interpolant       = " + li + "");
    }

    private static void test02(int prob, int nd)

        //*****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests LAGRANGE_VALUE_1D with Chebyshev spaced data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int ND, the number of data points to use.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Use Chebyshev spacing for data points.");
        Console.WriteLine("  Number of data points = " + nd + "");

        const double a = 0.0;
        const double b = 1.0;

        double[] xd = typeMethods.r8vec_cheby_extreme_new(nd, a, b);

        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = Lagrange1D.lagrange_value_1d(nd, xd, yd, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
        //  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
        //  (YMAX-YMIN).
        //
        double ymin = typeMethods.r8vec_min(nd, yd);
        double ymax = typeMethods.r8vec_max(nd, yd);

        ni = 501;
        xi = typeMethods.r8vec_linspace_new(ni, a, b);
        yi = Lagrange1D.lagrange_value_1d(nd, xd, yd, ni, xi);

        double ld = 0.0;
        for (i = 0; i < nd - 1; i++)
        {
            ld += Math.Sqrt(Math.Pow((xd[i + 1] - xd[i]) / (b - a), 2)
                            + Math.Pow((yd[i + 1] - yd[i]) / (ymax - ymin), 2));
        }

        double li = 0.0;
        for (i = 0; i < ni - 1; i++)
        {
            li += Math.Sqrt(Math.Pow((xi[i + 1] - xi[i]) / (b - a), 2)
                            + Math.Pow((yi[i + 1] - yi[i]) / (ymax - ymin), 2));
        }

        Console.WriteLine("");
        Console.WriteLine("  Normalized length of piecewise linear interpolant = " + ld + "");
        Console.WriteLine("  Normalized length of polynomial interpolant       = " + li + "");
    }
}