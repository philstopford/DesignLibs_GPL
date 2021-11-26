using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace BarycentricInterp1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BARYCENTRIC_INTERP_1D_TEST.
        //
        //  Discussion:
        //
        //    BARYCENTRIC_INTERP_1D_TEST tests the BARYCENTRIC_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int nd;
        int[] nd_test =  {
                4, 8, 16, 32, 64, 1000
            }
            ;
        const int nd_test_num = 6;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("BARYCENTRIC_INTERP_1D_TEST:");
        Console.WriteLine("  Test the BARYCENTRIC_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The tests need the TEST_INTERP_1D library.");

        int prob_num = Data_1D.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (i = 0; i < nd_test_num; i++)
            {
                nd = nd_test[i];
                lagcheby1_interp_1d_test(prob, nd);
            }
        }

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (i = 0; i < nd_test_num; i++)
            {
                nd = nd_test[i];
                lagcheby2_interp_1d_test(prob, nd);
            }
        }

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (i = 0; i < nd_test_num; i++)
            {
                nd = nd_test[i];
                lageven_interp_1d_test(prob, nd);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("BARYCENTRIC_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void lagcheby1_interp_1d_test(int prob, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGCHEBY1_INTERP_1D_TEST tests LAGCHEBY1_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
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
        Console.WriteLine("");
        Console.WriteLine("LAGCHEBY1_INTERP_1D_TEST:");
        Console.WriteLine("  LAGCHEBY1_INTERP_1D uses Chebyshev Type 1 spacing for data points.");
        Console.WriteLine("  Interpolate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Number of data points = " + nd + "");
        //
        //  Define the data.
        //
        double a = 0.0;
        double b = +1.0;
        double[] xd = typeMethods.r8vec_cheby1space_new(nd, a, b);
        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does the interpolant match the function at the interpolation points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = BarycentricInterp1D.lagcheby1_interp_1d(nd, xd, yd, ni, xi);

        double int_error = typeMethods.r8vec_norm_affine(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

    }

    private static void lagcheby2_interp_1d_test(int prob, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGCHEBY2_INTERP_1D_TEST tests LAGCHEBY2_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
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
        Console.WriteLine("");
        Console.WriteLine("LAGCHEBY2_INTERP_1D_TEST:");
        Console.WriteLine("  LAGCHEBY2_INTERP_1D uses Chebyshev Type 2 spacing for data points.");
        Console.WriteLine("  Interpolate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Number of data points = " + nd + "");
        //
        //  Define the data.
        //
        const double a = 0.0;
        const double b = +1.0;
        double[] xd = typeMethods.r8vec_cheby2space_new(nd, a, b);
        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does the interpolant match the function at the interpolation points?
        //
        double[] xi = typeMethods.r8vec_copy_new(nd, xd);
        double[] yi = BarycentricInterp1D.lagcheby2_interp_1d(nd, xd, yd, nd, xi);

        double int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

    }

    private static void lageven_interp_1d_test(int prob, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGEVEN_INTERP_1D_TEST tests LAGEVEN_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
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
        Console.WriteLine("");
        Console.WriteLine("LAGEVEN_INTERP_1D_TEST:");
        Console.WriteLine("  LAGEVEN_INTERP_1D uses even spacing for data points.");
        Console.WriteLine("  Interpolate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Number of data points = " + nd + "");
        //
        //  Define the data.
        //
        const double a = 0.0;
        const double b = +1.0;
        double[] xd = typeMethods.r8vec_midspace_new(nd, a, b);
        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does the interpolant match the function at the interpolation points?
        //
        double[] xi = typeMethods.r8vec_copy_new(nd, xd);
        double[] yi = BarycentricInterp1D.lageven_interp_1d(nd, xd, yd, nd, xi);

        double int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

    }
}