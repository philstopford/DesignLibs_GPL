using System;
using Burkardt.PiecewiseLinear;
using Burkardt.Probability;
using Burkardt.Types;
using InterpTest;

namespace PWLApprox1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PWL_APPROX_1D_TEST.
        //
        //  Discussion:
        //
        //    PWL_APPROX_1D_TEST tests the PWL_APPROX_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] nc_test = { 2, 4, 8, 16 };
        const int nc_test_num = 4;
        int[] nd_test = { 16, 64 };
        const int nd_test_num = 2;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("PWL_APPROX_1D_TEST:");
        Console.WriteLine("  Test the PWL_APPROX_1D library.");
        Console.WriteLine("  The QR_SOLVE library is needed.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The test also needs the TEST_INTERP_1D library.");

        int prob_num = ProbabilityFunctions.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            int j;
            for (j = 0; j < nc_test_num; j++)
            {
                int nc = nc_test[j];
                int k;
                for (k = 0; k < nd_test_num; k++)
                {
                    int nd = nd_test[k];
                    test01(prob, nc, nd);
                }
            }
        }

        Console.WriteLine("");
        Console.WriteLine("PWL_APPROX_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int nc, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PWL_APPROX_1D.
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
        //    Input, int PROB, the problem index.
        //
        //    Input, int NC, the number of control points.
        //
        //    Input, int ND, the number of data points.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Approximate data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Number of control points = " + nc + "");
        Console.WriteLine("  Number of data points = " + nd + "");

        const double xmin = 0.0;
        const double xmax = 1.0;

        double[] xd = typeMethods.r8vec_linspace_new(nd, xmin, xmax);
        double[] yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  Determine control values.
        // 
        double[] xc = typeMethods.r8vec_linspace_new(nc, xmin, xmax);
        double[] yc = Approx1D.pwl_approx_1d(nd, xd, yd, nc, xc);
        //
        //  #1:  Does approximant come close to function at data points?
        //
        double[] xi = typeMethods.r8vec_copy_new(nd, xd);
        double[] yi = Approx1D.pwl_interp_1d(nc, xc, yc, nd, xi);

        double app_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per data node = " + app_error + "");
    }
}