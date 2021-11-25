using System;
using Burkardt.PiecewiseLinear;
using Burkardt.Types;
using InterpTest;

namespace PWLInterp2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PWL_INTERP_2D_TEST.
        //
        //  Discussion:
        //
        //    PWL_INTERP_2D_TEST tests the PWL_INTERP_2D library.
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
        int[] n_test = { 2, 3, 4, 5, 9 };
        const int n_test_num = 5;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_2D_TEST:");
        Console.WriteLine("  Test the PWL_INTERP_2D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The test needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();
        //
        //  Numerical tests.
        //
        for (prob = 1; prob <= prob_num; prob++)
        {
            int i;
            for (i = 0; i < n_test_num; i++)
            {
                int n = n_test[i];
                test01(prob, n);
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_INTERP_2D_TEST01 tests PWL_INTERP_2D.
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
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int N, the grid size in each dimension.
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("PWL_INTERP_2D_TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem # " + prob + "");
        Console.WriteLine("  Using polynomial interpolant of product degree " + n + " x " + n + "");

        int nd = n * n;
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xd_1d = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
        double[] yd_1d = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);

        double[] xd = new double[n * n];
        double[] yd = new double[n * n];
        double[] zd = new double[n * n];

        int ij = 0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                xd[ij] = xd_1d[i];
                yd[ij] = yd_1d[j];
                ij += 1;
            }
        }

        Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

        switch (nd)
        {
            case <= 20:
                typeMethods.r8vec3_print(nd, xd, yd, zd, "  X, Y, Z data:");
                break;
        }

        //
        //  #1:  Does interpolant match function at data points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = typeMethods.r8vec_copy_new(ni, yd);

        double[] zi = Interp2D.pwl_interp_2d(n, n, xd_1d, yd_1d, zd, ni, xi, yi);

        switch (ni)
        {
            case <= 20:
                typeMethods.r8vec3_print(ni, xi, yi, zi, "  X, Y, Z interpolation:");
                break;
        }

        double int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  RMS data interpolation error = " + int_error + "");

        switch (nd)
        {
            //
            //  #2:  Does interpolant approximate data at midpoints?
            //
            case > 1:
            {
                double[] xi_1d = new double[n - 1];
                double[] yi_1d = new double[n - 1];

                for (i = 0; i < n - 1; i++)
                {
                    xi_1d[i] = 0.5 * (xd_1d[i] + xd_1d[i + 1]);
                }

                for (i = 0; i < n - 1; i++)
                {
                    yi_1d[i] = 0.5 * (yd_1d[i] + yd_1d[i + 1]);
                }

                ni = (n - 1) * (n - 1);

                xi = new double[ni];
                yi = new double[ni];
                double[] zdm = new double[ni];

                ij = 0;
                for (j = 0; j < n - 1; j++)
                {
                    for (i = 0; i < n - 1; i++)
                    {
                        xi[ij] = xi_1d[i];
                        yi[ij] = yi_1d[j];
                        ij += 1;
                    }
                }

                Data_2D.f00_f0(prob, ni, xi, yi, ref zdm);

                zi = Interp2D.pwl_interp_2d(n, n, xd_1d, yd_1d, zd, ni, xi, yi);

                double app_error = typeMethods.r8vec_norm_affine(ni, zi, zdm) / ni;

                Console.WriteLine("  RMS data approximation error = " + app_error + "");
                break;
            }
        }
    }
}