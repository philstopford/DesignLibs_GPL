using System;
using Burkardt.Lagrange;
using Burkardt.Types;
using InterpTest;

namespace LagrangeInterp2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGRANGE_INTERP_2D_TEST.
        //
        //  Discussion:
        //
        //    LAGRANGE_INTERP_2D_TEST tests the LAGRANGE_INTERP_2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] m_test =  {
                1, 2, 3, 4, 8
            }
            ;
        int m_test_num = 5;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_2D_TEST:");
        Console.WriteLine("  Test the LAGRANGE_INTERP_2D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();
        //
        //  Numerical tests.
        //
        for (prob = 1; prob <= prob_num; prob++)
        {
            int i;
            for (i = 0; i < m_test_num; i++)
            {
                int m = m_test[i];
                test01(prob, m);
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_2D_TEST01 tests LAGRANGE_INTERP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int M, the polynomial degree in each dimension.
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_INTERP_2D_TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
        Console.WriteLine("  Using polynomial interpolant of product degree " + m + " x " + m + "");

        int nd = (m + 1) * (m + 1);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xd_1d = typeMethods.r8vec_cheby_extreme_new(m + 1, 0.0, 1.0);
        double[] yd_1d = typeMethods.r8vec_cheby_extreme_new(m + 1, 0.0, 1.0);

        double[] xd = new double[(m + 1) * (m + 1)];
        double[] yd = new double[(m + 1) * (m + 1)];
        double[] zd = new double[(m + 1) * (m + 1)];

        int ij = 0;
        for (j = 0; j < m + 1; j++)
        {
            for (i = 0; i < m + 1; i++)
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

        double[] xi = new double[ni];
        double[] yi = new double[ni];

        for (i = 0; i < ni; i++)
        {
            xi[i] = xd[i];
            yi[i] = yd[i];
        }

        double[] zi = Lagrange2D.lagrange_interp_2d(m, m, xd_1d, yd_1d, zd, ni, xi, yi);

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
                double[] xi_1d = new double[m];
                double[] yi_1d = new double[m];

                for (i = 0; i < m; i++)
                {
                    xi_1d[i] = 0.5 * (xd_1d[i] + xd_1d[i + 1]);
                }

                for (i = 0; i < m; i++)
                {
                    yi_1d[i] = 0.5 * (yd_1d[i] + yd_1d[i + 1]);
                }

                ni = m * m;

                xi = new double[ni];
                yi = new double[ni];
                double[] zdm = new double[ni];

                ij = 0;
                for (j = 0; j < m; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        xi[ij] = xi_1d[i];
                        yi[ij] = yi_1d[j];
                        ij += 1;
                    }
                }

                Data_2D.f00_f0(prob, ni, xi, yi, ref zdm);

                zi = Lagrange2D.lagrange_interp_2d(m, m, xd_1d, yd_1d, zd, ni, xi, yi);

                double app_error = typeMethods.r8vec_norm_affine(ni, zi, zdm) / ni;

                Console.WriteLine("");
                Console.WriteLine("  RMS data approximation error = " + app_error + "");
                break;
            }
        }
    }
}