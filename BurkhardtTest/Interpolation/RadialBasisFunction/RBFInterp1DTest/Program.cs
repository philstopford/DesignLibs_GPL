using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace RBFInterp1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for RBF_INTERP_1D_TEST.
        //
        //  Discussion:
        //
        //    RBF_INTERP_1D_TEST tests the RBF_INTERP_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int nd;
        int prob;
        int prob_num;
        double r0;
        double[] xd;
        double xmax;
        double xmin;
        double[] xy;

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_1D_TEST:");
        Console.WriteLine("  Test the RBF_INTERP_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The test needs the TEST_INTERP library.");

        prob_num = TestInterp.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            //
            //  Determine an appropriate value of R0, the spacing parameter.
            //
            nd = TestInterp.p00_data_num(prob);
            xy = TestInterp.p00_data(prob, 2, nd);
            xd = new double[nd];
            for (i = 0; i < nd; i++)
            {
                xd[i] = xy[0 + i * 2];
            }

            xmax = typeMethods.r8vec_max(nd, xd);
            xmin = typeMethods.r8vec_min(nd, xd);
            r0 = (xmax - xmin) / (nd - 1);

            test01(prob, RadialBasisFunctions.phi1, "phi1", r0);
            test01(prob, RadialBasisFunctions.phi2, "phi2", r0);
            test01(prob, RadialBasisFunctions.phi3, "phi3", r0);
            test01(prob, RadialBasisFunctions.phi4, "phi4", r0);
        }

        /*
        Terminate.
        */
        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, Func<int, double[], double, double[], double[]> phi,
            string phi_name, double r0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests RBF_INTERP_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the index of the problem.
        //
        //    Input, double PHI ( int n, double r[], double r0, double v[] ), 
        //    the name of the radial basis function.
        //
        //    Input, string PHI_NAME, the name of the radial basis function.
        //
        //    Input, double R0, the scale factor.  Typically, this might be
        //    a small multiple of the average distance between points.
        //
    {
        bool debug = false;
        int i;
        double int_error;
        double ld;
        double li;
        int m;
        int nd;
        int ni;
        double[] w;
        double[] xd;
        double[] xi;
        double xmax;
        double xmin;
        double[] xy;
        double[] yd;
        double[] yi;
        double ymax;
        double ymin;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP problem #" + prob + "");
        Console.WriteLine("  using radial basis function \"" + phi_name + "\".");
        Console.WriteLine("  Scale factor R0 = " + r0 + "");

        nd = TestInterp.p00_data_num(prob);
        Console.WriteLine("  Number of data points = " + nd + "");

        xy = TestInterp.p00_data(prob, 2, nd);

        switch (debug)
        {
            case true:
                typeMethods.r8mat_transpose_print(2, nd, xy, "  Data array:");
                break;
        }

        xd = new double[nd];
        yd = new double[nd];
        for (i = 0; i < nd; i++)
        {
            xd[i] = xy[0 + i * 2];
            yd[i] = xy[1 + i * 2];
        }

        m = 1;
        w = RadialBasisFunctions.rbf_weight(m, nd, xd, r0, phi, yd);
        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8vec_copy_new(ni, xd);
        yi = RadialBasisFunctions.rbf_interp(m, nd, xd, r0, phi, w, ni, xi);
        int_error = typeMethods.r8vec_norm_affine(ni, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        //
        //  #2: Compare estimated curve length to piecewise linear (minimal) curve length.
        //  Assume data is sorted, and normalize X and Y dimensions by (XMAX-XMIN) and
        //  (YMAX-YMIN).
        //
        xmax = typeMethods.r8vec_max(nd, xd);
        xmin = typeMethods.r8vec_min(nd, xd);
        ymax = typeMethods.r8vec_max(nd, yd);
        ymin = typeMethods.r8vec_min(nd, yd);

        ni = 501;
        xi = typeMethods.r8vec_linspace_new(ni, xmin, xmax);
        yi = RadialBasisFunctions.rbf_interp(m, nd, xd, r0, phi, w, ni, xi);

        ld = 0.0;
        for (i = 0; i < nd - 1; i++)
        {
            ld += Math.Sqrt(Math.Pow((xd[i + 1] - xd[i]) / (xmax - xmin), 2)
                            + Math.Pow((yd[i + 1] - yd[i]) / (ymax - ymin), 2));
        }

        li = 0.0;
        for (i = 0; i < ni - 1; i++)
        {
            li += Math.Sqrt(Math.Pow((xi[i + 1] - xi[i]) / (xmax - xmin), 2)
                            + Math.Pow((yi[i + 1] - yi[i]) / (ymax - ymin), 2));
        }

        Console.WriteLine("  Normalized length of piecewise linear interpolant = " + ld + "");
        Console.WriteLine("  Normalized length of polynomial interpolant       = " + li + "");

    }
}