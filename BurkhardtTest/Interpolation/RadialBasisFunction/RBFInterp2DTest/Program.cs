using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace RBFInterp2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for RBF_INTERP_2D_TEST.
        //
        //  Discussion:
        //
        //    RBF_INTERP_2D_TEST tests the RBF_INTERP_2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int prob;

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_2D_TEST:");
        Console.WriteLine("  Test the RBF_INTERP_2D library.");
        Console.WriteLine("  The R8LIB library is required.");
        Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();
        const int g = 1;

        for (prob = 1; prob <= prob_num; prob++)
        {
            test01(prob, g, RadialBasisFunctions.phi1, "phi1");
            test01(prob, g, RadialBasisFunctions.phi2, "phi2");
            test01(prob, g, RadialBasisFunctions.phi3, "phi3");
            test01(prob, g, RadialBasisFunctions.phi4, "phi4");
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int g,
            Func<int, double[], double, double[], double[]> phi, string phi_name)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RBF_INTERP_2D_TEST01 tests RBF_INTERP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the index of the problem.
        //
        //    Input, int G, the index of the grid.
        //
        //    Input, void PHI ( int n, double r[], double r0, double v[] ), the 
        //    radial basis function.
        //
        //    Input, string PHI_NAME, the name of the radial basis function.
        //
    {
        const bool debug = false;
        int i;

        Console.WriteLine("");
        Console.WriteLine("RBF_INTERP_2D_TEST01:");
        Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
        Console.WriteLine("  using grid #" + g + "");
        Console.WriteLine("  using radial basis function \"" + phi_name + "\".");

        int nd = Data_2D.g00_size(g);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xd = new double[nd];
        double[] yd = new double[nd];
        Data_2D.g00_xy(g, nd, ref xd, ref yd);

        double[] zd = new double[nd];
        Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

        switch (debug)
        {
            case true:
                typeMethods.r8vec3_print(nd, xd, yd, zd, "  X, Y, Z data:");
                break;
        }

        const int m = 2;
        double[] xyd = new double[2 * nd];

        for (i = 0; i < nd; i++)
        {
            xyd[0 + i * 2] = xd[i];
            xyd[1 + i * 2] = yd[i];
        }

        double xmax = typeMethods.r8vec_max(nd, xd);
        double xmin = typeMethods.r8vec_min(nd, xd);
        double ymax = typeMethods.r8vec_max(nd, yd);
        double ymin = typeMethods.r8vec_min(nd, yd);
        double volume = (xmax - xmin) * (ymax - ymin);

        const double e = 1.0 / m;
        double r0 = Math.Pow(volume / nd, e);

        Console.WriteLine("  Setting R0 = " + r0 + "");

        double[] w = RadialBasisFunctions.rbf_weight(m, nd, xyd, r0, phi, zd);
        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        double[] xyi = typeMethods.r8mat_copy_new(2, nd, xyd);

        double[] zi = RadialBasisFunctions.rbf_interp(m, nd, xyd, r0, phi, w, nd, xyi);

        double int_error = typeMethods.r8vec_norm_affine(nd, zi, zd) / nd;

        Console.WriteLine("");
        Console.WriteLine("  L2 interpolation error averaged per interpolant node = "
                          + int_error + "");

    }
}