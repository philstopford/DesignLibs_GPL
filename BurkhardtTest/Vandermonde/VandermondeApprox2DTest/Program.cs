using System;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.SolveNS;
using Burkardt.Types;
using InterpTest;

namespace VandermondeApprox2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for VANDERMONDE_APPROX_2D_TEST.
        //
        //  Discussion:
        //
        //    VANDERMONDE_APPROX_2D_TEST tests the VANDERMONDE_APPROX_2D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] m_test = {0, 1, 2, 4, 8};
        const int m_test_num = 5;
        int prob;

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_APPROX_2D_TEST:");
        Console.WriteLine("  Test the VANDERMONDE_APPROX_2D library.");
        Console.WriteLine("  The QR_SOLVE library is needed.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

        int prob_num = Data_2D.f00_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            int grid = 1;
            int j;
            for (j = 0; j < m_test_num; j++)
            {
                int m = m_test[j];
                test01(prob, grid, m);
            }
        }

        Console.WriteLine("");
        Console.WriteLine("VANDERMONDE_APPROX_2D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int prob, int grd, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_APPROX_2D_TEST01 tests VANDERMONDE_APPROX_2D_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int PROB, the problem number.
        //
        //    Input, int GRD, the grid number.
        //    (Can't use GRID as the name because that's also a plotting function.)
        //
        //    Input, int M, the total polynomial degree.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Approximate data from TEST_INTERP_2D problem #" + prob + "");
        Console.WriteLine("  Use grid from TEST_INTERP_2D with index #" + grd + "");
        Console.WriteLine("  Using polynomial approximant of total degree " + m + "");

        int nd = Data_2D.g00_size(grd);
        Console.WriteLine("  Number of data points = " + nd + "");

        double[] xd = new double[nd];
        double[] yd = new double[nd];
        Data_2D.g00_xy(grd, nd, ref xd, ref yd);

        double[] zd = new double[nd];
        Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec3_print(nd, xd, yd, zd, "  X, Y, Z data:");
                break;
        }

        //
        //  Compute the Vandermonde matrix.
        //
        int tm = typeMethods.triangle_num(m + 1);
        double[] a = VandermondeMatrix.vandermonde_approx_2d_matrix(nd, m, tm, xd, yd);
        //
        //  Solve linear system.
        //
        double[] c = QRSolve.qr_solve(nd, tm, a, zd);
        //
        //  #1:  Does approximant match function at data points?
        //
        int ni = nd;
        double[] xi = typeMethods.r8vec_copy_new(ni, xd);
        double[] yi = typeMethods.r8vec_copy_new(ni, yd);
        double[] zi = Polynomial.r8poly_values_2d(m, c, ni, xi, yi);

        double app_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 data approximation error = " + app_error + "");

    }
}