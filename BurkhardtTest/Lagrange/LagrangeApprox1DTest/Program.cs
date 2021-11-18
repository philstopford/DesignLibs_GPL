using System;
using Burkardt.Lagrange;
using Burkardt.Types;
using InterpTest;

namespace Lagrange1DTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGRANGE_APPROX_1D_TEST.
        //
        //  Discussion:
        //
        //    LAGRANGE_APPROX_1D_TEST tests the LAGRANGE_APPROX_1D library.
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
        int j;
        int k;
        int m;
        int[] m_test =  {
                0, 1, 2, 3, 4, 8, 16
            }
            ;
        int m_test_num = 7;
        int nd;
        int[] nd_test =  {
                16, 64, 1000
            }
            ;
        int nd_test_num = 3;
        int prob;
        int prob_num;

        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_APPROX_1D_TEST:");
        Console.WriteLine("  Test the LAGRANGE_APPROX_1D library.");
        Console.WriteLine("  The R8LIB library is needed.");
        Console.WriteLine("  The QR_SOLVE library is needed.");
        Console.WriteLine("  These tests need the TEST_INTERP_1D library.");

        prob_num = Data_1D.p00_prob_num();

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (j = 0; j < m_test_num; j++)
            {
                m = m_test[j];
                for (k = 0; k < nd_test_num; k++)
                {
                    nd = nd_test[k];
                    test02(prob, m, nd);
                }
            }
        }

        for (prob = 1; prob <= prob_num; prob++)
        {
            for (j = 0; j < m_test_num; j++)
            {
                m = m_test[j];
                for (k = 0; k < nd_test_num; k++)
                {
                    nd = nd_test[k];
                    test03(prob, m, nd);
                }
            }
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGRANGE_APPROX_1D_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test02(int prob, int m, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests LAGRANGE_APPROX_1D with evenly spaced data
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
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int M, the polynomial approximant degree.
        //
        //    Input, int ND, the number of data points.
        //
    {
        double a;
        double b;
        double int_error;
        int ni;
        double[] xd;
        double[] xi;
        double[] yd;
        double[] yi;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Approximate evenly spaced data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Use polynomial approximant of degree " + m + "");
        Console.WriteLine("  Number of data points = " + nd + "");

        a = 0.0;
        b = 1.0;
        xd = typeMethods.r8vec_linspace_new(nd, a, b);

        yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does approximant come close to function at data points?
        //
        ni = nd;
        xi = typeMethods.r8vec_copy_new(ni, xd);
        yi = Lagrange1D.lagrange_approx_1d(m, nd, xd, yd, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per data node = " + int_error + "");
    }

    private static void test03(int prob, int m, int nd)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests LAGRANGE_APPROX_1D with Chebyshev spaced data.
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
        //  Parameters:
        //
        //    Input, int PROB, the problem index.
        //
        //    Input, int M, the polynomial approximant degree.
        //
        //    Input, int ND, the number of data points.
        //
    {
        double a;
        double b;
        double int_error;
        int ni;
        double[] xd;
        double[] xi;
        double[] yd;
        double[] yi;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  Approximate Chebyshev-spaced data from TEST_INTERP_1D problem #" + prob + "");
        Console.WriteLine("  Use polynomial approximant of degree " + m + "");
        Console.WriteLine("  Number of data points = " + nd + "");

        a = 0.0;
        b = 1.0;
        xd = typeMethods.r8vec_cheby_extreme_new(nd, a, b);

        yd = Data_1D.p00_f(prob, nd, xd);

        switch (nd)
        {
            case < 10:
                typeMethods.r8vec2_print(nd, xd, yd, "  Data array:");
                break;
        }

        //
        //  #1:  Does interpolant match function at interpolation points?
        //
        ni = nd;
        xi = typeMethods.r8vec_copy_new(ni, xd);
        yi = Lagrange1D.lagrange_approx_1d(m, nd, xd, yd, ni, xi);

        int_error = typeMethods.r8vec_norm_affine(nd, yi, yd) / ni;

        Console.WriteLine("");
        Console.WriteLine("  L2 approximation error averaged per data node = " + int_error + "");
    }
}