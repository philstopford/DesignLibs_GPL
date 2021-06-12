using System;
using Burkardt.ShepardInterpolation;
using Burkardt.Types;
using InterpTest;

namespace Shepard2DTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SHEPARD_INTERP_2D_TEST.
//
//  Discussion:
//
//    SHEPARD_INTERP_2D_TEST tests the SHEPARD_INTERP_2D library.
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
            int g;
            int j;
            double p;
            double[] p_test =  {
                1.0, 2.0, 4.0, 8.0
            }
            ;
            int p_test_num = 4;
            int prob;
            int prob_num;

            Console.WriteLine("");
            Console.WriteLine("SHEPARD_INTERP_2D_TEST:");
            Console.WriteLine("  Test the SHEPARD_INTERP_2D library.");
            Console.WriteLine("  The R8LIB library is needed.");
            Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

            prob_num = Data_2D.f00_num();
            g = 1;

            for (prob = 1; prob <= prob_num; prob++)
            {
                for (j = 0; j < p_test_num; j++)
                {
                    p = p_test[j];
                    test01(prob, g, p);
                }
            }

            Console.WriteLine("");
            Console.WriteLine("SHEPARD_INTERP_2D_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int prob, int g, double p)

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests SHEPARD_INTERP_2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int G, the grid number.
//
//    Input, double P, the power used in the distance weighting.
//
        {
            bool debug = false;
            double int_error;
            int nd;
            int ni;
            double[] xd;
            double[] xi;
            double[] yd;
            double[] yi;
            double[] zd;
            double[] zi;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
            Console.WriteLine("  using grid #" + g + "");
            Console.WriteLine("  using Shepard interpolation with P = " + p + "");

            nd = Data_2D.g00_size(g);
            Console.WriteLine("  Number of data points = " + nd + "");

            xd = new double[nd];
            yd = new double[nd];
            Data_2D.g00_xy(g, nd, ref xd, ref yd);

            zd = new double[nd];
            Data_2D.f00_f0(prob, nd, xd, yd, ref zd);

            if (debug)
            {
                typeMethods.r8vec3_print(nd, xd, yd, zd, "  X, Y, Z data:");
            }

//
//  #1:  Does interpolant match function at interpolation points?
//
            ni = nd;
            xi = typeMethods.r8vec_copy_new(ni, xd);
            yi = typeMethods.r8vec_copy_new(ni, yd);

            zi = Shepard.shepard_interp_2d(nd, xd, yd, zd, p, ni, xi, yi);

            int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / (double) (ni);

            Console.WriteLine("");
            Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

        }
    }
}