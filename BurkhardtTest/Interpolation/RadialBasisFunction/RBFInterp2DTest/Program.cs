using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using InterpTest;

namespace RBFInterp2DTest
{
    class Program
    {
        static void Main(string[] args)
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
            int g;
            int prob;
            int prob_num;

            Console.WriteLine("");
            Console.WriteLine("RBF_INTERP_2D_TEST:");
            Console.WriteLine("  Test the RBF_INTERP_2D library.");
            Console.WriteLine("  The R8LIB library is required.");
            Console.WriteLine("  This test also needs the TEST_INTERP_2D library.");

            prob_num = Data_2D.f00_num();
            g = 1;

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

        static void test01(int prob, int g,
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
            bool debug = false;
            double e;
            int i;
            double int_error;
            int m;
            int nd;
            int ni;
            double r0;
            double volume;
            double[] w;
            double[] xd;
            double xmax;
            double xmin;
            double[] xyd;
            double[] xyi;
            double[] yd;
            double ymax;
            double ymin;
            double[] zd;
            double[] zi;

            Console.WriteLine("");
            Console.WriteLine("RBF_INTERP_2D_TEST01:");
            Console.WriteLine("  Interpolate data from TEST_INTERP_2D problem #" + prob + "");
            Console.WriteLine("  using grid #" + g + "");
            Console.WriteLine("  using radial basis function \"" + phi_name + "\".");

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

            m = 2;
            xyd = new double[2 * nd];

            for (i = 0; i < nd; i++)
            {
                xyd[0 + i * 2] = xd[i];
                xyd[1 + i * 2] = yd[i];
            }

            xmax = typeMethods.r8vec_max(nd, xd);
            xmin = typeMethods.r8vec_min(nd, xd);
            ymax = typeMethods.r8vec_max(nd, yd);
            ymin = typeMethods.r8vec_min(nd, yd);
            volume = (xmax - xmin) * (ymax - ymin);

            e = 1.0 / (double) (m);
            r0 = Math.Pow(volume / nd, e);

            Console.WriteLine("  Setting R0 = " + r0 + "");

            w = RadialBasisFunctions.rbf_weight(m, nd, xyd, r0, phi, zd);
            //
            //  #1:  Does interpolant match function at interpolation points?
            //
            ni = nd;
            xyi = typeMethods.r8mat_copy_new(2, ni, xyd);

            zi = RadialBasisFunctions.rbf_interp(m, nd, xyd, r0, phi, w, ni, xyi);

            int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / (double) (ni);

            Console.WriteLine("");
            Console.WriteLine("  L2 interpolation error averaged per interpolant node = "
                              + int_error + "");

        }
    }
}