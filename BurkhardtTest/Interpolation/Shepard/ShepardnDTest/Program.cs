using System;
using Burkardt.Interpolation;
using Burkardt.Types;
using Burkardt.Uniform;
using InterpTest;

namespace ShepardnDTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SHEPARD_INTERP_ND_TEST.
            //
            //  Discussion:
            //
            //    SHEPARD_INTERP_ND_TEST tests the SHEPARD_INTERP_ND library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int j;
            int m;
            int n1d;
            int nd;
            double p;
            double[] p_test = {1.0, 2.0, 4.0, 8.0};
            int p_test_num = 4;
            int prob;
            int prob_num;

            Console.WriteLine("");
            Console.WriteLine("SHEPARD_INTERP_ND_TEST:");
            Console.WriteLine("  Test the SHEPARD_INTERP_ND library.");
            Console.WriteLine("  The R8LIB library is needed.");
            Console.WriteLine("  This test also needs the TEST_INTERP_ND library.");
            //
            //  Look at Shepard interpolant on an irregular grid.
            //
            nd = 25;

            prob_num = Data_nD.p00_prob_num();

            for (prob = 1; prob <= prob_num; prob++)
            {
                for (m = 2; m <= 5; m = m + 3)
                {
                    for (j = 0; j < p_test_num; j++)
                    {
                        p = p_test[j];
                        test01(prob, p, m, nd);
                    }

                }
            }

            //
            //  Look at Shepard interpolant on a regular N1D^M grid.
            //
            n1d = 5;

            for (prob = 1; prob <= prob_num; prob++)
            {
                for (m = 2; m <= 5; m = m + 3)
                {
                    for (j = 0; j < p_test_num; j++)
                    {
                        p = p_test[j];
                        test02(prob, p, m, n1d);
                    }
                }
            }

            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("SHEPARD_INTERP_ND_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01(int prob, double p, int m, int nd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests SHEPARD_INTERP on an irregular grid.
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
            //    Input, double P, the power used in the distance weighting.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int ND, the number of data points.
            //
        {
            double app_error;
            double[] c;
            double int_error;
            int ni;
            int seed;
            double[] w;
            double[] xd;
            double[] xi;
            double[] zd;
            double[] ze;
            double[] zi;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Interpolate data from TEST_INTERP_ND problem #" + prob + "");
            Console.WriteLine("  using Shepard interpolation with P = " + p + "");
            Console.WriteLine("  spatial dimension M = " + m + "");
            Console.WriteLine("  and an irregular grid of ND = " + nd + " data points.");
            //
            //  Set problem parameters:
            //
            seed = 123456789;
            c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
            w = UniformRNG.r8vec_uniform_01_new(m, ref seed);

            xd = UniformRNG.r8mat_uniform_01_new(m, nd, ref seed);

            zd = Data_nD.p00_f(prob, m, c, w, nd, xd);
            //
            //  #1:  Does interpolant match function at interpolation points?
            //
            ni = nd;
            xi = typeMethods.r8mat_copy_new(m, ni, xd);
            zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

            int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / (double) (ni);

            Console.WriteLine("");
            Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

            //
            //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
            //
            ni = 1000;
            ni = 50;
            xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);
            zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);
            ze = Data_nD.p00_f(prob, m, c, w, ni, xi);

            app_error = typeMethods.r8vec_norm_affine(ni, zi, ze) / (double) (ni);

            Console.WriteLine("  L2 approximation error averaged per 1000 samples =     " + app_error + "");

        }

        static void test02(int prob, double p, int m, int n1d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests SHEPARD_INTERP_ND on a regular N1D^M grid.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 October 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int PROB, the problem number.
            //
            //    Input, double P, the power used in the distance weighting.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N1D, the number of points in 1D.
            //
        {
            double a;
            double app_error;
            double b;
            double[] c;
            int i;
            double int_error;
            int nd;
            int ni;
            int seed;
            double[] w;
            double[] x1d;
            double[] xd;
            double[] xi;
            double[] zd;
            double[] ze;
            double[] zi;
            //
            //  Set problem parameters:
            //
            seed = 123456789;
            c = UniformRNG.r8vec_uniform_01_new(m, ref seed);
            w = UniformRNG.r8vec_uniform_01_new(m, ref seed);

            nd = (int) Math.Pow(n1d, m);

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Interpolate data from TEST_INTERP_ND problem #" + prob + "");
            Console.WriteLine("  using Shepard interpolation with P = " + p + "");
            Console.WriteLine("  spatial dimension M = " + m + "");
            Console.WriteLine("  and a regular grid of N1D^M = " + nd + " data points.");

            a = 0.0;
            b = 1.0;

            x1d = typeMethods.r8vec_linspace_new(n1d, a, b);

            xd = new double[m * nd];
            for (i = 0; i < m; i++)
            {
                typeMethods.r8vecDPData data = new typeMethods.r8vecDPData();
                typeMethods.r8vec_direct_product(ref data, i, n1d, x1d, m, nd, ref xd);
            }

            zd = Data_nD.p00_f(prob, m, c, w, nd, xd);
            //
            //  #1:  Does interpolant match function at interpolation points?
            //
            ni = nd;
            xi = typeMethods.r8mat_copy_new(m, nd, xd);
            zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

            int_error = typeMethods.r8vec_norm_affine(ni, zi, zd) / (double) (ni);

            Console.WriteLine("");
            Console.WriteLine("  L2 interpolation error averaged per interpolant node = " + int_error + "");

            //
            //  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
            //
            ni = 1000;
            xi = UniformRNG.r8mat_uniform_01_new(m, ni, ref seed);

            zi = Shepard.shepard_interp_nd(m, nd, xd, zd, p, ni, xi);

            ze = Data_nD.p00_f(prob, m, c, w, ni, xi);

            app_error = typeMethods.r8vec_norm_affine(ni, zi, ze) / (double) (ni);

            Console.WriteLine("  L2 approximation error averaged per 1000 samples =     " + app_error + "");

        }
    }
}