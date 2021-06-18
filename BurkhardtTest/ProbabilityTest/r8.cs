using System;
using Burkardt.Probability;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ProbabilityTest
{
    partial class Program
    {
        static void r8_beta_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_BETA_TEST tests R8_BETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double x = 0;
            double y = 0;
            double fxy1 = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_BETA_TEST:");
            Console.WriteLine("  R8_BETA evaluates the Beta function.");
            Console.WriteLine("");
            Console.WriteLine("      X              Y         BETA(X,Y)         R8_BETA(X,Y)");
            Console.WriteLine("                               tabulated         computed.");
            Console.WriteLine("");

            var n_data = 0;

            for (;;)
            {
                Beta.beta_values(ref n_data, ref x, ref y, ref fxy1);

                if (n_data == 0)
                {
                    break;
                }

                var fxy2 = typeMethods.r8_beta(x, y);

                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                       + "  " + y.ToString().PadLeft(12)
                                       + "  " + fxy1.ToString("0.################").PadLeft(24)
                                       + "  " + fxy2.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void r8_ceiling_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_CEILING_TEST tests R8_CEILING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
        {
            Console.WriteLine("");
            Console.WriteLine("R8_CEILING_TEST");
            Console.WriteLine("  R8_CEILING rounds an R8 up.");
            Console.WriteLine("");
            Console.WriteLine("       X           R8_CEILING(X)");
            Console.WriteLine("");

            for (int i = -6; i <= 6; i++)
            {
                var rval = (double) (i) / 5.0;
                var ival = (int)Math.Ceiling(rval);
                Console.WriteLine("  "
                                  + rval.ToString().PadLeft(14) + "  "
                                  + ival.ToString().PadLeft(6) + "");
            }
        }

        static void r8_error_f_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERROR_F_TEST tests R8_ERROR_F.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int i;
            int seed;
            double x;
            double y;
            double z;

            Console.WriteLine("");
            Console.WriteLine("R8_ERROR_F_TEST");
            Console.WriteLine("  R8_ERROR_F evaluates ERF(X).");
            Console.WriteLine("");
            Console.WriteLine("X   -> Y = R8_ERROR_F(X) -> Z = R8_ERROR_F_INVERSE(Y)");
            Console.WriteLine("");

            seed = 123456789;

            x = 1.0;

            for (i = 1; i <= 20; i++)
            {
                x = Normal.normal_01_sample(ref seed);
                y = typeMethods.r8_error_f(x);
                z = typeMethods.r8_error_f_inverse(y);
                Console.WriteLine("  " + x.ToString().PadLeft(14)
                                  + "  " + y.ToString().PadLeft(14)
                                  + "  " + z.ToString().PadLeft(14) + "");
            }

            return;
        }

        static void r8_factorial_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_TEST tests R8_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double f;
            int i;

            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL_TEST");
            Console.WriteLine("  R8_FACTORIAL evaluates the factorial function.");
            Console.WriteLine("");
            Console.WriteLine("    I                R8_FACTORIAL(I)");
            Console.WriteLine("");

            for (i = 0; i <= 20; i++)
            {
                f = typeMethods.r8_factorial(i);

                Console.WriteLine("  "
                                  + i.ToString().PadLeft(4) + "  "
                                  + f.ToString().PadLeft(24) + "");
            }

            return;
        }

        static void r8_gamma_inc_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_INC_TEST tests R8_GAMMA_INC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
        {
            double a = 0;
            double fx = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_GAMMA_INC_TEST:");
            Console.WriteLine("  R8_GAMMA_INC evaluates the normalized incomplete Gamma");
            Console.WriteLine("  function.");
            Console.WriteLine("");
            Console.WriteLine("   A      X       Exact F       R8_GAMMA_INC(A,X)");
            Console.WriteLine("");

            var n_data = 0;

            for (;;)
            {
                Gamma.gamma_inc_values(ref n_data, ref a, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                var fx2 = typeMethods.r8_gamma_inc(a, x);

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(16) + "  "
                                  + fx2.ToString().PadLeft(16) + "");
            }

            return;
        }

        static void r8_gamma_log_int_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA_LOG_INT_TEST tests R8_GAMMA_LOG_INT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            Console.WriteLine("");
            Console.WriteLine("R8_GAMMA_LOG_INT_TEST");
            Console.WriteLine("  R8_GAMMA_LOG_INT evaluates the logarithm of the gamma function");
            Console.WriteLine("  for integer argument.");

            Console.WriteLine("");
            Console.WriteLine("       I     R8_GAMMA_LOG_INT(I)");
            Console.WriteLine("");

            for (int i = 1; i <= 20; i++)
            {
                var g = typeMethods.r8_gamma_log_int(i);

                Console.WriteLine("  "
                                  + i.ToString().PadLeft(6) + "  "
                                  + g.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void r8_uniform_01_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
        {
            int N = 1000;

            int i;
            double max;
            double mean;
            double min;
            int seed = 123456789;
            double[] x = new double[N];
            double variance;

            Console.WriteLine("");
            Console.WriteLine("R8_UNIFORM_01_TEST");
            Console.WriteLine("  R8_UNIFORM_01 samples a uniform random distribution in [0,1].");
            Console.WriteLine("  distributed random numbers.");
            Console.WriteLine("  Using initial random number seed = " + seed + "");

            for (i = 0; i < N; i++)
            {
                x[i] = UniformRNG.r8_uniform_01(ref seed);
            }

            Console.WriteLine("");
            Console.WriteLine("  First few values:");
            Console.WriteLine("");
            for (i = 0; i < 10; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(6)
                                  + "  " + x[i].ToString().PadLeft(14) + "");
            }

            min = typeMethods.r8vec_min(N, x);
            max = typeMethods.r8vec_max(N, x);
            mean = typeMethods.r8vec_mean(N, x);
            variance = typeMethods.r8vec_variance(N, x);

            Console.WriteLine("");
            Console.WriteLine("  Number of samples was " + N + "");
            Console.WriteLine("  Minimum value was " + min + "");
            Console.WriteLine("  Maximum value was " + max + "");
            Console.WriteLine("  Average value was " + mean + "");
            Console.WriteLine("  Variance was      " + variance + "");
        }

        static void r8_zeta_test()

//****************************************************************************80
//
//  Purpose:
//
//    R8_ZETA_TEST tests R8_ZETA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int i;
            double p;
            double v;

            Console.WriteLine("");
            Console.WriteLine("R8_ZETA_TEST");
            Console.WriteLine("  R8_ZETA estimates the Zeta function.");

            Console.WriteLine("");
            Console.WriteLine("       P     R8_Zeta(P)");
            Console.WriteLine("");
            for (i = 1; i <= 25; i++)
            {
                p = (double) (i);
                v = typeMethods.r8_zeta(p);
                Console.WriteLine("  " + p.ToString().PadLeft(6)
                                  + "  " + v.ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            for (i = 0; i <= 8; i++)
            {
                p = 3.0 + (double) (i) / 8.0;
                v = typeMethods.r8_zeta(p);
                Console.WriteLine("  " + p.ToString().PadLeft(6)
                                  + "  " + v.ToString().PadLeft(14) + "");
            }
        }
    }
}