using System;
using Burkardt;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace BlackScholesTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for BLACK_SCHOLES_TEST.
            //
            //  Discussion:
            //
            //    BLACK_SCHOLES_TEST tests the BLACK_SCHOLES library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("BLACK_SCHOLES_TEST");
            Console.WriteLine("  Test the BLACK_SCHOLES library.");

            asset_path_test();
            binomial_test();
            bsf_test();
            forward_test();
            mc_test();

            Console.WriteLine("");
            Console.WriteLine("BLACK_SCHOLES_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void asset_path_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ASSET_PATH_TEST tests ASSET_PATH.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 June 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double mu;
            int n = 100;
            string output_filename;
            double[] s;
            double s0;
            int seed;
            double sigma;
            double t1;

            Console.WriteLine("");
            Console.WriteLine("ASSET_PATH_TEST:");
            Console.WriteLine("  Demonstrate the simulated of an asset price path.");

            s0 = 2.0;
            mu = 0.1;
            sigma = 0.3;
            t1 = 1.0;
            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            Console.WriteLine("");
            Console.WriteLine("  The asset price at time 0      S0    = " + s0 + "");
            Console.WriteLine("  The asset expected growth rate MU    = " + mu + "");
            Console.WriteLine("  The asset volatility           SIGMA = " + sigma + "");
            Console.WriteLine("  The expiry date                T1    = " + t1 + "");
            Console.WriteLine("  The number of time steps       N     = " + n + "");
            Console.WriteLine("  The random number seed was     SEED  = " + seed + "");

            s = BlackScholes.asset_path(s0, mu, sigma, t1, n, ref data, ref seed);

            typeMethods.r8vec_print_part(n + 1, s, 10, "  Partial results:");

            output_filename = "asset_path.txt";
            typeMethods.r8vec_write(output_filename, n + 1, s);

            Console.WriteLine("");
            Console.WriteLine("  Full results written to \"" + output_filename + "\".");

        }

        static void binomial_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_TEST tests BINOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c;
            double e;
            int m;
            double r;
            double s0;
            double sigma;
            double t1;

            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_TEST:");
            Console.WriteLine("  A demonstration of the binomial method");
            Console.WriteLine("  for option valuation.");

            s0 = 2.0;
            e = 1.0;
            r = 0.05;
            sigma = 0.25;
            t1 = 3.0;
            m = 256;

            Console.WriteLine("");
            Console.WriteLine("  The asset price at time 0 S0    = " + s0 + "");
            Console.WriteLine("  The exercise price        E     = " + e + "");
            Console.WriteLine("  The interest rate         R     = " + r + "");
            Console.WriteLine("  The asset volatility      SIGMA = " + sigma + "");
            Console.WriteLine("  The expiry date           T1    = " + t1 + "");
            Console.WriteLine("  The number of intervals   M     = " + m + "");

            c = BlackScholes.binomial(s0, e, r, sigma, t1, m);

            Console.WriteLine("");
            Console.WriteLine("  The option value is " + c + "");

        }

        static void bsf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BSF_TEST tests BSF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double c;
            double e;
            double r;
            double s0;
            double sigma;
            double t0;
            double t1;

            Console.WriteLine("");
            Console.WriteLine("BSF_TEST:");
            Console.WriteLine("  A demonstration of the Black-Scholes formula");
            Console.WriteLine("  for option valuation.");

            s0 = 2.0;
            t0 = 0.0;
            e = 1.0;
            r = 0.05;
            sigma = 0.25;
            t1 = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  The asset price at time T0 S0    = " + s0 + "");
            Console.WriteLine("  The time                   T0    = " + t0 + "");
            Console.WriteLine("  The exercise price         E     = " + e + "");
            Console.WriteLine("  The interest rate          R     = " + r + "");
            Console.WriteLine("  The asset volatility       SIGMA = " + sigma + "");
            Console.WriteLine("  The expiry date            T1    = " + t1 + "");

            c = BlackScholes.bsf(s0, t0, e, r, sigma, t1);

            Console.WriteLine("");
            Console.WriteLine("  The option value C = " + c + "");

        }

        static void forward_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FORWARD_TEST tests FORWARD.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double e;
            int i;
            int nt;
            int nx;
            double r;
            double s;
            double sigma;
            double smax;
            double smin;
            double t1;
            double[] u;

            Console.WriteLine("");
            Console.WriteLine("FORWARD_TEST:");
            Console.WriteLine("  A demonstration of the forward difference method");
            Console.WriteLine("  for option valuation.");

            e = 4.0;
            r = 0.03;
            sigma = 0.50;
            t1 = 1.0;
            nx = 11;
            nt = 29;
            smax = 10.0;

            Console.WriteLine("");
            Console.WriteLine("  The exercise price        E =     " + e + "");
            Console.WriteLine("  The interest rate         R =     " + r + "");
            Console.WriteLine("  The asset volatility      SIGMA = " + sigma + "");
            Console.WriteLine("  The expiry date           T1 =    " + t1 + "");
            Console.WriteLine("  The number of space steps NX =    " + nx + "");
            Console.WriteLine("  The number of time steps  NT =    " + nt + "");
            Console.WriteLine("  The value of              SMAX =  " + smax + "");

            u = BlackScholes.forward(e, r, sigma, t1, nx, nt, smax);

            Console.WriteLine("");
            Console.WriteLine("         Initial          Option");
            Console.WriteLine("           Value           Value");
            Console.WriteLine("");

            smin = 0.0;
            for (i = 0; i < nx - 1; i++)
            {
                s = ((nx - i - 2) * smin + (i + 1) * smax) / (double) (nx - 1);
                Console.WriteLine("  " + s.ToString().PadLeft(12)
                    + "  " + u[i + nt * (nx - 1)].ToString().PadLeft(14) + "");
            }
        }

        static void mc_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MC_TEST tests MC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 June 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] conf;
            double e;
            int m;
            double r;
            double s0;
            int seed;
            double sigma;
            double t1;

            Console.WriteLine("");
            Console.WriteLine("MC_TEST:");
            Console.WriteLine("  A demonstration of the Monte Carlo method");
            Console.WriteLine("  for option valuation.");

            s0 = 2.0;
            e = 1.0;
            r = 0.05;
            sigma = 0.25;
            t1 = 3.0;
            m = 1000000;
            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            Console.WriteLine("");
            Console.WriteLine("  The asset price at time 0, S0    = " + s0 + "");
            Console.WriteLine("  The exercise price         E     = " + e + "");
            Console.WriteLine("  The interest rate          R     = " + r + "");
            Console.WriteLine("  The asset volatility       SIGMA = " + sigma + "");
            Console.WriteLine("  The expiry date            T1    = " + t1 + "");
            Console.WriteLine("  The number of simulations  M     = " + m + "");
            Console.WriteLine("  The random number seed was SEED  = " + seed + "");

            conf = BlackScholes.mc(s0, e, r, sigma, t1, m, ref data, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  The confidence interval is [" + conf[0] + ", " + conf[1] + "].");

        }
    }
}