using System;
using System.Globalization;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace BlackScholesTest;

internal static class Program
{
    private static void Main()
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

    private static void asset_path_test()

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
        const int n = 100;

        Console.WriteLine("");
        Console.WriteLine("ASSET_PATH_TEST:");
        Console.WriteLine("  Demonstrate the simulated of an asset price path.");

        const double s0 = 2.0;
        const double mu = 0.1;
        const double sigma = 0.3;
        const double t1 = 1.0;
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("  The asset price at time 0      S0    = " + s0 + "");
        Console.WriteLine("  The asset expected growth rate MU    = " + mu + "");
        Console.WriteLine("  The asset volatility           SIGMA = " + sigma + "");
        Console.WriteLine("  The expiry date                T1    = " + t1 + "");
        Console.WriteLine("  The number of time steps       N     = " + n + "");
        Console.WriteLine("  The random number seed was     SEED  = " + seed + "");

        double[] s = BlackScholes.asset_path(s0, mu, sigma, t1, n, ref data, ref seed);

        typeMethods.r8vec_print_part(n + 1, s, 10, "  Partial results:");

        string output_filename = "asset_path.txt";
        typeMethods.r8vec_write(output_filename, n + 1, s);

        Console.WriteLine("");
        Console.WriteLine("  Full results written to \"" + output_filename + "\".");

    }

    private static void binomial_test()

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
        Console.WriteLine("");
        Console.WriteLine("BINOMIAL_TEST:");
        Console.WriteLine("  A demonstration of the binomial method");
        Console.WriteLine("  for option valuation.");

        const double s0 = 2.0;
        const double e = 1.0;
        const double r = 0.05;
        const double sigma = 0.25;
        const double t1 = 3.0;
        const int m = 256;

        Console.WriteLine("");
        Console.WriteLine("  The asset price at time 0 S0    = " + s0 + "");
        Console.WriteLine("  The exercise price        E     = " + e + "");
        Console.WriteLine("  The interest rate         R     = " + r + "");
        Console.WriteLine("  The asset volatility      SIGMA = " + sigma + "");
        Console.WriteLine("  The expiry date           T1    = " + t1 + "");
        Console.WriteLine("  The number of intervals   M     = " + m + "");

        double c = BlackScholes.binomial(s0, e, r, sigma, t1, m);

        Console.WriteLine("");
        Console.WriteLine("  The option value is " + c + "");

    }

    private static void bsf_test()

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
        Console.WriteLine("");
        Console.WriteLine("BSF_TEST:");
        Console.WriteLine("  A demonstration of the Black-Scholes formula");
        Console.WriteLine("  for option valuation.");

        const double s0 = 2.0;
        const double t0 = 0.0;
        const double e = 1.0;
        const double r = 0.05;
        const double sigma = 0.25;
        const double t1 = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  The asset price at time T0 S0    = " + s0 + "");
        Console.WriteLine("  The time                   T0    = " + t0 + "");
        Console.WriteLine("  The exercise price         E     = " + e + "");
        Console.WriteLine("  The interest rate          R     = " + r + "");
        Console.WriteLine("  The asset volatility       SIGMA = " + sigma + "");
        Console.WriteLine("  The expiry date            T1    = " + t1 + "");

        double c = BlackScholes.bsf(s0, t0, e, r, sigma, t1);

        Console.WriteLine("");
        Console.WriteLine("  The option value C = " + c + "");

    }

    private static void forward_test()

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
        int i;

        Console.WriteLine("");
        Console.WriteLine("FORWARD_TEST:");
        Console.WriteLine("  A demonstration of the forward difference method");
        Console.WriteLine("  for option valuation.");

        const double e = 4.0;
        const double r = 0.03;
        const double sigma = 0.50;
        const double t1 = 1.0;
        const int nx = 11;
        const int nt = 29;
        const double smax = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  The exercise price        E =     " + e + "");
        Console.WriteLine("  The interest rate         R =     " + r + "");
        Console.WriteLine("  The asset volatility      SIGMA = " + sigma + "");
        Console.WriteLine("  The expiry date           T1 =    " + t1 + "");
        Console.WriteLine("  The number of space steps NX =    " + nx + "");
        Console.WriteLine("  The number of time steps  NT =    " + nt + "");
        Console.WriteLine("  The value of              SMAX =  " + smax + "");

        double[] u = BlackScholes.forward(e, r, sigma, t1, nx, nt, smax);

        Console.WriteLine("");
        Console.WriteLine("         Initial          Option");
        Console.WriteLine("           Value           Value");
        Console.WriteLine("");

        double smin = 0.0;
        for (i = 0; i < nx - 1; i++)
        {
            double s = ((nx - i - 2) * smin + (i + 1) * smax) / (nx - 1);
            Console.WriteLine("  " + s.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + u[i + nt * (nx - 1)].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void mc_test()

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
        Console.WriteLine("");
        Console.WriteLine("MC_TEST:");
        Console.WriteLine("  A demonstration of the Monte Carlo method");
        Console.WriteLine("  for option valuation.");

        const double s0 = 2.0;
        const double e = 1.0;
        const double r = 0.05;
        const double sigma = 0.25;
        const double t1 = 3.0;
        const int m = 1000000;
        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("  The asset price at time 0, S0    = " + s0 + "");
        Console.WriteLine("  The exercise price         E     = " + e + "");
        Console.WriteLine("  The interest rate          R     = " + r + "");
        Console.WriteLine("  The asset volatility       SIGMA = " + sigma + "");
        Console.WriteLine("  The expiry date            T1    = " + t1 + "");
        Console.WriteLine("  The number of simulations  M     = " + m + "");
        Console.WriteLine("  The random number seed was SEED  = " + seed + "");

        double[] conf = BlackScholes.mc(s0, e, r, sigma, t1, m, ref data, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  The confidence interval is [" + conf[0] + ", " + conf[1] + "].");

    }
}