using System;
using Burkardt;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace TruncatedNormalRuleTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRUNCATED_NORMAL_RULE.
        //
        //  Discussion:
        //
        //    This program computes a truncated normal quadrature rule
        //    and writes it to a file.
        //
        //    The user specifies:
        //    * option: 0/1/2/3 for none, lower, upper, double truncation.
        //    * N, the number of points in the rule;
        //    * MU, the mean of the original normal distribution;
        //    * SIGMA, the standard deviation of the original normal distribution,
        //    * A, the left endpoint (for options 1 or 3)
        //    * B, the right endpoint (for options 2 or 3);
        //    * FILENAME, the root name of the output files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        string filename;
        int iarg;
        double[] moment = new double[1];
        double mu;
        int n;
        int option;
        double[] r;
        double sigma;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_RULE");
        Console.WriteLine("");
        Console.WriteLine("  For the (truncated) Gaussian probability density function");
        Console.WriteLine("    pdf(x) = exp(-0.5*((x-MU)/SIGMA)^2) / SIGMA / sqrt ( 2 * pi )");
        Console.WriteLine("  compute an N-point quadrature rule for approximating");
        Console.WriteLine("    Integral ( A <= x <= B ) f(x) pdf(x) dx");
        Console.WriteLine("");
        Console.WriteLine("  The value of OPTION determines the truncation interval [A,B]:");
        Console.WriteLine("  0: (-oo,+oo)");
        Console.WriteLine("  1: [A,+oo)");
        Console.WriteLine("  2: (-oo,B]");
        Console.WriteLine("  3: [A,B]");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies OPTION, N, MU, SIGMA, A, B and FILENAME.");
        Console.WriteLine("");
        Console.WriteLine("  FILENAME is used to generate 3 files:");
        Console.WriteLine("");
        Console.WriteLine("    filename_w.txt - the weight file");
        Console.WriteLine("    filename_x.txt - the abscissa file.");
        Console.WriteLine("    filename_r.txt - the region file, listing A and B.");

        iarg = 0;
        //
        //  Get OPTION.
        //
        try
        {
            option = Convert.ToInt32(args[iarg]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of OPTION 0/1/2/3:  ");
            option = Convert.ToInt32(Console.ReadLine());
        }

        switch (option)
        {
            case < 0:
            case > 3:
                Console.WriteLine("");
                Console.WriteLine("TRUNCATED_NORMAL_RULE - Fatal error!");
                Console.WriteLine("  0 <= OPTION <= 3 was required.");
                return;
        }

        //
        //  Get N.
        //
        iarg += 1;
        try
        {
            n = Convert.ToInt32(args[iarg]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of N (1 or greater)");
            n = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get MU.
        //
        iarg += 1;
        try
        {
            mu = Convert.ToDouble(args[iarg]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter MU, the mean value of the normal distribution:");
            mu = Convert.ToDouble(Console.ReadLine());
        }

        //
        //  Get SIGMA.
        //
        iarg += 1;
        try
        {
            sigma = Convert.ToDouble(args[iarg]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter SIGMA, the standard deviation of the normal distribution:");
            sigma = Convert.ToDouble(Console.ReadLine());
        }

        sigma = Math.Abs(sigma);
        switch (option)
        {
            //
            //  Get A.
            //
            case 1:
            case 3:
                iarg += 1;
                try
                {
                    a = Convert.ToDouble(args[iarg]);
                }
                catch
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Enter the left endpoint A:");
                    a = Convert.ToDouble(Console.ReadLine());
                }

                break;
            default:
                a = -typeMethods.r8_huge();
                break;
        }

        switch (option)
        {
            //
            //  Get B.
            //
            case 2:
            case 3:
                iarg += 1;
                try
                {
                    b = Convert.ToDouble(args[iarg]);
                }
                catch
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Enter the right endpoint B:");
                    b = Convert.ToDouble(Console.ReadLine());
                }

                break;
            default:
                b = typeMethods.r8_huge();
                break;
        }

        if (b <= a)
        {
            Console.WriteLine("");
            Console.WriteLine("TRUNCATED_NORMAL_RULE - Fatal error!");
            Console.WriteLine("  A < B required!");
            return;
        }

        //
        //  Get FILENAME:
        //
        iarg += 1;
        try
        {
            filename = args[iarg];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter FILENAME, the \"root name\" of the quadrature files).");
            filename = Console.ReadLine();
        }

        //
        //  Input summary.
        //
        Console.WriteLine("");
        Console.WriteLine("  OPTION = " + option + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("  MU = " + mu + "");
        Console.WriteLine("  SIGMA = " + sigma + "");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("  FILENAME = \"" + filename + "\".");
        moment = option switch
        {
            //
            //  Compute the moments.
            //
            0 => Moments.moments_normal(2 * n + 1, mu, sigma),
            1 => Moments.moments_truncated_normal_a(2 * n + 1, mu, sigma, a),
            2 => Moments.moments_truncated_normal_b(2 * n + 1, mu, sigma, b),
            3 => Moments.moments_truncated_normal_ab(2 * n + 1, mu, sigma, a, b),
            _ => moment
        };

        //
        //  Construct the rule from the moments.
        //
        w = new double[n];
        x = new double[n];

        Moments.moment_method(n, moment, ref x, ref w);
        //
        //  Output the rule.
        //
        r = new double[2];
        r[0] = a;
        r[1] = b;

        QuadratureRule.rule_write(n, filename, x, w, r);

        Console.WriteLine("");
        Console.WriteLine("TRUNCATED_NORMAL_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}