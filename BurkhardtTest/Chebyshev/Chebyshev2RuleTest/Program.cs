using System;
using Burkardt.Quadrature;

namespace Chebyshev2RuleTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CHEBYSHEV2_RULE.
        //
        //  Discussion:
        //
        //    This program computes a standard Gauss-Chebyshev type 2 quadrature rule
        //    and writes it to a file.
        //
        //    The user specifies:
        //    * the ORDER (number of points) in the rule;
        //    * A, the left endpoint.
        //    * B, the right endpoint.
        //    * FILENAME, which defines the names of the output files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        string filename;
        int order;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV2_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Compute a Gauss-Chebyshev type 2 rule for approximating");
        Console.WriteLine("    Integral ( A <= x <= B ) f(x) * sqrt ( ( x - A ) * ( B - x ) ) dx");
        Console.WriteLine("  of order ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies ORDER, A, B, and FILENAME.");
        Console.WriteLine("");
        Console.WriteLine("  ORDER is the number of points:");
        Console.WriteLine("  A is the left endpoint:");
        Console.WriteLine("  B is the right endpoint:");
        Console.WriteLine("  FILENAME is used to generate 3 files:");
        Console.WriteLine("    filename_w.txt - the weight file");
        Console.WriteLine("    filename_x.txt - the abscissa file.");
        Console.WriteLine("    filename_r.txt - the region file.");
        //
        //  Initialize parameters;
        //
        const double alpha = 0.0;
        const double beta = 0.0;
        //
        //  Get ORDER.
        //
        try
        {
            order = Convert.ToInt32(args[0]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of ORDER (1 or greater)");
            order = Convert.ToInt32(Console.ReadLine());
        }

        //
        //  Get A.
        //
        try
        {
            a = Convert.ToDouble(args[1]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the left endpoint A:");
            a = Convert.ToDouble(Console.ReadLine());
        }

        //
        //  Get B.
        //
        try
        {
            b = Convert.ToDouble(args[2]);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the right endpoint B:");
            b = Convert.ToDouble(Console.ReadLine());
        }

        //
        //  Get FILENAME:
        //
        try
        {
            filename = args[3];
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter FILENAME, the \"root name\" of the quadrature files).");
            filename = Console.ReadLine();
        }

        //
        //  Input summary.
        //
        Console.WriteLine("");
        Console.WriteLine("  ORDER = " + order + "");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("  FILENAME \"" + filename + "\".");
        //
        //  Construct the rule.
        //
        double[] w = new double[order];
        double[] x = new double[order];

        const int kind = 9;
        CGQF.cgqf(order, kind, alpha, beta, a, b, ref x, ref w);
        //
        //  Write the rule.
        //
        double[] r = new double[2];
        r[0] = a;
        r[1] = b;

        QuadratureRule.rule_write(order, filename, x, w, r);

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV2_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}