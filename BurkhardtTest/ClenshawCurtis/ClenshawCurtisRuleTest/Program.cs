using System;
using Burkardt.ClenshawCurtisNS;

namespace ClenshawCurtisRuleTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CLENSHAW_CURTIS_RULE.
        //
        //  Discussion:
        //
        //    This program computes a standard Clenshaw Curtis quadrature rule
        //    and writes it to a file.
        //
        //    The user specifies:
        //    * the ORDER (number of points) in the rule;
        //    * A, the left endpoint;
        //    * B, the right endpoint;
        //    * FILENAME, which defines the output filenames.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 February 2010
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
        Console.WriteLine("CLENSHAW_CURTIS_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Compute a Clenshaw Curtis rule for approximating");
        Console.WriteLine("");
        Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx");
        Console.WriteLine("");
        Console.WriteLine("  of order ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies ORDER, A, B and FILENAME.");
        Console.WriteLine("");
        Console.WriteLine("  ORDER is the number of points.");
        Console.WriteLine("");
        Console.WriteLine("  A is the left endpoint.");
        Console.WriteLine("");
        Console.WriteLine("  B is the right endpoint.");
        Console.WriteLine("");
        Console.WriteLine("  FILENAME is used to generate 3 files:");
        Console.WriteLine("");
        Console.WriteLine("    filename_w.txt - the weight file");
        Console.WriteLine("    filename_x.txt - the abscissa file.");
        Console.WriteLine("    filename_r.txt - the region file.");
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
        Console.WriteLine("  FILENAME = \"" + filename + "\".");
        //
        //  Construct the rule.
        //
        double[] r = new double[2];
        double[] w = new double[order];
        double[] x = new double[order];

        r[0] = a;
        r[1] = b;

        ClenshawCurtis.clenshaw_curtis_compute(order, ref x, ref w);
        //
        //  Rescale the rule.
        //
        ClenshawCurtis.rescale(a, b, order, ref x, ref w);
        //
        //  Output the rule.
        //
        ClenshawCurtis.rule_write(order, filename, x, w, r);

        Console.WriteLine("");
        Console.WriteLine("CLENSHAW_CURTIS_RULE:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}