using System;
using Burkardt.ClenshawCurtisNS;
using Burkardt.OrderNS;
using Burkardt.Quadrature;

namespace PattersonRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PATTERSON_RULE.
        //
        //  Discussion:
        //
        //    This program computes a standard Gauss-Patterson quadrature rule
        //    and writes it to a file.
        //
        //  Usage:
        //
        //    patterson_rule order a b output
        //
        //    where
        //
        //    * ORDER is the number of points in the rule, and must be
        //      1, 3, 7, 15, 31, 63, 127, 255 or 511.
        //    * A is the left endpoint;
        //    * B is the right endpoint;
        //    * FILENAME is the "root name" of the output files.
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
        double[] r;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("PATTERSON_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Compute a Gauss-Patterson rule for approximating");
        Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx");
        Console.WriteLine("  of order ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies ORDER, A, B, and FILENAME.");
        Console.WriteLine("");
        Console.WriteLine("  ORDER is 1, 3, 7, 15, 31, 63, 127, 255 or 511.");
        Console.WriteLine("  A is the left endpoint.");
        Console.WriteLine("  B is the right endpoint.");
        Console.WriteLine("  FILENAME is used to generate 3 files:");
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
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of ORDER.");
            order = Convert.ToInt32(Console.ReadLine());
        }

        if (!Order.order_check(order))
        {
            Console.WriteLine("");
            Console.WriteLine("PATTERSON_RULE:");
            Console.WriteLine("  ORDER is illegal.");
            Console.WriteLine("  Abnormal end of execution.");
            return;
        }

        //
        //  Get A.
        //
        try
        {
            a = Convert.ToDouble(args[1]);
        }
        catch
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
        catch
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
        Console.WriteLine("  ORDER = " + order + "");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("  FILENAME = \"" + filename + "\".");
        //
        //  Construct the rule.
        //
        r = new double[2];
        w = new double[order];
        x = new double[order];

        r[0] = a;
        r[1] = b;

        PattersonQuadrature.patterson_set(order, ref x, ref w);
        //
        //  Rescale the rule.
        //
        ClenshawCurtis.rescale(a, b, order, ref x, ref w);
        //
        //  Output the rule.
        //
        QuadratureRule.rule_write(order, filename, x, w, r);

        Console.WriteLine("");
        Console.WriteLine("PATTERSON_RULE:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}