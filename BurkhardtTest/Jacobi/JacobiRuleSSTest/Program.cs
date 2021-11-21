using System;
using Burkardt.Quadrature;

namespace JacobiRuleSSTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for JACOBI_RULE_SS.
        //
        //  Discussion:
        //
        //    This program computes a standard Gauss-Jacobi quadrature rule
        //    and writes it to a file.
        //
        //    The user specifies:
        //    * the ORDER (number of points) in the rule
        //    * ALPHA and BETA parameters,
        //    * the root name of the output files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha;
        double beta;
        int order;
        string output;

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("JACOBI_RULE_SS");
        Console.WriteLine("");
        Console.WriteLine("  Compute a Gauss-Jacobi quadrature rule for approximating");
        Console.WriteLine("");
        Console.WriteLine("    Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
        Console.WriteLine("");
        Console.WriteLine("  of order ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies ORDER, ALPHA, BETA, and OUTPUT.");
        Console.WriteLine("");
        Console.WriteLine("  OUTPUT is:");
        Console.WriteLine("");
        Console.WriteLine("  \"C++\" for printed C++ output;");
        Console.WriteLine("  \"F77\" for printed Fortran77 output;");
        Console.WriteLine("  \"F90\" for printed Fortran90 output;");
        Console.WriteLine("  \"MAT\" for printed MATLAB output;");
        Console.WriteLine("");
        Console.WriteLine("  or:");
        Console.WriteLine("");
        Console.WriteLine("  \"filename\" to generate 3 files:");
        Console.WriteLine("");
        Console.WriteLine("    filename_w.txt - the weight file");
        Console.WriteLine("    filename_x.txt - the abscissa file.");
        Console.WriteLine("    filename_r.txt - the region file.");
        //
        //  Get the order.
        //
        try
        {
            order = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of ORDER (1 or greater)");
            order = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested order of the rule is = " + order + "");
        //
        //  Get ALPHA.
        //
        try
        {
            alpha = Convert.ToDouble(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  ALPHA is the exponent of (1-x) in the integral:");
            Console.WriteLine("");
            Console.WriteLine("      Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
            Console.WriteLine("");
            Console.WriteLine("  Note that -1.0 < ALPHA is required.");
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of ALPHA:");
            alpha = Convert.ToDouble(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested value of ALPHA = " + alpha + "");
        //
        //  Get BETA.
        //
        try
        {
            beta = Convert.ToDouble(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  BETA is the exponent of (1+x) in the integral:");
            Console.WriteLine("");
            Console.WriteLine("      Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
            Console.WriteLine("");
            Console.WriteLine("  Note that -1.0 < BETA is required.");
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of BETA:");
            beta = Convert.ToDouble(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested value of BETA = " + beta + "");
        //
        //  Get the output option or quadrature file root name:
        //
        try
        {
            output = args[3];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter OUTPUT (one of C++, F77, F90, MAT");
            Console.WriteLine("  or else the \"root name\" of the quadrature files).");
            output = Console.ReadLine();
        }

        Console.WriteLine("");
        Console.WriteLine("  OUTPUT option is \"" + output + "\".");
        //
        //  Construct the rule and output it.
        //
        JacobiQuadrature.jacobi_handle(order, alpha, beta, output);

        Console.WriteLine("");
        Console.WriteLine("JACOBI_RULE_SS:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }
}