using System;
using Burkardt.Laguerre;

namespace LaguerreRuleStroudSecrestTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGUERRE_RULE_SS.
        //
        //  Discussion:
        //
        //    This program computes a standard or exponentially weighted 
        //    Gauss-Laguerre quadrature rule and writes it to a file.
        //
        //    The user specifies:
        //    * the ORDER (number of points) in the rule
        //    * the OPTION (standard or reweighted rule)
        //    * the root name of the output files.
        //
        //    The parameter A is meant to represent the left endpoint of the interval
        //    of integration, and is currently fixed at 0.  It might be useful to allow
        //    the user to input a nonzero value of A.  In that case, the weights and
        //    abscissa would need to be appropriately adjusted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0.0;
        int option;
        int order;
        string output;

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_RULE_SS");
        Console.WriteLine("");
        Console.WriteLine("  Compute a Gauss-Laguerre rule for approximating");
        Console.WriteLine("");
        Console.WriteLine("    Integral ( A <= x < +oo ) exp(-x) f(x) dx");
        Console.WriteLine("");
        Console.WriteLine("  of order ORDER.");
        Console.WriteLine("");
        Console.WriteLine("  For now, A is fixed at 0.0.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies ORDER, OPTION, and OUTPUT.");
        Console.WriteLine("");
        Console.WriteLine("  OPTION is:");
        Console.WriteLine("");
        Console.WriteLine("    0 to get the standard rule for handling:");
        Console.WriteLine("      Integral ( A <= x < +oo ) exp(-x) f(x) dx");
        Console.WriteLine("");
        Console.WriteLine("    1 to get the modified rule for handling:");
        Console.WriteLine("      Integral ( A <= x < +oo )         f(x) dx");
        Console.WriteLine("");
        Console.WriteLine("    For OPTION = 1, the weights of the standard rule");
        Console.WriteLine("    are multiplied by exp(+x).");
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
        //  Get the option.
        //
        try
        {
            option = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  OPTION = 0 to get the standard rule for handling:");
            Console.WriteLine("      Integral ( A <= x < +oo ) exp(-x) f(x) dx");
            Console.WriteLine("");
            Console.WriteLine("  OPTION = 1 to get the modified rule for handling:");
            Console.WriteLine("      Integral ( A <= x < +oo )         f(x) dx");
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of OPTION:");
            option = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested value of OPTION = " + option + "");
        //
        //  Get the output option or quadrature file root name:
        //
        try
        {
            output = args[2];
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
        QuadratureRule.laguerre_handle(order, a, option, output);
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_RULE_SS:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}