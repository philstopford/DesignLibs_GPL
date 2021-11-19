using System;
using Burkardt.PyramidNS;

namespace PyramidRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PYRAMID_RULE.
        //
        //  Discussion:
        //
        //    This program computes a quadrature rule for a pyramid
        //    and writes it to a file.
        //
        //    The user specifies:
        //    * the LEGENDRE_ORDER (number of points in the X and Y dimensions)
        //    * the JACOBI_ORDER (number of points in the Z dimension)
        //    * FILENAME, the root name of the output files.
        //
        //    The integration region is:
        //
        //      - ( 1 - Z ) <= X <= 1 - Z
        //      - ( 1 - Z ) <= Y <= 1 - Z
        //                0 <= Z <= 1.
        //
        //    When Z is zero, the integration region is a square lying in the (X,Y) 
        //    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
        //    radius of the square diminishes, and when Z reaches 1, the square has 
        //    contracted to the single point (0,0,1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string filename;
        int jacobi_order;
        int legendre_order;

        Console.WriteLine("");
        Console.WriteLine("PYRAMID_RULE");
        Console.WriteLine("");
        Console.WriteLine("  Compute a quadrature rule for approximating");
        Console.WriteLine("  the integral of a function over a pyramid.");
        Console.WriteLine("");
        Console.WriteLine("  The user specifies:");
        Console.WriteLine("");
        Console.WriteLine("  LEGENDRE_ORDER, the order of the Legendre rule for X and Y.");
        Console.WriteLine("  JACOBI_ORDER, the order of the Jacobi rule for Z,");
        Console.WriteLine("  FILENAME, the prefix of the three output files:");
        Console.WriteLine("");
        Console.WriteLine("    filename_w.txt - the weight file");
        Console.WriteLine("    filename_x.txt - the abscissa file.");
        Console.WriteLine("    filename_r.txt - the region file.");
        //
        //  Get the Legendre order.
        //
        try
        {
            legendre_order = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the Legendre rule order:");
            legendre_order = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested Legendre order of the rule is " + legendre_order + "");
        //
        //  Get the Jacobi order.
        //
        try
        {
            jacobi_order = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the Jacobi rule order:");
            jacobi_order = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  The requested Jacobi order of the rule is " + jacobi_order + "");
        //
        //  Get the output option or quadrature file root name:
        //
        try
        {
            filename = args[2];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter FILENAME, the root name of the quadrature files.");
            filename = Console.ReadLine();
        }

        QuadratureRule.pyramid_handle(legendre_order, jacobi_order, filename);

        Console.WriteLine("");
        Console.WriteLine("PYRAMID_RULE:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}