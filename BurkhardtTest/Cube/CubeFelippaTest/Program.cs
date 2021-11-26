using System;
using Burkardt.Cube;

namespace CubeFelippaTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CUBE_FELIPPA_RULE_TEST.
        //
        //  Discussion:
        //
        //    CUBE_FELIPPA_RULE_TEST tests the CUBE_FELIPPA_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CUBE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Test the CUBE_FELIPPA_RULE library.");

        int degree_max = 4;
        Integrals.cube_monomial_test(degree_max);

        degree_max = 6;
        QuadratureRule.cube_quad_test(degree_max);

        Console.WriteLine("");
        Console.WriteLine("CUBE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}