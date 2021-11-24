using System;
using Burkardt.Square;

namespace SquareFelippaRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SQUARE_FELIPPA_RULE_TEST.
        //
        //  Discussion:
        //
        //    SQUARE_FELIPPA_RULE_TEST tests the SQUARE_FELIPPA_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SQUARE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Test the SQUARE_FELIPPA_RULE library.");

        int degree_max = 4;
        FelippaRule.square_monomial_test(degree_max);

        degree_max = 5;
        FelippaRule.square_quad_test(degree_max);

        Console.WriteLine("");
        Console.WriteLine("SQUARE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}