﻿using System;
using Burkardt.Line;

namespace LineFelippaTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LINE_FELIPPA_RULE_TEST.
            //
            //  Discussion:
            //
            //    LINE_FELIPPA_RULE_TEST tests the LINE_FELIPPA_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 September 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree_max;

            Console.WriteLine("");
            Console.WriteLine("LINE_FELIPPA_RULE_TEST");
            Console.WriteLine("  Test the LINE_FELIPPA_RULE library.");

            degree_max = 4;
            Felippa.line_monomial_test(degree_max);

            degree_max = 10;
            Felippa.line_quad_test(degree_max);

            Console.WriteLine("");
            Console.WriteLine("LINE_FELIPPA_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}