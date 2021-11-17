using System;
using Burkardt.AppliedStatistics;

namespace ASA121Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA121_TEST.
        //
        //  Discussion:
        //
        //    ASA121_TEST tests the ASA121 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA121_TEST:");
        Console.WriteLine("  Test the ASA121 library.");

        test01 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA121_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of TRIGAMMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data = 0;
        double x = 0;
        int ifault = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  TRIGAMMA computes the trigamma function. ");
        Console.WriteLine("  We compare the result to tabulated values.");
        Console.WriteLine("");
        Console.WriteLine("          X                     "
                          + "FX                        FX2");
        Console.WriteLine("                                "
                          + "(Tabulated)               (TRIGAMMA)                DIFF");
        Console.WriteLine("");

        n_data = 0;

        for ( ; ; )
        {
            Algorithms.trigamma_values ( ref n_data, ref x, ref fx );

            if ( n_data == 0 )
            {
                break;
            }

            double fx2 = Algorithms.trigamma ( x, ref ifault );

            Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                   + "  " + fx.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + Math.Abs (fx - fx2).ToString("0.####").PadLeft(10) + "");
        }
    }
        
}