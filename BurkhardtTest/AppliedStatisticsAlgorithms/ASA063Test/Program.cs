using System;
using Burkardt;
using Burkardt.AppliedStatistics;

namespace ASA063Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA063_TEST.
        //
        //  Discussion:
        //
        //    ASA063_TEST tests the ASA063 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA063_TEST:");
            Console.WriteLine("  Test the ASA063 library.");

            test01 ( );

            Console.WriteLine("");
            Console.WriteLine("ASA063_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
        
        
        static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of BETAIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            int ifault = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  BETAIN computes the incomplete Beta function.");
            Console.WriteLine("  Compare to tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("           A           B           X      "
                + "    FX                        FX2");
            Console.WriteLine("                                          "
                + "    (Tabulated)               (BETAIN)            DIFF");
            Console.WriteLine("");

            int n_data = 0;

            for ( ; ; )
            {
                Algorithms.beta_inc_values ( ref n_data, ref a, ref b, ref x, ref fx );

                if ( n_data == 0 )
                {
                    break;
                }

                double beta_log = Helpers.LogGamma ( a )
                                  + Helpers.LogGamma ( b )
                                  - Helpers.LogGamma ( a + b );

                double fx2 = Algorithms.betain ( x, a, b, beta_log, ref ifault );

                Console.WriteLine("  " + a.ToString("0.####").PadLeft(10)
                    + "  " + b.ToString("0.####").PadLeft(10)
                    + "  " + x.ToString("0.####").PadLeft(10)
                    + "  " + fx.ToString("0.################").PadLeft(24)
                    + "  " + fx2.ToString("0.################").PadLeft(24)
                    + "  " + Math.Abs(( fx - fx2 )).ToString("0.####").PadLeft(10) + "");
            }
        }
        
    }
}