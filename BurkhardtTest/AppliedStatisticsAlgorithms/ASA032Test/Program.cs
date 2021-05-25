using System;
using Burkardt.AppliedStatistics;

namespace Burkardt.ASA032Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA032_TEST.
        //
        //  Discussion:
        //
        //    ASA032_TEST tests the ASA032 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA032_TEST:");
            Console.WriteLine("  Test the ASA032 library.");

            test01 ( );

            Console.WriteLine("");
            Console.WriteLine("ASA032_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
        
        
        static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of ALNGAM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double a = 0;
            double fx = 0;
            int ifault = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  GAMAIN computes the incomplete Gamma function.");
            Console.WriteLine("  Compare the result to tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("          A               X           " 
                + "FX                        FX2");
            Console.WriteLine("                                      "
                + "(Tabulated)               (GAMAIN)                DIFF");
            Console.WriteLine("");

            int n_data = 0;

            for ( ; ; )
            {
                Algorithms.gamma_inc_values ( ref n_data, ref a, ref x, ref fx );

                if ( n_data == 0 )
                {
                    break;
                }

                double fx2 = Algorithms.gamain ( x, a, ref ifault );

                Console.WriteLine("  " + a.ToString("0.########").PadLeft(12) 
                               + "  " + x.ToString("0.########").PadLeft(12)
                               + "  " + fx.ToString("0.################").PadLeft(24)
                               + "  " + fx2.ToString("0.################").PadLeft(24)
                               + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10)
                );
            }
        }
    }
}