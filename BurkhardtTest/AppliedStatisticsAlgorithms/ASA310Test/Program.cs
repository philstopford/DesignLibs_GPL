using System;
using Burkardt.AppliedStatistics;

namespace ASA310Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA310_PRB.
        //
        //  Discussion:
        //
        //    ASA310_PRB tests the ASA310 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA310_PRB:");
        Console.WriteLine("  Test the ASA310 library.");

        test01 ( );

        Console.WriteLine("");
        Console.WriteLine("ASA310_PRB:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
            
    }

    private static void test01 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of NCBETA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double errmax = 1.0E-10;
        double fx = 0;
        int ifault = 0;
        double lambda = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  NCBETA computes the noncentral incomplete Beta function.");
        Console.WriteLine("  Compare to tabulated values.");
        Console.WriteLine("");
        Console.WriteLine("      A        B     LAMBDA        X      "
                          + "    FX                        FX2");
        Console.WriteLine("                                          "
                          + "    (Tabulated)               (NCBETA)            DIFF");
        Console.WriteLine("");

        int n_data = 0;

        for ( ; ; )
        {
            Algorithms.beta_noncentral_cdf_values ( ref n_data, ref a, ref b, ref lambda, ref x, ref fx );

            if ( n_data == 0 )
            {
                break;
            }

            double fx2 = Algorithms.ncbeta ( a, b, lambda, x, errmax, ref ifault );

            Console.WriteLine("  " + a.ToString("0.##").PadLeft(7)
                                   + "  " + b.ToString("0.##").PadLeft(7)
                                   + "  " + lambda.ToString("0.###").PadLeft(7)
                                   + "  " + x.ToString("0.####").PadLeft(10)
                                   + "  " + fx.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + Math.Abs ( fx - fx2 ).ToString("0.####").PadLeft(10) + "");
        }
    }
        
}