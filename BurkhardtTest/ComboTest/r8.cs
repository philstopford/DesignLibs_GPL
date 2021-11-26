using System;
using System.Globalization;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void r8_choose_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_CHOOSE_TEST tests R8_CHOOSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("R8_CHOOSE_TEST");
        Console.WriteLine("  R8_CHOOSE evaluates C(N,K).");
        Console.WriteLine("");
        Console.WriteLine("         N         K       CNK");
 
        for ( n = 0; n <= 5; n++ )
        {
            Console.WriteLine("");
            int k;
            for ( k = 0; k <= n; k++ )
            {
                double cnk = typeMethods.r8_choose ( n, k );
                Console.WriteLine(n.ToString().PadLeft(10) + "  "
                                                           + k.ToString().PadLeft(8) + "  "
                                                           + cnk.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void r8_gamma_log_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_GAMMA_LOG_TEST tests R8_GAMMA_LOG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("R8_GAMMA_LOG_TEST:");
        Console.WriteLine("   R8_GAMMA_LOG evaluates the logarithm of the Gamma function.");
        Console.WriteLine("");
        Console.WriteLine("      X            GAMMA_LOG(X)     R8_GAMMA_LOG(X)");
        Console.WriteLine("");

        int n_data = 0;

        for ( ; ; )
        {
            Burkardt.Values.Gamma.gamma_log_values ( ref n_data, ref x, ref fx1 );

            if ( n_data == 0 )
            {
                break;
            }
            double fx2 = typeMethods.r8_gamma_log ( x );
            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx1.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24) + "");
        }
    }
}