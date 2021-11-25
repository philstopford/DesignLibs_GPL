using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class ZetaTest
{
    public static void zeta_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZETA_VALUES_TEST tests ZETA_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 0;
        double zeta = 0;
        Console.WriteLine("");
        Console.WriteLine("ZETA_VALUES_TEST:");
        Console.WriteLine("  ZETA_VALUES returns values of ");
        Console.WriteLine("  the Riemann Zeta function.");
        Console.WriteLine("");
        Console.WriteLine("     N        ZETA(N)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            Zeta.zeta_values ( ref n_data, ref n, ref zeta );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6)                 + "  "
                                                                      + zeta.ToString("0.################").PadLeft(24) + "");
        }
    }
    public static void zeta_m1_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZETA_M1_VALUES_TEST tests ZETA_M1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double p = 0;
        double zeta_m1 = 0;
        Console.WriteLine("");
        Console.WriteLine("ZETA_M1_VALUES_TEST:");
        Console.WriteLine("  ZETA_M1_VALUES returns values of ");
        Console.WriteLine("  the Riemann Zeta Minus One function.");
        Console.WriteLine("");
        Console.WriteLine("     N        ZETA_M1(N)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            Zeta.zeta_m1_values ( ref n_data, ref p, ref zeta_m1 );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine(p.ToString(CultureInfo.InvariantCulture).PadLeft(6)                         + "  "
                                                                              + zeta_m1.ToString("0.################").PadLeft(24) + "");
        }
    }

}