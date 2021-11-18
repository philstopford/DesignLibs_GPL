using System;
using Burkardt.Function;

namespace PolPakTest;

public static class lerchTest
{
    public static void lerch_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LERCH_TEST tests LERCH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx = 0;
        double fx2 = 0;
        int n_data = 0;
        int s = 0;
        double z = 0;

        Console.WriteLine("");
        Console.WriteLine("LERCH_TEST:");
        Console.WriteLine("  LERCH evaluates the Lerch function.");
        Console.WriteLine("");
        Console.WriteLine("       Z       S       A         Lerch           Lerch");
        Console.WriteLine("                             Tabulated        Computed");
        Console.WriteLine("");

        n_data = 0;

        for ( ; ; )
        {
            Burkardt.Values.Lerch.lerch_values ( ref n_data, ref z, ref s, ref a, ref fx );

            if ( n_data == 0 )
            {
                break;
            }

            fx2 = Lerch.lerch ( z, s, a );

            Console.WriteLine("  "
                              + z.ToString(CultureInfo.InvariantCulture).PadLeft(8)   + "  "
                              + s.ToString(CultureInfo.InvariantCulture).PadLeft(4)   + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)   + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)  + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

}