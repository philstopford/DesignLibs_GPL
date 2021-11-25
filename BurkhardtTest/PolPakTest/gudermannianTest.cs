using System;
using System.Globalization;
using Burkardt.Function;

namespace PolPakTest;

public static class gudermannianTest
{
    public static void gud_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUD_TEST tests GUD.
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
        double fx = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("GUD_TEST:");
        Console.WriteLine("  GUD evaluates the Gudermannian function.");
        Console.WriteLine("");
        Console.WriteLine("     X      Exact F       GUD(X)");
        Console.WriteLine("");

        int n_data = 0;

        for ( ; ; )
        {
            Burkardt.Values.Gudermannian.gud_values ( ref n_data, ref x, ref fx );

            if ( n_data == 0 )
            {
                break;
            }

            double fx2 = Gudermannian.gud ( x );

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)   + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(10)  + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void agud_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AGUD_TEST tests AGUD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("AGUD_TEST");
        Console.WriteLine("  AGUD computes the inverse Gudermannian;");
        Console.WriteLine("");
        Console.WriteLine("         X     GUD(X)     AGUD(GUD(X))");
        Console.WriteLine("");

        for ( i = 0; i <= 10; i++ )
        {
            double x = 1.0 + i / 5.0;
            double g = Gudermannian.gud ( x );
            double x2 = Gudermannian.agud ( g );

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + g.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)    + "");
        }

    }

}