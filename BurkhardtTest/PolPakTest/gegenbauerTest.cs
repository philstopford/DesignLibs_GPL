using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace PolPakTest;

public static class gegenbauerTest
{
    public static void gegenbauer_poly_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_POLY_TEST tests GEGENBAUER_POLY.
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
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("GEGENBAUER_POLY_TEST");
        Console.WriteLine("  GEGENBAUER_POLY evaluates the Gegenbauer polynomials.");
        Console.WriteLine("");
        Console.WriteLine("        N       A       X       GPV      GEGENBAUER");
        Console.WriteLine("");

        int n_data = 0;

        for ( ; ; )
        {

            Burkardt.Values.Gegenbauer.gegenbauer_poly_values ( ref n_data, ref n, ref a, ref x, ref fx );

            if ( n_data == 0 )
            {
                break;
            }

            double[] c = new double[n+1];

            GegenbauerPolynomial.gegenbauer_poly ( n, a, x, c );
            double fx2 = c[n];

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6)   + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)   + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(10)   + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)  + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }

    }

}