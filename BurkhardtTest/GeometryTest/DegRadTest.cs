using System;
using System.Globalization;
using Burkardt;

namespace GeometryTest;

public static class DegRadTest
{
    public static void test0205 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0205 tests DEGREES_TO_RADIANS and RADIANS_TO_DEGREES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST0205");
        Console.WriteLine("  DEGREES_TO_RADIANS converts an angle from degrees");
        Console.WriteLine("    to radians;");
        Console.WriteLine("  RADIANS_TO_DEGREES converts an angle from radians");
        Console.WriteLine("    to degrees;");
        Console.WriteLine("");
        Console.WriteLine("  Degrees     Radians     Degrees");
        Console.WriteLine("");

        for ( i = -2; i <= 14; i++ )
        {
            double angle_deg = 30 * i;
            double angle_rad = Helpers.degrees_to_radians ( angle_deg );
            double angle_deg2 = Helpers.radians_to_degrees ( angle_rad );
            Console.WriteLine("  " + angle_deg.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + angle_rad.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + angle_deg2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

}