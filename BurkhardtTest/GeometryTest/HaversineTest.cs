using System;
using System.Globalization;
using Burkardt;

namespace GeometryTest;

public static class HaversineTest
{
    public static void test031 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST031 tests HAVERSINE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("TEST031");
        Console.WriteLine("  HAVERSINE computes the haversine of an angle.");
        Console.WriteLine("");
        Console.WriteLine("  Degrees  Radians  Haversine");
        Console.WriteLine("");

        for ( test = 0; test <= test_num; test++ )
        {
            double x = test * 2.0 * Math.PI / test_num;
            double d = Helpers.radians_to_degrees ( x );
            double hx = Helpers.haversine ( x );
            Console.WriteLine("  " + d.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + hx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}