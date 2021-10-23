using System;
using Burkardt;

namespace GeometryTest
{
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
            double d;
            double hx;
            int test;
            int test_num = 12;
            double x;

            Console.WriteLine("");
            Console.WriteLine("TEST031");
            Console.WriteLine("  HAVERSINE computes the haversine of an angle.");
            Console.WriteLine("");
            Console.WriteLine("  Degrees  Radians  Haversine");
            Console.WriteLine("");

            for ( test = 0; test <= test_num; test++ )
            {
                x = ( double ) ( test ) * 2.0 * Math.PI / ( double ) ( test_num );
                d = Helpers.radians_to_degrees ( x );
                hx = Helpers.haversine ( x );
                Console.WriteLine("  " + d.ToString().PadLeft(8)
                                  + "  " + x.ToString().PadLeft(8)
                                  + "  " + hx.ToString().PadLeft(12) + "");
            }
        }
    }
}