using System;
using Burkardt;

namespace GeometryTest
{
    public static class dmsradTest
    {
        public static void test0235 ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0235 tests DMS_TO_RADIANS and RADIANS_TO_DMS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int angle_deg = 0;
            int angle_min = 0;
            double angle_rad = 0;
            double angle_rad2 = 0;
            int angle_sec = 0;
            int i;

            Console.WriteLine("");
            Console.WriteLine("TEST0235");
            Console.WriteLine("  DMS_TO_RADIANS converts an angle from");
            Console.WriteLine("    degrees/minutes/seconds to radians;");
            Console.WriteLine("  RADIANS_TO_DEGREES converts an angle from radians");
            Console.WriteLine("    to degrees/minutes/seconds;");
            Console.WriteLine("");
            Console.WriteLine("  Radians     DMS     Radians");
            Console.WriteLine("");

            for ( i = -2; i <= 15; i++ )
            {
                angle_rad = Math.PI * ( double ) ( i ) / 7.0;

                Helpers.radians_to_dms ( angle_rad, ref angle_deg, ref angle_min, ref angle_sec );

                angle_rad2 = Helpers.dms_to_radians ( angle_deg, angle_min, angle_sec );

                Console.WriteLine("  " + angle_rad.ToString().PadLeft(10)
                                  + "  " + angle_deg.ToString().PadLeft(4)
                                  + "  " + angle_min.ToString().PadLeft(3)
                                  + "  " + angle_sec.ToString().PadLeft(3)
                                  + "  " + angle_rad2.ToString().PadLeft(10) + "");
            }
        }

    }
}