using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class SphereTest
    {
        public static void sphere_unit_area_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_AREA_VALUES_TEST tests SPHERE_UNIT_AREA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data;
            int n = 0;
            Console.WriteLine("");
            Console.WriteLine("SPHERE_UNIT_AREA_VALUES_TEST:");
            Console.WriteLine("  SPHERE_UNIT_AREA_VALUES stores values of");
            Console.WriteLine("  the area of the unit sphere in various dimensions.");
            Console.WriteLine("");
            Console.WriteLine("      N           AREA");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sphere.sphere_unit_area_values(ref n_data, ref n, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void sphere_unit_volume_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERE_UNIT_VOLUME_VALUES_TEST tests SPHERE_UNIT_VOLUME_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data;
            int n = 0;
            Console.WriteLine("");
            Console.WriteLine("SPHERE_UNIT_VOLUME_VALUES_TEST:");
            Console.WriteLine("  SPHERE_UNIT_VOLUME_VALUES stores values of");
            Console.WriteLine("  the volume of the unit sphere in various dimensions.");
            Console.WriteLine("");
            Console.WriteLine("      N           VOLUME");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sphere.sphere_unit_volume_values(ref n_data, ref n, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}