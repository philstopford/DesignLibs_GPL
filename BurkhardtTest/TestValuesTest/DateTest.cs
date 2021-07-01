using System;
using TestValues;

namespace TestValuesTest
{
    public static class DateTest
    {
        public static void datenum_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DATENUM_VALUES_TEST tests DATENUM_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 December 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d = 0;
            double date_num = 0;
            int m = 0;
            int n_data;
            int y = 0;
            Console.WriteLine("");
            Console.WriteLine("DATENUM_VALUES_TEST:");
            Console.WriteLine("  DATENUM_VALUES returns values of ");
            Console.WriteLine("  the MATLAB datenum for a given Y/M/D date.");
            Console.WriteLine("");
            Console.WriteLine("     Y     M     D  DATENUM");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.datenum_values(ref n_data, ref y, ref m, ref d, ref date_num);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + y.ToString().PadLeft(4)
                                       + "  " + m.ToString().PadLeft(4)
                                       + "  " + d.ToString().PadLeft(4)
                                       + "  " + date_num.ToString().PadLeft(11) + "");
            }
        }

        public static void easter_gregorian_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EASTER_GREGORIAN_VALUES_TEST tests EASTER_GREGORIAN_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d = 0;
            int m = 0;
            int n_data;
            int y = 0;
            Console.WriteLine("");
            Console.WriteLine("EASTER_GREGORIAN_VALUES_TEST:");
            Console.WriteLine("  EASTER_GREGORIAN_VALUES returns values of ");
            Console.WriteLine("  the date of Easter in the Gregorian calendar.");
            Console.WriteLine("");
            Console.WriteLine("   D   M      Y");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.easter_gregorian_values(ref n_data, ref d, ref m, ref y);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + m.ToString().PadLeft(2)
                                       + "  " + y.ToString().PadLeft(4) + "");
            }
        }

        public static void easter_julian_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EASTER_JULIAN_VALUES_TEST tests EASTER_JULIAN_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 January 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d = 0;
            int m = 0;
            int n_data;
            int y = 0;
            Console.WriteLine("");
            Console.WriteLine("EASTER_JULIAN_VALUES_TEST:");
            Console.WriteLine("  EASTER_JULIAN_VALUES returns values of ");
            Console.WriteLine("  the date of Easter in the Julian calendar.");
            Console.WriteLine("");
            Console.WriteLine("   D   M      Y");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.easter_julian_values(ref n_data, ref d, ref m, ref y);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + d.ToString().PadLeft(2)
                                       + "  " + m.ToString().PadLeft(2)
                                       + "  " + y.ToString().PadLeft(4) + "");
            }
        }
    }
}