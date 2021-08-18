using System;
using Burkardt.Values;

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

        public static void jed_ce_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_CE_VALUES_TEST tests JED_CE_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 May 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d = 0;
            double f = 0;
            double jed = 0;
            int n_data;
            int m = 0;
            int y = 0;
            Console.WriteLine("");
            Console.WriteLine("JED_CE_VALUES_TEST:");
            Console.WriteLine("  JED_CE_VALUES returns:");
            Console.WriteLine("  JED, ref a Julian Ephemeris Date, ref and");
            Console.WriteLine("  YMDF, the corresponding year, ref month, ref day, ref fraction.");
            Console.WriteLine("");
            Console.WriteLine("        JED          Y   M   D    F");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.jed_ce_values(ref n_data, ref jed, ref y, ref m, ref d, ref f);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + jed.ToString().PadLeft(12)
                                       + "  " + y.ToString().PadLeft(6)
                                       + "  " + m.ToString().PadLeft(2)
                                       + "  " + d.ToString().PadLeft(2)
                                       + "  " + f.ToString().PadLeft(6) + "");
            }
        }

        public static void jed_mjd_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_MJD_VALUES_TEST tests JED_MJD_VALUES.
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
            double jed = 0;
            int n_data;
            double mjd = 0;
            Console.WriteLine("");
            Console.WriteLine("JED_MJD_VALUES_TEST:");
            Console.WriteLine("  JED_MJD_VALUES returns:");
            Console.WriteLine("  JED, ref a Julian Ephemeris Date, ref and");
            Console.WriteLine("  MJD, the corresponding Modified Julian Day count.");
            Console.WriteLine("");
            Console.WriteLine("   JED      MJD");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.jed_mjd_values(ref n_data, ref jed, ref mjd);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + jed.ToString().PadLeft(12) + "  "
                                  + mjd.ToString().PadLeft(12) + "");
            }
        }

        public static void jed_rd_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_RD_VALUES_TEST tests JED_RD_VALUES.
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
            double jed = 0;
            int n_data;
            double rd = 0;
            Console.WriteLine("");
            Console.WriteLine("JED_RD_VALUES_TEST:");
            Console.WriteLine("  JED_RD_VALUES returns:");
            Console.WriteLine("  JED, ref a Julian Ephemeris Date, ref and");
            Console.WriteLine("  RD, the corresponding Reingold Dershowitz Day count.");
            Console.WriteLine("");
            Console.WriteLine("   JED      RD");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.jed_rd_values(ref n_data, ref jed, ref rd);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + jed.ToString().PadLeft(12) + "  "
                                  + rd.ToString().PadLeft(12) + "");
            }
        }

        public static void jed_weekday_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_WEEKDAY_VALUES_TEST tests JED_WEEKDAY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 September 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double jed = 0;
            int n_data;
            int weekday = 0;
            string[] weekday_name =
                {
                    "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday",
                    "Friday", "Saturday"
                }
                ;
            Console.WriteLine("");
            Console.WriteLine("JED_WEEKDAY_VALUES_TEST:");
            Console.WriteLine("  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates ");
            Console.WriteLine("  (JED) and the corresponding weekday");
            Console.WriteLine("");
            Console.WriteLine("   JED      #  Weekday");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.jed_weekday_values(ref n_data, ref jed, ref weekday);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + jed.ToString().PadLeft(12) + "  "
                                  + weekday.ToString().PadLeft(1) + "  "
                                  + weekday_name[weekday - 1] + "");
            }
        }
        
        public static void weekday_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEEKDAY_VALUES_TEST tests WEEKDAY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d = 0;
            int m = 0;
            int n_data;
            int w = 0; int y = 0;
            Console.WriteLine("");
            Console.WriteLine("WEEKDAY_VALUES_TEST:");
            Console.WriteLine("  WEEKDAY_VALUES returns values of ");
            Console.WriteLine("  the weekday for a given Y/M/D date.");
            Console.WriteLine("");
            Console.WriteLine("     Y     M     D     W");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Date.weekday_values(ref n_data, ref y, ref m, ref d, ref w);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + y.ToString().PadLeft(4)
                                       + "  " + m.ToString().PadLeft(4)
                                       + "  " + d.ToString().PadLeft(4)
                                       + "  " + w.ToString().PadLeft(4) + "");
            }
        }
    }
}