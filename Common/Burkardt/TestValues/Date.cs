namespace TestValues
{
    public static class Date
    {
        public static void datenum_values(ref int n_data, ref int y, ref int m, ref int d, ref double date_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DATENUM_VALUES returns the MATLAB datenum for various dates.
            //
            //  Discussion:
            //
            //    The MATLAB datenum function returns a numeric value for a given date,
            //    which is 1 for the (fictitious) date 1 January 0.
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
            //    before the first call.  On each call, the routine increments N_DATA by 1,
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int Y, &M, &D, the Common Era date.
            //
            //    Output, ref double DATE_NUM, the datenum.
            //
        {
            int N_MAX = 11;

            int[] d_vec =
            {
                1,
                1,
                1,
                1,
                17,
                9,
                10,
                12,
                6,
                25,
                1
            };

            double[] date_num_vec =
            {
                1.0,
                367.0,
                36526.0,
                365244.0,
                708434.0,
                710284.0,
                713023.0,
                718199.0,
                723186.0,
                729080.0,
                730486.0
            };

            int[] m_vec =
            {
                1,
                1,
                1,
                1,
                8,
                9,
                3,
                5,
                1,
                2,
                1
            };

            int[] y_vec =
            {
                0,
                1,
                100,
                1000,
                1939,
                1944,
                1952,
                1966,
                1980,
                1996,
                2000
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                y = 0;
                m = 0;
                d = 0;
                date_num = 0;
            }
            else
            {
                y = y_vec[n_data - 1];
                m = m_vec[n_data - 1];
                d = d_vec[n_data - 1];
                date_num = date_num_vec[n_data - 1];
            }
        }

        public static void easter_gregorian_values(ref int n_data, ref int d, ref int m, ref int y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EASTER_GREGORIAN_VALUES: dates of Easter in Gregorian calendar.
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int D, &M, &Y, day, month and year of Easter.
            //
        {
            int N_MAX = 10;

            int[] d_vec =
            {
                30, 12, 4, 23, 15, 31, 20, 11, 27, 16
            };

            int[] m_vec =
            {
                3, 4, 4, 4, 4, 3, 4, 4, 3, 4
            };

            int[] y_vec =
            {
                1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                d = 0;
                m = 0;
                y = 0;
            }
            else
            {
                d = d_vec[n_data - 1];
                m = m_vec[n_data - 1];
                y = y_vec[n_data - 1];
            }
        }

        public static void easter_julian_values(ref int n_data, ref int d, ref int m, ref int y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EASTER_JULIAN_VALUES: dates of Easter in Julian calendar.
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int D, &M, &Y, day, month and year of Easter.
            //
        {
            int N_MAX = 10;

            int[] d_vec =
            {
                27, 19, 11, 30, 15, 5, 27, 11, 1, 23
            };

            int[] m_vec =
            {
                4, 4, 4, 4, 4, 5, 4, 4, 5, 4
            };

            int[] y_vec =
            {
                1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                d = 0;
                m = 0;
                y = 0;
            }
            else
            {
                d = d_vec[n_data - 1];
                m = m_vec[n_data - 1];
                y = y_vec[n_data - 1];
            }
        }


        public static void jed_ce_values(ref int n_data, ref double jed, ref int y, ref int m, ref int d,
                ref double f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_CE_VALUES returns the Common Era dates for Julian Ephemeris Dates.
            //
            //  Discussion:
            //
            //    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
            //    starting from noon on 1 January 4713 BCE.
            //
            //    The CE or Common Era is the day, month and year under the
            //    hybrid Julian/Gregorian Calendar, with a transition from Julian
            //    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
            //
            //    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
            //    years BC/BCE are indicated by a negative year value.
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
            //  Reference:
            //
            //    Edward Reingold and Nachum Dershowitz,
            //    Calendrical Calculations: The Millennium Edition,
            //    Cambridge University Press, 2001,
            //    ISBN: 0 521 77752 6
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double JED, the Julian Ephemeris Date.
            //
            //    Output, ref int Y, &M, &D, the Common Era date.
            //
            //    Output, ref double F, the fractional part of the day.
            //
        {
            int N_MAX = 51;

            int[] d_vec =
            {
                1,
                2,
                26,
                8,
                6,
                18,
                8,
                9,
                1,
                26,
                26,
                1,
                1,
                29,
                31,
                1,
                3,
                3,
                29,
                24,
                24,
                29,
                3,
                11,
                12,
                24,
                19,
                15,
                16,
                16,
                21,
                17,
                9,
                4,
                15,
                4,
                13,
                14,
                18,
                22,
                21,
                24,
                17,
                31,
                1,
                6,
                25,
                1,
                9,
                23,
                1
            };
            double[] f_vec =
            {
                0.50,
                0.50,
                0.50,
                0.00,
                0.00,
                0.25,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.50,
                0.50,
                0.00,
                0.50,
                0.50,
                0.00,
                0.00,
                0.00,
                0.00,
                0.00,
                0.81,
                0.00,
                0.00,
                0.00,
                0.00,
                0.33,
                0.00,
                0.50
            };
            double[] jed_vec =
            {
                0.00,
                1.00,
                259261.00,
                347998.50,
                584282.50,
                588465.75,
                758325.50,
                1438178.50,
                1446389.50,
                1448637.50,
                1448637.50,
                1607708.50,
                1607738.50,
                1713262.50,
                1721422.50,
                1721423.50,
                1721425.50,
                1721425.50,
                1724220.50,
                1741959.50,
                1749994.50,
                1825029.50,
                1862836.50,
                1922867.50,
                1936747.50,
                1940351.50,
                1948320.50,
                1948438.50,
                1948439.50,
                1952062.50,
                1952067.50,
                2114872.50,
                2289425.50,
                2299160.00,
                2299161.00,
                2333269.50,
                2361221.00,
                2361222.00,
                2372547.50,
                2375839.50,
                2394646.50,
                2394710.50,
                2400000.50,
                2415020.31,
                2440587.50,
                2444244.50,
                2450138.50,
                2451544.50,
                2453073.83,
                2456284.50,
                2913943.00
            };
            int[] m_vec =
            {
                1,
                1,
                10,
                10,
                9,
                2,
                3,
                7,
                1,
                2,
                2,
                9,
                10,
                8,
                12,
                1,
                1,
                1,
                8,
                3,
                3,
                8,
                3,
                7,
                7,
                5,
                3,
                7,
                7,
                6,
                6,
                3,
                2,
                10,
                10,
                3,
                9,
                9,
                9,
                9,
                3,
                5,
                11,
                12,
                1,
                1,
                2,
                1,
                3,
                12,
                1
            };
            int[] y_vec =
            {
                -4713,
                -4713,
                -4004,
                -3761,
                -3114,
                -3102,
                -2637,
                -776,
                -753,
                -747,
                -747,
                -312,
                -312,
                -23,
                -1,
                1,
                1,
                1,
                8,
                57,
                79,
                284,
                388,
                552,
                590,
                600,
                622,
                622,
                622,
                632,
                632,
                1078,
                1556,
                1582,
                1582,
                1676,
                1752,
                1752,
                1783,
                1792,
                1844,
                1844,
                1858,
                1899,
                1970,
                1980,
                1996,
                2000,
                2004,
                2012,
                3266
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                jed = 0.0;
                y = 0;
                m = 0;
                d = 0;
                f = 0.0;
            }
            else
            {
                jed = jed_vec[n_data - 1];
                y = y_vec[n_data - 1];
                m = m_vec[n_data - 1];
                d = d_vec[n_data - 1];
                f = f_vec[n_data - 1];
            }
        }

        public static void jed_mjd_values(ref int n_data, ref double jed, ref double mjd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_MJD_VALUES returns the MJD for Julian Ephemeris Dates.
            //
            //  Discussion:
            //
            //    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
            //    starting from noon on 1 January 4713 BCE.
            //
            //    The MJD (Modified Julian Day) counts days starting from midnight,
            //    17 November 1858.  This essentially subtracts 2400000.5 days from the JED.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Reingold and Nachum Dershowitz,
            //    Calendrical Calculations: The Millennium Edition,
            //    Cambridge University Press, 2001,
            //    ISBN: 0 521 77752 6
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double JED, the Julian Ephemeris Date.
            //
            //    Output, ref double MJD, the Modified Julian Ephemeris Date.
            //
        {
            int N_MAX = 33;

            double[] jed_vec =
            {
                1507231.5E+00,
                1660037.5E+00,
                1746893.5E+00,
                1770641.5E+00,
                1892731.5E+00,
                1931579.5E+00,
                1974851.5E+00,
                2091164.5E+00,
                2121509.5E+00,
                2155779.5E+00,
                2174029.5E+00,
                2191584.5E+00,
                2195261.5E+00,
                2229274.5E+00,
                2245580.5E+00,
                2266100.5E+00,
                2288542.5E+00,
                2290901.5E+00,
                2323140.5E+00,
                2334848.5E+00,
                2348020.5E+00,
                2366978.5E+00,
                2385648.5E+00,
                2392825.5E+00,
                2416223.5E+00,
                2425848.5E+00,
                2430266.5E+00,
                2430833.5E+00,
                2431004.5E+00,
                2448698.5E+00,
                2450138.5E+00,
                2465737.5E+00,
                2486076.5E+00
            };

            double[] mjd_vec =
            {
                -892769.0E+00,
                -739963.0E+00,
                -653107.0E+00,
                -629359.0E+00,
                -507269.0E+00,
                -468421.0E+00,
                -425149.0E+00,
                -308836.0E+00,
                -278491.0E+00,
                -244221.0E+00,
                -225971.0E+00,
                -208416.0E+00,
                -204739.0E+00,
                -170726.0E+00,
                -154420.0E+00,
                -133900.0E+00,
                -111458.0E+00,
                -109099.0E+00,
                -76860.0E+00,
                -65152.0E+00,
                -51980.0E+00,
                -33022.0E+00,
                -14352.0E+00,
                -7175.0E+00,
                16223.0E+00,
                25848.0E+00,
                30266.0E+00,
                30833.0E+00,
                31004.0E+00,
                48698.0E+00,
                50138.0E+00,
                65737.0E+00,
                86076.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                jed = 0.0;
                mjd = 0.0;
            }
            else
            {
                jed = jed_vec[n_data - 1];
                mjd = mjd_vec[n_data - 1];
            }
        }

        public static void jed_rd_values(ref int n_data, ref double jed, ref double rd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_RD_VALUES returns the RD for Julian Ephemeris Dates.
            //
            //  Discussion:
            //
            //    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
            //    starting from noon on 1 January 4713 BCE.
            //
            //    The RD is the Reingold Dershowitz Date, which counts days from
            //    midnight, 1 January year 1 in the Gregorian calendar.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Reingold and Nachum Dershowitz,
            //    Calendrical Calculations: The Millennium Edition,
            //    Cambridge University Press, 2001,
            //    ISBN: 0 521 77752 6
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double JED, the Julian Ephemeris Date.
            //
            //    Output, ref double RD, the Reingold Dershowitz Date.
            //
        {
            int N_MAX = 33;

            double[] jed_vec =
            {
                1507231.5E+00,
                1660037.5E+00,
                1746893.5E+00,
                1770641.5E+00,
                1892731.5E+00,
                1931579.5E+00,
                1974851.5E+00,
                2091164.5E+00,
                2121509.5E+00,
                2155779.5E+00,
                2174029.5E+00,
                2191584.5E+00,
                2195261.5E+00,
                2229274.5E+00,
                2245580.5E+00,
                2266100.5E+00,
                2288542.5E+00,
                2290901.5E+00,
                2323140.5E+00,
                2334848.5E+00,
                2348020.5E+00,
                2366978.5E+00,
                2385648.5E+00,
                2392825.5E+00,
                2416223.5E+00,
                2425848.5E+00,
                2430266.5E+00,
                2430833.5E+00,
                2431004.5E+00,
                2448698.5E+00,
                2450138.5E+00,
                2465737.5E+00,
                2486076.5E+00
            };

            double[] rd_vec =
            {
                -214193.0E+00,
                -61387.0E+00,
                25469.0E+00,
                49217.0E+00,
                171307.0E+00,
                210155.0E+00,
                253427.0E+00,
                369740.0E+00,
                400085.0E+00,
                434355.0E+00,
                452605.0E+00,
                470160.0E+00,
                473837.0E+00,
                507850.0E+00,
                524156.0E+00,
                544676.0E+00,
                567118.0E+00,
                569477.0E+00,
                601716.0E+00,
                613424.0E+00,
                626596.0E+00,
                645554.0E+00,
                664224.0E+00,
                671401.0E+00,
                694799.0E+00,
                704424.0E+00,
                708842.0E+00,
                709409.0E+00,
                709580.0E+00,
                727274.0E+00,
                728714.0E+00,
                744313.0E+00,
                764652.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                jed = 0.0;
                rd = 0.0;
            }
            else
            {
                jed = jed_vec[n_data - 1];
                rd = rd_vec[n_data - 1];
            }
        }

        public static void jed_weekday_values(ref int n_data, ref double jed, ref int weekday)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JED_WEEKDAY_VALUES returns the day of the week for Julian Ephemeris Dates.
            //
            //  Discussion:
            //
            //    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
            //    starting from noon on 1 January 4713 BCE.
            //
            //    Weekdays are numbered as follows:
            //
            //    1  Sunday
            //    2  Monday
            //    3  Tuesday
            //    4  Wednesday
            //    5  Thursday
            //    6  Friday
            //    7  Saturday
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Reingold and Nachum Dershowitz,
            //    Calendrical Calculations: The Millennium Edition,
            //    Cambridge University Press, 2001,
            //    ISBN: 0 521 77752 6
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double JED, the Julian Ephemeris Date.
            //
            //    Output, ref int WEEKDAY, the day of the week.
            //
        {
            int N_MAX = 33;

            double[] jed_vec =
            {
                1507231.5E+00,
                1660037.5E+00,
                1746893.5E+00,
                1770641.5E+00,
                1892731.5E+00,
                1931579.5E+00,
                1974851.5E+00,
                2091164.5E+00,
                2121509.5E+00,
                2155779.5E+00,
                2174029.5E+00,
                2191584.5E+00,
                2195261.5E+00,
                2229274.5E+00,
                2245580.5E+00,
                2266100.5E+00,
                2288542.5E+00,
                2290901.5E+00,
                2323140.5E+00,
                2334848.5E+00,
                2348020.5E+00,
                2366978.5E+00,
                2385648.5E+00,
                2392825.5E+00,
                2416223.5E+00,
                2425848.5E+00,
                2430266.5E+00,
                2430833.5E+00,
                2431004.5E+00,
                2448698.5E+00,
                2450138.5E+00,
                2465737.5E+00,
                2486076.5E+00
            };

            int[] weekday_vec =
            {
                1, 4, 4, 1, 4,
                2, 7, 1, 1, 6,
                7, 6, 1, 1, 4,
                7, 7, 7, 4, 1,
                6, 1, 2, 4, 1,
                1, 2, 2, 5, 3,
                1, 4, 1
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                jed = 0.0;
                weekday = 0;
            }
            else
            {
                jed = jed_vec[n_data - 1];
                weekday = weekday_vec[n_data - 1];
            }
        }

        public static void weekday_values(ref int n_data, ref int y, ref int m, ref int d, ref int w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WEEKDAY_VALUES returns the day of the week for various dates.
            //
            //  Discussion:
            //
            //    The CE or Common Era calendar is used, under the
            //    hybrid Julian/Gregorian Calendar, with a transition from Julian
            //    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
            //
            //    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
            //    years BC/BCE are indicated by a negative year value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 May 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Reingold, Nachum Dershowitz,
            //    Calendrical Calculations: The Millennium Edition,
            //    Cambridge University Press, 2001,
            //    ISBN: 0 521 77752 6
            //    LC: CE12.R45.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
            //    before the first call.  On each call, the routine increments N_DATA by 1,
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int Y, &M, &D, the Common Era date.
            //
            //    Output, ref int W, the day of the week.  Sunday = 1.
            //
        {
            int N_MAX = 34;

            int[] d_vec =
            {
                30,
                8,
                26,
                3,
                7,
                18,
                7,
                19,
                14,
                18,
                16,
                3,
                26,
                20,
                4,
                25,
                31,
                9,
                24,
                10,
                30,
                24,
                19,
                2,
                27,
                19,
                25,
                29,
                19,
                7,
                17,
                25,
                10,
                18
            };
            int[] m_vec =
            {
                7,
                12,
                9,
                10,
                1,
                5,
                11,
                4,
                10,
                5,
                3,
                3,
                3,
                4,
                6,
                1,
                3,
                9,
                2,
                6,
                6,
                7,
                6,
                8,
                3,
                4,
                8,
                9,
                4,
                10,
                3,
                2,
                11,
                7
            };
            int[] w_vec =
            {
                1,
                4,
                4,
                1,
                4,
                2,
                7,
                1,
                7,
                1,
                6,
                7,
                6,
                1,
                1,
                4,
                7,
                7,
                7,
                4,
                1,
                6,
                1,
                2,
                4,
                1,
                1,
                2,
                2,
                5,
                3,
                1,
                4,
                1
            };
            int[] y_vec =
            {
                -587,
                -169,
                70,
                135,
                470,
                576,
                694,
                1013,
                1066,
                1096,
                1190,
                1240,
                1288,
                1298,
                1391,
                1436,
                1492,
                1553,
                1560,
                1648,
                1680,
                1716,
                1768,
                1819,
                1839,
                1903,
                1929,
                1941,
                1943,
                1943,
                1992,
                1996,
                2038,
                2094
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                y = 0;
                m = 0;
                d = 0;
                w = 0;
            }
            else
            {
                y = y_vec[n_data - 1];
                m = m_vec[n_data - 1];
                d = d_vec[n_data - 1];
                w = w_vec[n_data - 1];
            }
        }


    }
}