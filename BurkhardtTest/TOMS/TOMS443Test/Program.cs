﻿using System;
using System.Globalization;
using Burkardt.WFunction;

namespace TOMS443Test;

internal static class Program
{
    private static void Main()
            
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS443_TEST.
        //
        //  Discussion:
        //
        //    TOMS443_TEST tests the TOMS443 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOMS443_TEST");
        Console.WriteLine("  Test the TOMS443 library.");

        test01();
        test02();
        Console.WriteLine("");
        Console.WriteLine("TOMS433_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests WEW_A
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double en = 0;
        double w1 = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test WEW_A to evaluate");
        Console.WriteLine("  Lambert's W function.");
        Console.WriteLine("");
        Console.WriteLine("          X             Exact             Computed      Error");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Lambert.lambert_w_values(ref n_data, ref x, ref w1);

            if (n_data <= 0)
            {
                break;
            }

            double w2 = x switch
            {
                0.0 => 0.0,
                _ => Lambert.wew_a(x, ref en)
            };

            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                       + w1.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "  "
                                                       + w2.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "  "
                                                       + Math.Abs(w1 - w2).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests WEW_B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double en = 0;
        double w1 = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test WEW_B to evaluate");
        Console.WriteLine("  Lambert's W function.");
        Console.WriteLine("");
        Console.WriteLine("          X             Exact             Computed      Error");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Lambert.lambert_w_values(ref n_data, ref x, ref w1);

            if (n_data <= 0)
            {
                break;
            }

            double w2 = x switch
            {
                0.0 => 0.0,
                _ => Lambert.wew_b(x, ref en)
            };

            Console.WriteLine(x.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                       + w1.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "  "
                                                       + w2.ToString(CultureInfo.InvariantCulture).PadLeft(16) + "  "
                                                       + Math.Abs(w1 - w2).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }
}

internal class TestValues
{
}