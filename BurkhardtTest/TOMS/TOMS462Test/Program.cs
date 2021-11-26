﻿using System;
using System.Globalization;
using Burkardt.CDFLib;
using Burkardt.Values;

namespace TOMS462Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS462_TEST.
        //
        //  Discussion:
        //
        //    TOMS462_TEST tests the TOMS462 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOMS462_TEST");
        Console.WriteLine("  Test the TOMS462 library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("TOMS462_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************807
        //
        //  Purpose:
        //
        //    TEST01 tests BIVNOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bivariatenormal.BivnorData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Compare BIVNOR with some simple data");
        Console.WriteLine("  with 3 digit accuracy.");
        Console.WriteLine("");
        Console.WriteLine("       X         Y          R          P               P");
        Console.WriteLine("                                      (Tabulated)     (BIVNOR)");
        Console.WriteLine("");

        double x = 0.8;
        double y = -1.5;
        double r = -0.9;
        double expect = 0.148;
        double cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        x = 0.6;
        y = -1.4;
        r = -0.7;
        expect = 0.208;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        x = 0.2;
        y = -1.0;
        r = -0.5;
        expect = 0.304;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        x = -1.2;
        y = 0.1;
        r = 0.0;
        expect = 0.407;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        x = -1.2;
        y = -0.1;
        r = 0.3;
        expect = 0.501;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        x = -0.4;
        y = -0.9;
        r = 0.6;
        expect = 0.601;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                               + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + expect.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                               + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BIVNOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fxy1 = 0;
        double r = 0;
        double x = 0;
        double y = 0;
        bivariatenormal.BivnorData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Compare BIVNOR with some tabulated data.");
        Console.WriteLine("");
        Console.WriteLine("      X          Y          " +
                          "R           P                         P" +
                          "                      DIFF" +
                          "                                " +
                          "       (Tabulated)               (BIVNOR)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Bivariate.bivariate_normal_cdf_values(ref n_data, ref x, ref y, ref r, ref fxy1);

            if (n_data == 0)
            {
                break;
            }

            //
            //  BIVNOR computes the "tail" of the probability, and we want the
            //  initial part//
            //
            double fxy2 = bivariatenormal.bivnor(ref data, -x, -y, r);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fxy1.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + fxy2.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + Math.Abs(fxy1 - fxy2).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }
}