using System;
using Burkardt.CDFLib;
using Burkardt.Values;

namespace TOMS462Test;

internal class Program
{
    private static void Main(string[] args)
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
        double cdf;
        double expect;
        double r;
        double x;
        double y;
        bivariatenormal.BivnorData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Compare BIVNOR with some simple data");
        Console.WriteLine("  with 3 digit accuracy.");
        Console.WriteLine("");
        Console.WriteLine("       X         Y          R          P               P");
        Console.WriteLine("                                      (Tabulated)     (BIVNOR)");
        Console.WriteLine("");

        x = 0.8;
        y = -1.5;
        r = -0.9;
        expect = 0.148;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

        x = 0.6;
        y = -1.4;
        r = -0.7;
        expect = 0.208;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

        x = 0.2;
        y = -1.0;
        r = -0.5;
        expect = 0.304;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

        x = -1.2;
        y = 0.1;
        r = 0.0;
        expect = 0.407;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

        x = -1.2;
        y = -0.1;
        r = 0.3;
        expect = 0.501;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

        x = -0.4;
        y = -0.9;
        r = 0.6;
        expect = 0.601;
        cdf = bivariatenormal.bivnor(ref data, x, y, r);
        Console.WriteLine("  " + x.ToString().PadLeft(9)
                               + "  " + y.ToString().PadLeft(8)
                               + "  " + r.ToString().PadLeft(8)
                               + "  " + expect.ToString().PadLeft(14)
                               + "  " + cdf.ToString().PadLeft(14) + "");

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
        double fxy2 = 0;
        int n_data = 0;
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

        n_data = 0;

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
            fxy2 = bivariatenormal.bivnor(ref data, -x, -y, r);

            Console.WriteLine("  " + x.ToString().PadLeft(8)
                                   + "  " + y.ToString().PadLeft(8)
                                   + "  " + r.ToString().PadLeft(8)
                                   + "  " + fxy1.ToString().PadLeft(24)
                                   + "  " + fxy2.ToString().PadLeft(24)
                                   + "  " + Math.Abs(fxy1 - fxy2).ToString().PadLeft(10) + "");
        }
    }
}