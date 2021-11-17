using System;
using Burkardt;
using Burkardt.Function;

namespace HilbertTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HILBERT_CURVE_TEST.
        //
        //  Discussion:
        //
        //    HILBERT_CURVE_TEST tests HILBERT_CURVE.
        //
        //  Modified:
        //
        //    02 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HILBERT_CURVE_TEST:");
        Console.WriteLine("  Test the HILBERT_CURVE library.");

        d2xy_test();
        rot_test();
        xy2d_test();

        Console.WriteLine("");
        Console.WriteLine("HILBERT_CURVE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void d2xy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    D2XY_TEST tests D2XY.
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int m;
        int n;
        int x = 0;
        int y = 0;

        Console.WriteLine("");
        Console.WriteLine("D2XY_TEST:");
        Console.WriteLine("  D2XY converts a Hilbert linear D coordinate to an (X,Y) 2D coordinate.");

        m = 3;
        n = (int) Math.Pow(2, m);

        Console.WriteLine("");
        Console.WriteLine("    D    X    Y");
        Console.WriteLine("");
        for (d = 0; d < n * n; d++)
        {
            Hilbert.d2xy(m, d, ref x, ref y);
            Console.WriteLine("  " + d.ToString().PadLeft(3)
                                   + "  " + x.ToString().PadLeft(3)
                                   + "  " + y.ToString().PadLeft(3) + "");
        }
    }

    private static void rot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROT_TEST tests ROT.
        //
        //  Modified:
        //
        //    02 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        int rx;
        int ry;
        int x;
        int x0;
        int x1;
        int y;
        int y0;
        int y1;

        Console.WriteLine("");
        Console.WriteLine("ROT_TEST:");
        Console.WriteLine("  ROT rotates and flips a quadrant appropriately.");
        Console.WriteLine("");
        Console.WriteLine("   X   Y  X0  Y0  X1  Y1");
        Console.WriteLine("");

        m = 3;
        n = (int) Math.Pow(2, m);
        ry = 0;

        for (y = 0; y < 8; y++)
        {
            for (x = 0; x < 8; x++)
            {
                rx = 0;
                x0 = x;
                y0 = y;
                Hilbert.rot(n, ref x0, ref y0, rx, ry);
                rx = 1;
                x1 = x;
                y1 = y;
                Hilbert.rot(n, ref x1, ref y1, rx, ry);
                Console.WriteLine("  " + x.ToString().PadLeft(2)
                                       + "  " + y.ToString().PadLeft(2)
                                       + "  " + x0.ToString().PadLeft(2)
                                       + "  " + y0.ToString().PadLeft(2)
                                       + "  " + x1.ToString().PadLeft(2)
                                       + "  " + y1.ToString().PadLeft(2) + "");
            }
        }
    }

    private static void xy2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY2D_TEST tests XY2D.
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int m;
        int n;
        int x;
        int y;

        Console.WriteLine("");
        Console.WriteLine("XY2D_TEST:");
        Console.WriteLine("  XY2D converts an (X,Y) 2D coordinate to a Hilbert linear D coordinate.");

        m = 3;
        n = (int) Math.Pow(2, m);

        Console.WriteLine("");
        string cout = "        ";
        for (x = 0; x < n; x++)
        {
            cout += x.ToString().PadLeft(3);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        for (y = n - 1; 0 <= y; y--)
        {
            cout = "  " + y.ToString().PadLeft(3) + ":  ";
            for (x = 0; x < n; x++)
            {
                d = Hilbert.xy2d(m, x, y);
                cout += d.ToString().PadLeft(3);
            }

            Console.WriteLine(cout);
        }
    }
}