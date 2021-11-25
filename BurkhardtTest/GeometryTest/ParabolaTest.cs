using System;
using Burkardt.ParabolaNS;

namespace GeometryTest;

public static class ParabolaTest
{
    public static void test0493()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0493 tests PARABOLA_EX and PARABOLA_EX2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double xmin = 0;
        double ymin = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST0493");
        Console.WriteLine("  PARABOLA_EX finds the extreme value of a parabola");
        Console.WriteLine("  determined by three points.");
        Console.WriteLine("  PARABOLA_EX2 finds the extreme value of a parabola");
        Console.WriteLine("  determined by three points.");

        double a = 2.0;
        double b = -4.0;
        double c = 10.0;

        const double x1 = 1.0;
        double y1 = a * x1 * x1 + b * x1 + c;
        const double x2 = 2.0;
        double y2 = a * x2 * x2 + b * x2 + c;
        const double x3 = 3.0;
        double y3 = a * x3 * x3 + b * x3 + c;

        Console.WriteLine("");
        Console.WriteLine("  Parabolic coefficients (A,B,C) = "
                          + a + "  " + b + "  " + c + "");
        Console.WriteLine("");
        Console.WriteLine("  X, Y data:");
        Console.WriteLine("");
        Console.WriteLine("  X1, Y1 = " + x1 + "  " + y1 + "");
        Console.WriteLine("  X2, Y2 = " + x2 + "  " + y2 + "");
        Console.WriteLine("  X3, Y3 = " + x3 + "  " + y3 + "");

        a = 0.0;
        b = 0.0;
        c = 0.0;

        int ierror = Geometry.parabola_ex(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

        switch (ierror)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("  PARABOLA_EX returns (XMIN,YMIN) = "
                                  + xmin + "  " + ymin + "");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  PARABOLA_EX returns error code " + ierror + "");
                break;
        }

        ierror = Geometry.parabola_ex2(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin, ref a, ref b, ref c);

        switch (ierror)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("  PARABOLA_EX2 returns (XMIN,YMIN) = " + xmin
                                                                          + "  " + ymin + "");
                Console.WriteLine("  and (A,B,C) = " + a + "  " + b + "  " + c + "");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  PARABOLA_EX2 returns error code " + ierror + "");
                break;
        }

    }
}