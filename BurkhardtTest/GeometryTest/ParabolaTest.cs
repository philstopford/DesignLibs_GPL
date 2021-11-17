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
        double a;
        double b;
        double c;
        int ierror = 0;
        double x1;
        double x2;
        double x3;
        double xmin = 0;
        double y1;
        double y2;
        double y3;
        double ymin = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST0493");
        Console.WriteLine("  PARABOLA_EX finds the extreme value of a parabola");
        Console.WriteLine("  determined by three points.");
        Console.WriteLine("  PARABOLA_EX2 finds the extreme value of a parabola");
        Console.WriteLine("  determined by three points.");

        a = 2.0;
        b = -4.0;
        c = 10.0;

        x1 = 1.0;
        y1 = a * x1 * x1 + b * x1 + c;
        x2 = 2.0;
        y2 = a * x2 * x2 + b * x2 + c;
        x3 = 3.0;
        y3 = a * x3 * x3 + b * x3 + c;

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

        ierror = Geometry.parabola_ex(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

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