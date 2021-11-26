using System;
using System.Globalization;
using Burkardt.Types;

namespace LineGridTest;

using Grid = Burkardt.LineNS.Grid;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LINE_GRID_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LINE_GRID_TEST:");
        Console.WriteLine("  Test the LINE_GRID library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("LINE_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests LINE_GRID using simple parameters.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double a = -1.0;
        const double b = +1.0;
        const int c = 1;
        const int n = 11;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Create a grid using LINE_GRID.");
        Console.WriteLine("  Use simple parameters.");
        Console.WriteLine("");
        Console.WriteLine("     N     C      A         B");
        Console.WriteLine("");
        Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                  + c.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                  + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                  + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        double[] x = Grid.line_grid(n, a, b, c);
        typeMethods.r8vec_print(n, x, "  Grid points:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 changes the number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double a = 0.0;
        const double b = 1.0;
        const int c = 2;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Create a grid using LINE_GRID.");
        Console.WriteLine("  Try an increasing number of points.");

        int n = 4;

        for (test = 1; test <= 3; test++)
        {
            n = 2 * n + 1;

            Console.WriteLine("");
            Console.WriteLine("     N     C      A         B");
            Console.WriteLine("");
            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + c.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                      + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            double[] x = Grid.line_grid(n, a, b, c);
            typeMethods.r8vec_print(n, x, "  Grid points:");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tries all the centering options.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const double a = 0.0;
        const double b = 100.0;
        int c;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Create a grid using LINE_GRID.");
        Console.WriteLine("  Try the different centering options.");

        const int n = 5;

        for (c = 1; c <= 5; c++)
        {
            Console.WriteLine("");
            Console.WriteLine("     N     C      A         B");
            Console.WriteLine("");
            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + c.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                      + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            double[] x = Grid.line_grid(n, a, b, c);
            typeMethods.r8vec_print(n, x, "  Grid points:");
        }
    }
}