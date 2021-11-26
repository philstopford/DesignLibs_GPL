using System;
using Burkardt.Types;

namespace EllipseGridTest;

using Grid = Burkardt.Ellipse.Grid;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ELLIPSE_GRID_TEST.
        //
        //  Discussion:
        //
        //    ELLIPSE_GRID_TEST tests the ELLIPSE_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_GRID_TEST:");
        Console.WriteLine("  Test the ELLIPSE_GRID library.");

        ellipse_grid_test01();

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void ellipse_grid_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_GRID_TEST01 tests ELLIPSE_GRID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[2];
        const string filename = "ellipse_grid_test01.xy";
        double[] r = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  ELLIPSE_GRID can define a grid of points");
        Console.WriteLine("  with N+1 points on the minor half axis,");
        Console.WriteLine("  based on any ellipse.");

        const int n = 8;
        r[0] = 2.0;
        r[1] = 1.0;
        c[0] = 1.0;
        c[1] = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  We use N = " + n + "");
        Console.WriteLine("  Radius R = (" + r[0] + "," + r[1] + ")");
        Console.WriteLine("  Center C = (" + c[0] + "," + c[1] + ")");

        int ng = Grid.ellipse_grid_count(n, r, c);

        Console.WriteLine("");
        Console.WriteLine("  Number of grid points will be " + ng + "");

        double[] xy = Grid.ellipse_grid(n, r, c, ng);

        typeMethods.r82vec_print_part(ng, xy, 20, "  Part of the grid point array:");

        typeMethods.r8mat_write(filename, 2, ng, xy);

        Console.WriteLine("");
        Console.WriteLine("  Data written to the file \"" + filename + "\".");
    }
}