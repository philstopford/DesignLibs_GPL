using System;
using Burkardt.Ball;
using Burkardt.Types;

namespace BallGridTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BALL_GRID_TEST.
        //
        //  Discussion:
        //
        //    BALL_GRID_TEST tests the BALL_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BALL_GRID_TEST:");
        Console.WriteLine("  Test the BALL_GRID library.");

        ball_grid_test01();

        Console.WriteLine("");
        Console.WriteLine("BALL_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void ball_grid_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL_GRID_TEST01 tests BALL_GRID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] bg;
        double[] c = new double[3];
        string filename = "ball_grid_test01.xyz";
        int n;
        int ng;
        double r;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  BALL_GRID can define a grid of points");
        Console.WriteLine("  with N+1 points on a horizontal or vertical radius,");
        Console.WriteLine("  based on any ball.");

        n = 10;
        r = 2.0;
        c[0] = 1.0;
        c[1] = 5.0;
        c[2] = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  We use N = " + n + "");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center C = (" + c[0] + "," + c[1] + "," + c[2] + ")");

        ng = Grid.ball_grid_count(n, r, c);

        Console.WriteLine("");
        Console.WriteLine("  Number of grid points will be " + ng + "");

        bg = Grid.ball_grid(n, r, c, ng);

        typeMethods.r83vec_print_part(ng, bg, 20, "  Part of the grid point array:");

        typeMethods.r8mat_write(filename, 3, ng, bg);

        Console.WriteLine("");
        Console.WriteLine("  Data written to the file \"" + filename + "\".");

    }
}