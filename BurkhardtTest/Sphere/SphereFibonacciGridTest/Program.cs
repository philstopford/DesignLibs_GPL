using System;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereFibonacciGridTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_FIBONACCI_GRID_TEST tests the SPHERE_FIBONACCI_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_FIBONACCI_GRID_TEST");
        Console.WriteLine("  Test the SPHERE_FIBONACCI_GRID library.");

        sphere_fibonacci_grid_points_test();
        sphere_fibonacci_grid_display_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SPHERE_FIBONACCI_GRID_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    public static void sphere_fibonacci_grid_points_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_FIBONACCI_GRID_POINTS_TEST tests SPHERE_FIBONACCI_GRID_POINTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_FIBONACCI_GRID_POINTS_TEST");
        Console.WriteLine("  SPHERE_FIBONACCI_GRID_POINTS computes points on a sphere");
        Console.WriteLine("  that lie on a Fibonacci spiral.");

        const int ng = 1000;
        Console.WriteLine("");
        Console.WriteLine("  Number of points NG = " + ng + "");

        double[] xg = Grid_Fibonacci.sphere_fibonacci_grid_points(ng);

        typeMethods.r8mat_transpose_print_some(3, ng, xg, 1, 1, 3, 10,
            "  Part of the grid array:");
        //
        //  Write the nodes to a file.
        //
        string filename = "sphere_fibonacci_grid_n1000.xyz";

        typeMethods.r8mat_write(filename, 3, ng, xg);
    }

    public static void sphere_fibonacci_grid_display_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_FIBONACCI_GRID_DISPLAY_TEST tests SPHERE_FIBONACCI_GRID_DISPLAY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_FIBONACCI_GRID_DISPLAY_TEST");
        Console.WriteLine("  SPHERE_FIBONACCI_GRID_DISPLAY displays points");
        Console.WriteLine("  on a sphere that lie on a Fibonacci spiral.");

        const int ng = 1000;
        Console.WriteLine("");
        Console.WriteLine("  Number of points NG = " + ng + "");

        double[] xg = Grid_Fibonacci.sphere_fibonacci_grid_points(ng);
        //
        //  Display the nodes on a sphere.
        //
        string prefix = "sphere_fibonacci_grid_n1000";

        Grid_Fibonacci.sphere_fibonacci_grid_display(ng, xg, prefix);
    }
}