using System;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetrahedronGridTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TETRAHEDRON_GRID_TEST.
        //
        //  Discussion:
        //
        //    TETRAHEDRON_GRID_TEST tests the TETRAHEDRON_GRID library.
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
        Console.WriteLine("TETRAHEDRON_GRID_TEST:");
        Console.WriteLine("  Test the TETRAHEDRON_GRID library.");

        tetrahedron_grid_test01();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void tetrahedron_grid_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_GRID_TEST01 tests TETRAHEDRON_GRID.
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
        string filename = "tetrahedron_grid_test01.xyz";
        int n;
        int ng;
        double[] t =
        {
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0
        };
        double[] tg;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  TETRAHEDRON_GRID can define a tetrahedral grid");
        Console.WriteLine("  with N+1 points on a side, based on any tetrahedron.");

        n = 10;
        Console.WriteLine("  N = " + n + "");

        ng = Grid.tetrahedron_grid_count(n);

        typeMethods.r8mat_print(3, 4, t, "  Tetrahedron vertices:");

        tg = Grid.tetrahedron_grid(n, t, ng);

        typeMethods.r83vec_print_part(ng, tg, 20, "  Part of the grid point array:");

        typeMethods.r8mat_write(filename, 3, ng, tg);

        Console.WriteLine("");
        Console.WriteLine("  Data written to the file \"" + filename + "\".");
    }
}