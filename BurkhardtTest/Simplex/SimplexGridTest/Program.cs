using System;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace SimplexGridTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_GRID_TEST tests the SIMPLEX_GRID library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_GRID_TEST:");
        Console.WriteLine("  Test the SIMPLEX_GRID library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SIMPLEX_GRID_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        int ng;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  SIMPLEX_GRID_SIZE counts the points in a regular grid");
        Console.WriteLine("  with N+1 points on a side, in an M-dimensional simplex.");
        Console.WriteLine("");
        Console.WriteLine("        M: 0     1     2     3     4     5");
        Console.WriteLine("    N:");
        for (n = 0; n <= 10; n++)
        {
            string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(3) + ":";
            for (m = 0; m <= 5; m++)
            {
                ng = Grid.simplex_grid_size(m, n);
                cout += ng.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests SIMPLEX_GRID_INDEX_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] g;
        int i;
        int j;
        int m = 3;
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  SIMPLEX_GRID_INDEX_NEXT lists, one by one, the indices");
        Console.WriteLine("  of a simplex grid that uses N+1 points on a side,");
        Console.WriteLine("  in an M-dimensional simplex.");
        Console.WriteLine("");
        Console.WriteLine("   #:  1  2  3  (*)");
        Console.WriteLine("");

        n = 3;

        j = 0;
        g = new int[m + 1];
        for (i = 0; i < m; i++)
        {
            g[i] = 0;
        }

        g[m] = n;

        while (true)
        {
            string cout = "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2);
            for (i = 0; i < m; i++)
            {
                cout += g[i].ToString(CultureInfo.InvariantCulture).PadLeft(3);
            }

            Console.WriteLine(cout + " (" + g[m].ToString(CultureInfo.InvariantCulture).PadLeft(3) + ")");

            if (g[0] == n)
            {
                break;
            }

            Grid.simplex_grid_index_next(m, n, g);

            j += 1;
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests SIMPLEX_GRID_INDEX_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] g;
        int i;
        int j;
        int m = 3;
        int n;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  SIMPLEX_GRID_INDEX_SAMPLE returns a randomly selected");
        Console.WriteLine("  index of a simplex grid that uses N+1 points on a side,");
        Console.WriteLine("  in an M-dimensional simplex.");
        Console.WriteLine("");
        Console.WriteLine("   #:  1  2  3  (*)");
        Console.WriteLine("");

        n = 3;
        seed = 123456789;

        for (j = 1; j <= 20; j++)
        {
            g = Grid.simplex_grid_index_sample(m, n, ref seed);

            string cout = "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2);
            for (i = 0; i < m; i++)
            {
                cout += g[i].ToString(CultureInfo.InvariantCulture).PadLeft(3);
            }

            Console.WriteLine(cout + " (" + g[m].ToString(CultureInfo.InvariantCulture).PadLeft(3) + ")");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests SIMPLEX_GRID_INDEX_TO_POINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] g;
        int i;
        int j;
        int m = 2;
        int n;
        int seed;
        double[] v =
        {
            20.0, 0.0,
            30.0, 40.0,
            10.0, 20.0
        };
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point");
        Console.WriteLine("  corresponding to a grid index of a simplex grid that");
        Console.WriteLine("  that uses N+1 points on a side,");
        Console.WriteLine("  in an M-dimensional simplex.");

        n = 5;

        typeMethods.r8mat_transpose_print(m, m + 1, v, "  Simplex vertices:");

        Console.WriteLine("");
        Console.WriteLine("  Choosing random simplex indices to convert:");
        Console.WriteLine("   #:  1  2  3     X        Y");
        Console.WriteLine("");

        seed = 123456789;

        for (j = 1; j <= 10; j++)
        {
            g = Grid.simplex_grid_index_sample(m, n, ref seed);
            x = Grid.simplex_grid_index_to_point(m, n, 1, g, v);

            string cout = "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2) + ":";
            for (i = 0; i <= m; i++)
            {
                cout += g[i].ToString(CultureInfo.InvariantCulture).PadLeft(3);
            }

            cout += "  ";
            for (i = 0; i < m; i++)
            {
                cout += x[i].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout);

        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests SIMPLEX_GRID_INDEX_ALL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] grid;
        int m;
        int n;
        int ng;

        Console.WriteLine("");
        Console.WriteLine("TEST05:");
        Console.WriteLine("  SIMPLEX_GRID_INDEX_ALL returns all the indices");
        Console.WriteLine("  of a simplex grid that uses N+1 points on a side,");
        Console.WriteLine("  in an M-dimensional simplex.");

        m = 3;
        n = 3;
        ng = Grid.simplex_grid_size(m, n);

        grid = Grid.simplex_grid_index_all(m, n, ng);

        typeMethods.i4mat_transpose_print(m + 1, ng, grid,
            "  Transposed Simplex Grid Index Matrix:");

    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests SIMPLEX_GRID_INDEX_TO_POINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] grid;
        int m = 2;
        int n;
        int ng;
        double[] v =
        {
            20.0, 0.0,
            30.0, 40.0,
            10.0, 20.0
        };
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST06:");
        Console.WriteLine("  SIMPLEX_GRID_INDEX_TO_POINT returns the physical point");
        Console.WriteLine("  corresponding to a grid index of a simplex grid that");
        Console.WriteLine("  that uses N+1 points on a side,");
        Console.WriteLine("  in an M-dimensional simplex.");

        n = 5;
        ng = Grid.simplex_grid_size(m, n);

        typeMethods.r8mat_transpose_print(m, m + 1, v, "  Simplex vertices:");

        grid = Grid.simplex_grid_index_all(m, n, ng);

        x = Grid.simplex_grid_index_to_point(m, n, ng, grid, v);

        typeMethods.r8mat_transpose_print(m, ng, x, "  Grid Point Coordinates:");
    }
}