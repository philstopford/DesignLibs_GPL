using System;
using Burkardt.Types;

namespace HypercubeGridTest;

using Grid = Burkardt.HyperGeometry.Hypercube.Grid;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERCUBE_GRID_TEST.
        //
        //  Discussion:
        //
        //    HYPERCUBE_GRID_TEST tests HYPERCUBE_GRID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HYPERCUBE_GRID_TEST:");
        Console.WriteLine("  Test the HYPERCUBE_GRID library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("HYPERCUBE_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests HYPERCUBE_GRID on a two dimensional example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 2;

        double[] a = {0.0, 0.0};
        double[] b = {1.0, 10.0};
        int[] c = {2, 4};
        int i;
        int m = M;
        int n;
        int[] ns = {4, 5};
        double[] x;

        n = typeMethods.i4vec_product(m, ns);

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Create a grid using HYPERCUBE_GRID.");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < m; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + ns[i].ToString().PadLeft(6) + "  "
                                                      + c[i].ToString().PadLeft(6) + "  "
                                                      + a[i].ToString().PadLeft(8) + "  "
                                                      + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.hypercube_grid(m, n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(m, n, x, "  Grid points:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests HYPERCUBE_GRID on a five dimensional example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 5;

        double[] a = {0.0, 0.0, 0.0, 0.0, 0.0};
        double[] b = {1.0, 1.0, 1.0, 1.0, 1.0};
        int[] c = {1, 2, 3, 4, 5};
        int i;
        int m = M;
        int n;
        int[] ns = {2, 2, 2, 2, 2};
        double[] x;

        n = typeMethods.i4vec_product(m, ns);

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Create a grid using HYPERCUBE_GRID.");
        Console.WriteLine("  Use a two point grid in each dimension.");
        Console.WriteLine("  Use a different centering option in each dimension.");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < m; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + ns[i].ToString().PadLeft(6) + "  "
                                                      + c[i].ToString().PadLeft(6) + "  "
                                                      + a[i].ToString().PadLeft(8) + "  "
                                                      + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.hypercube_grid(m, n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(m, n, x, "  Grid points:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests HYPERCUBE_GRID on a three dimensional example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 3;

        double[] a = {-1.0, -1.0, -1.0};
        double[] b = {+1.0, +1.0, +1.0};
        int[] c = {1, 1, 1};
        int i;
        int m = M;
        int n;
        int[] ns = {3, 3, 3};
        double[] x;

        n = typeMethods.i4vec_product(m, ns);

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Create a grid using HYPERCUBE_GRID.");
        Console.WriteLine("  Use the same parameters in every dimension.");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < m; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + ns[i].ToString().PadLeft(6) + "  "
                                                      + c[i].ToString().PadLeft(6) + "  "
                                                      + a[i].ToString().PadLeft(8) + "  "
                                                      + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.hypercube_grid(m, n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(m, n, x, "  Grid points:");
    }
}