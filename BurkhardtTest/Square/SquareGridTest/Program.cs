using System;
using Burkardt.Square;
using Burkardt.Types;

namespace SquareGridTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SQUARE_GRID_TEST.
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
        Console.WriteLine("");
        Console.WriteLine("SQUARE_GRID_TEST:");
            
        Console.WriteLine("  Test the SQUARE_GRID library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("SQUARE_GRID_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests SQUARE_GRID using the same parameters for all dimensions.
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
        double[] a = { -1.0, -1.0 };
        double[] b = { +1.0, +1.0 };
        int[] c = { 1, 1 };
        int i;
        int n;
        int[] ns = { 3, 3 };
        double[] x;

        n = ns[0] * ns[1];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Create a grid using SQUARE_GRID.");
        Console.WriteLine("  Use the same parameters in every dimension.");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < 2; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + ns[i].ToString().PadLeft(4) + "  "
                                                      + c[i].ToString().PadLeft(4) + "  "
                                                      + a[i].ToString().PadLeft(8) + "  "
                                                      + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.square_grid(n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(2, n, x, "  Grid points:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses a different number of points in each coordinate.
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
        double[] a = { 0.0, 0.0 };
        double[] b = { 1.0, 1.0 };
        int[] c = { 2, 2 };
        int i;
        int n;
        int[] ns = { 4, 2 };
        double[] x;

        n = ns[0] * ns[1];

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Create a grid using SQUARE_GRID.");
        Console.WriteLine("  se a different number of points in each dimension..");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < 2; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + ns[i].ToString().PadLeft(4) + "  "
                                                      + c[i].ToString().PadLeft(4) + "  "
                                                      + a[i].ToString().PadLeft(8) + "  "
                                                      + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.square_grid(n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(2, n, x, "  Grid points:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses a square with different sizes in each dimension.
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
        double[] a = { 0.0, -2.0 };
        double[] b = { +10.0, +2.0 };
        int[] c = { 3, 4 };
        int i;
        int n;
        int[] ns = { 3, 3 };
        double[] x;

        n = ns[0] * ns[1];

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Create a grid using SQUARE_GRID.");
        Console.WriteLine("  Use a different physical size in each dimension.");
        Console.WriteLine("  Number of grid points N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("     I    NS     C      A         B");
        Console.WriteLine("");
        for (i = 0; i < 2; i++)
        {
            Console.WriteLine(i + "  "
                                + ns[i].ToString().PadLeft(4) + "  "
                                + c[i].ToString().PadLeft(4) + "  "
                                + a[i].ToString().PadLeft(8) + "  "
                                + b[i].ToString().PadLeft(8) + "");
        }

        x = Grid.square_grid(n, ns, a, b, c);
        typeMethods.r8mat_transpose_print(2, n, x, "  Grid points:");
    }
}