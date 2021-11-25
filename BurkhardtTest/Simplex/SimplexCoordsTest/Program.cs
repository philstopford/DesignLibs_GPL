using System;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace SimplexCoordsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SIMPLEX_COORDINATED_TEST.
        //
        //  Discussion:
        //
        //    SIMPLEX_COORDINATES_TEST tests the SIMPLEX_COORDINATES library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_COORDINATES_TEST");
            
        Console.WriteLine("  Test the SIMPLEX_COORDINATES library.");

        int n = 3;
        simplex_coordinates1_test(n);
        simplex_coordinates2_test(n);

        n = 4;
        simplex_coordinates1_test(n);
        simplex_coordinates2_test(n);
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_COORDINATES_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void simplex_coordinates1_test(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_COORDINATES1_TEST tests SIMPLEX_COORDINATES1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_COORDINATES1_TEST");
        Console.WriteLine("  Call SIMPLEX_COORDINATES1");

        double[] x = Coordinates.simplex_coordinates1(n);

        typeMethods.r8mat_transpose_print(n, n + 1, x, "  Simplex vertex coordinates:");

        double side = 0.0;
        for (i = 0; i < n; i++)
        {
            side += Math.Pow(x[i + 0 * n] - x[i + 1 * n], 2);
        }

        side = Math.Sqrt(side);

        double volume = Coordinates.simplex_volume(n, x);

        double volume2 = Math.Sqrt(n + 1) / typeMethods.r8_factorial(n)
                                          / Math.Sqrt(Math.Pow(2.0, n)) * Math.Pow(side, n);

        Console.WriteLine("");
        Console.WriteLine("  Side length =     " + side + "");
        Console.WriteLine("  Volume =          " + volume + "");
        Console.WriteLine("  Expected volume = " + volume2 + "");

        double[] xtx = new double[(n + 1) * (n + 1)];

        for (j = 0; j < n + 1; j++)
        {
            for (i = 0; i < n + 1; i++)
            {
                xtx[i + j * (n + 1)] = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    xtx[i + j * (n + 1)] += x[k + i * n] * x[k + j * n];
                }
            }
        }

        typeMethods.r8mat_transpose_print(n + 1, n + 1, xtx, "  Dot product matrix:");

    }

    private static void simplex_coordinates2_test(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_COORDINATES2_TEST tests SIMPLEX_COORDINATES2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_COORDINATES2_TEST");
        Console.WriteLine("  Call SIMPLEX_COORDINATES2");

        double[] x = Coordinates.simplex_coordinates2(n);

        typeMethods.r8mat_transpose_print(n, n + 1, x, "  Simplex vertex coordinates:");

        double side = 0.0;
        for (i = 0; i < n; i++)
        {
            side += Math.Pow(x[i + 0 * n] - x[i + 1 * n], 2);
        }

        side = Math.Sqrt(side);

        double volume = Coordinates.simplex_volume(n, x);

        double volume2 = Math.Sqrt(n + 1) / typeMethods.r8_factorial(n)
                                          / Math.Sqrt(Math.Pow(2.0, n)) * Math.Pow(side, n);

        Console.WriteLine("");
        Console.WriteLine("  Side length =     " + side + "");
        Console.WriteLine("  Volume =          " + volume + "");
        Console.WriteLine("  Expected volume = " + volume2 + "");

        double[] xtx = new double[(n + 1) * (n + 1)];

        for (j = 0; j < n + 1; j++)
        {
            for (i = 0; i < n + 1; i++)
            {
                xtx[i + j * (n + 1)] = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    xtx[i + j * (n + 1)] += x[k + i * n] * x[k + j * n];
                }
            }
        }

        typeMethods.r8mat_transpose_print(n + 1, n + 1, xtx, "  Dot product matrix:");
    }
}