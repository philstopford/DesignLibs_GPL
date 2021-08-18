using System;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace SimplexCoordsTest
{
    class Program
    {
        static void Main(string[] args)
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
            int n;

            Console.WriteLine("");
            Console.WriteLine("SIMPLEX_COORDINATES_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the SIMPLEX_COORDINATES library.");

            n = 3;
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

        static void simplex_coordinates1_test(int n)

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
            int k;
            double side;
            double volume;
            double volume2;
            double[] x;
            double[] xtx;

            Console.WriteLine("");
            Console.WriteLine("SIMPLEX_COORDINATES1_TEST");
            Console.WriteLine("  Call SIMPLEX_COORDINATES1");

            x = Coordinates.simplex_coordinates1(n);

            typeMethods.r8mat_transpose_print(n, n + 1, x, "  Simplex vertex coordinates:");

            side = 0.0;
            for (i = 0; i < n; i++)
            {
                side = side + Math.Pow(x[i + 0 * n] - x[i + 1 * n], 2);
            }

            side = Math.Sqrt(side);

            volume = Coordinates.simplex_volume(n, x);

            volume2 = Math.Sqrt((double)(n + 1)) / typeMethods.r8_factorial(n)
                                                 / Math.Sqrt(Math.Pow(2.0, n)) * Math.Pow(side, n);

            Console.WriteLine("");
            Console.WriteLine("  Side length =     " + side + "");
            Console.WriteLine("  Volume =          " + volume + "");
            Console.WriteLine("  Expected volume = " + volume2 + "");

            xtx = new double[(n + 1) * (n + 1)];

            for (j = 0; j < n + 1; j++)
            {
                for (i = 0; i < n + 1; i++)
                {
                    xtx[i + j * (n + 1)] = 0.0;
                    for (k = 0; k < n; k++)
                    {
                        xtx[i + j * (n + 1)] = xtx[i + j * (n + 1)] + x[k + i * n] * x[k + j * n];
                    }
                }
            }

            typeMethods.r8mat_transpose_print(n + 1, n + 1, xtx, "  Dot product matrix:");

        }

        static void simplex_coordinates2_test(int n)

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
            int k;
            double side;
            double volume;
            double volume2;
            double[] x;
            double[] xtx;

            Console.WriteLine("");
            Console.WriteLine("SIMPLEX_COORDINATES2_TEST");
            Console.WriteLine("  Call SIMPLEX_COORDINATES2");

            x = Coordinates.simplex_coordinates2(n);

            typeMethods.r8mat_transpose_print(n, n + 1, x, "  Simplex vertex coordinates:");

            side = 0.0;
            for (i = 0; i < n; i++)
            {
                side = side + Math.Pow(x[i + 0 * n] - x[i + 1 * n], 2);
            }

            side = Math.Sqrt(side);

            volume = Coordinates.simplex_volume(n, x);

            volume2 = Math.Sqrt((double)(n + 1)) / typeMethods.r8_factorial(n)
                                                 / Math.Sqrt(Math.Pow(2.0, n)) * Math.Pow(side, n);

            Console.WriteLine("");
            Console.WriteLine("  Side length =     " + side + "");
            Console.WriteLine("  Volume =          " + volume + "");
            Console.WriteLine("  Expected volume = " + volume2 + "");

            xtx = new double[(n + 1) * (n + 1)];

            for (j = 0; j < n + 1; j++)
            {
                for (i = 0; i < n + 1; i++)
                {
                    xtx[i + j * (n + 1)] = 0.0;
                    for (k = 0; k < n; k++)
                    {
                        xtx[i + j * (n + 1)] = xtx[i + j * (n + 1)] + x[k + i * n] * x[k + j * n];
                    }
                }
            }

            typeMethods.r8mat_transpose_print(n + 1, n + 1, xtx, "  Dot product matrix:");
        }
    }
}