using System;
using Burkardt.Types;

namespace SubsetTestNS;

public static class DecMatTest
{

    public static void decmat_det_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DECMAT_DET_TEST tests DECMAT_DET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N3 = 3;
        const int N4 = 4;

        int a = 0;
        int b = 0;
        int[] a3 = new int[N3 * N3];
        int[] a4 = new int[N4 * N4];
        int[] b3 = new int[N3 * N3];
        int[] b4 = new int[N4 * N4];
        int i;
        int dbot = 0;
        int dtop = 0;
        int j;

        Console.WriteLine("");
        Console.WriteLine("DECMAT_DET_TEST");
        Console.WriteLine("  DECMAT_DET: determinant of a decimal matrix.");
        Console.WriteLine("");

        const int dec_digit = 5;

        int k = 0;
        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                k += 1;
                a3[i + j * N3] = k;
            }
        }

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                b3[i + j * N3] = 0;
            }
        }

        typeMethods.decmat_print(N3, N3, a3, b3, "  The 123/456/789 matrix:");

        typeMethods.decmat_det(N3, a3, b3, dec_digit, ref dtop, ref dbot);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the 123/456/789 matrix = "
                          + dtop + "* 10^("
                          + dbot + ")");

        for (i = 0; i < N4; i++)
        {
            for (j = 0; j < N4; j++)
            {
                double r = 1.0 / (i + j + 2);

                typeMethods.r8_to_dec(r, dec_digit, ref a, ref b);
                a4[i + j * N4] = a;
                b4[i + j * N4] = b;
            }
        }

        typeMethods.decmat_print(N4, N4, a4, b4, "  The Hilbert matrix:");

        typeMethods.decmat_det(N4, a4, b4, dec_digit, ref dtop, ref dbot);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the Hilbert matrix = "
                          + dtop + "* 10^("
                          + dbot + ")");

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                if (i == j)
                {
                    a3[i + j * N3] = 2;
                }
                else if (i == j + 1 || i == j - 1)
                {
                    a3[i + j * N3] = -1;
                }
                else
                {
                    a3[i + j * N3] = 0;
                }
            }
        }

        for (i = 0; i < N3; i++)
        {
            for (j = 0; j < N3; j++)
            {
                b3[i + j * N3] = 0;
            }
        }

        typeMethods.decmat_print(N3, N3, a3, b3, "  The -1,2,-1 matrix:");

        typeMethods.decmat_det(N3, a3, b3, dec_digit, ref dtop, ref dbot);

        Console.WriteLine("");
        Console.WriteLine("  Determinant of the -1,2,-1 matrix = "
                          + dtop + "* 10^("
                          + dbot + ")");
    }

    public static void decmat_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DECMAT_PRINT_TEST tests DECMAT_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int a = 0;
        int[] amat = new int[4 * 3];
        int b = 0;
        int[] bmat = new int[4 * 3];
        int i;
        const int m = 4;
        const int n = 3;

        Console.WriteLine("");
        Console.WriteLine("DECMAT_PRINT_TEST");
        Console.WriteLine("  DECMAT_PRINT prints a decimal matrix.");

        const int dec_digit = 5;

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                double r = 1.0 / (i + j + 2);

                typeMethods.r8_to_dec(r, dec_digit, ref a, ref b);
                amat[i + j * m] = a;
                bmat[i + j * n] = b;
            }
        }

        typeMethods.decmat_print(m, n, amat, bmat, "  The Hilbert matrix:");
    }

}