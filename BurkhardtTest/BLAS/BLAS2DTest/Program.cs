using System;
using BLASTestData;
using Burkardt.BLAS;
using Burkardt.Types;

namespace BLAS2DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BLAS2_D_TEST.
        //
        //  Discussion:
        //
        //    BLAS2_D_TEST tests the BLAS2_D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BLAS2_D_TEST");
        Console.WriteLine("  Test the BLAS2_D library.");

        dgemv_test();
        dger_test();
        dtrmv_test();

        Console.WriteLine("");
        Console.WriteLine("BLAS2_D_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void dgemv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEMV_TEST tests DGEMV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("DGEMV_TEST");
        Console.WriteLine("  DGEMV computes y := alpha * A * x + beta * y");
        Console.WriteLine("  or             y := alpha * A'' * x + beta * y,");
        Console.WriteLine("  for a general matrix A,");
        //
        //  y = alpha * A * x + beta * y
        //
        char trans = 'N';
        int m = 5;
        int n = 4;
        double alpha = 2.0;
        int lda = m;
        double[] a = BLASData.r8mat_test(trans, lda, m, n);
        double[] x = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = i + 1;
        }

        int incx = 1;
        double beta = 3.0;
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
        {
            y[i] = 10 * (i + 1);
        }

        int incy = 1;

        typeMethods.r8mat_print(m, n, a, "  Matrix A:");
        typeMethods.r8vec_print(n, x, "  Vector X:");
        typeMethods.r8vec_print(m, y, "  Vector Y:");

        BLAS2D.dgemv(trans, m, n, alpha, a, lda, x, incx, beta, ref y, incy);

        typeMethods.r8vec_print(m, y, "  Result Y = alpha * A  * x + beta * y");

        //
        //  y = alpha * A' * x + beta * y
        //
        trans = 'T';
        m = 5;
        n = 4;
        alpha = 2.0;
        lda = m;
        a = BLASData.r8mat_test(trans, lda, n, m);
        x = new double[m];
        for (int i = 0; i < m; i++)
        {
            x[i] = i + 1;
        }

        incx = 1;
        beta = 3.0;
        y = new double[n];
        for (int i = 0; i < n; i++)
        {
            y[i] = 10 * (i + 1);
        }

        incy = 1;

        typeMethods.r8mat_print(m, n, a, "  Matrix A:");
        typeMethods.r8vec_print(m, x, "  Vector X:");
        typeMethods.r8vec_print(n, y, "  Vector Y:");

        BLAS2D.dgemv(trans, m, n, alpha, a, lda, x, incx, beta, ref y, incy);

        typeMethods.r8vec_print(n, y, "  Result Y = alpha * A  * x + beta * y");
    }

    private static void dger_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGER_TEST tests DGER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("DGER_TEST");
        Console.WriteLine("  DGER computes A := A + alpha * x * y'");
        Console.WriteLine("  for a general matrix A.");

        int m = 5;
        int n = 4;
        double alpha = 2.0;
        char trans = 'N';
        int lda = m;
        double[] a = BLASData.r8mat_test(trans, lda, m, n);

        double[] x = new double[m];
        for (int i = 0; i < m; i++)
        {
            x[i] = i + 1;
        }

        int incx = 1;

        double[] y = new double[n];
        for (int i = 0; i < n; i++)
        {
            y[i] = 10 * (i + 1);
        }

        int incy = 1;

        typeMethods.r8mat_print(m, n, a, "  Matrix A:");
        typeMethods.r8vec_print(m, x, "  Vector X:");
        typeMethods.r8vec_print(n, y, "  Vector Y:");

        BLAS2D.dger(m, n, alpha, x, incx, y, incy, ref a, lda);

        typeMethods.r8mat_print(m, n, a, "  Result A = A + alpha * x * y");
    }

    private static void dtrmv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRMV_TEST tests DTRMV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int lda = 5;
        int m = 5;
        int n = 5;

        double[] a = new double[lda * n];
        double[] x = new double[n];

        Console.WriteLine("");
        Console.WriteLine("DTRMV_TEST");
        Console.WriteLine("  DTRMV computes y := A * x or y := A' * x");
        Console.WriteLine("  For a triangular matrix A.");

        for (int test = 1; test <= 2; test++)
        {
            char uplo = 'U';

            char trans = test switch
            {
                1 => 'N',
                _ => 'T'
            };

            char diag = 'N';

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i <= j; i++)
                {
                    a[i + j * m] = i + j + 2;
                }

                for (int i = j + 1; i < m; i++)
                {
                    a[i + j * m] = 0.0;
                }
            }

            int incx = 1;
            for (int i = 0; i < n; i++)
            {
                x[i] = i + 1;
            }

            BLAS2D.dtrmv(uplo, trans, diag, n, a, lda, ref x, incx);

            switch (trans)
            {
                case 'N':
                    typeMethods.r8vec_print(n, x, "  Result y = A * x");
                    break;
                default:
                    typeMethods.r8vec_print(n, x, "  Result y = A' * x");
                    break;
            }
        }

    }
}