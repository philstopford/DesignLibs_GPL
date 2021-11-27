using System;
using BLASTestData;
using Burkardt.BLAS;
using Burkardt.Types;

namespace BLAS3DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BLAS3_D_TEST.
        //
        //  Discussion:
        //
        //    BLAS3_D_TEST tests the BLAS3_D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 March 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BLAS3_D_TEST");
        Console.WriteLine("  Test the BLAS3_D library.");

        dgemm_test();
        dtrmm_test();
        dtrsm_test();

        Console.WriteLine("");
        Console.WriteLine("BLAS3_D_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void dgemm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEMM_TEST tests DGEMM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("DGEMM_TEST");
        Console.WriteLine("  DGEMM carries out matrix multiplications");
        Console.WriteLine("  for double precision real matrices.");
        Console.WriteLine("");
        Console.WriteLine("  1: C = alpha * A  * B  + beta * C;");
        Console.WriteLine("  2: C = alpha * A' * B  + beta * C;");
        Console.WriteLine("  3: C = alpha * A  * B' + beta * C;");
        Console.WriteLine("  4: C = alpha * A' * B' + beta * C;");
        Console.WriteLine("");
        Console.WriteLine("  We carry out all four calculations, but in each case,");
        Console.WriteLine("  we choose our input matrices so that we get the same result.");
        //
        //  C = alpha * A * B + beta * C.
        //
        char transa = 'N';
        char transb = 'N';
        char transc = 'N';
        int m = 4;
        int n = 5;
        int k = 3;
        double alpha = 2.0;
        int lda = m;
        double[] a = BLASData.r8mat_test(transa, lda, m, k);
        int ldb = k;
        double[] b = BLASData.r8mat_test(transb, ldb, k, n);
        double beta = 3.0;
        int ldc = m;
        double[] c = BLASData.r8mat_test(transc, ldc, m, n);

        BLAS3D.dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, ref c, ldc);

        typeMethods.r8mat_print(m, n, c, "  C = alpha * A * B + beta * C:");

        //
        //  C = alpha * A' * B + beta * C.
        //
        transa = 'T';
        transb = 'N';
        transc = 'N';
        m = 4;
        n = 5;
        k = 3;
        alpha = 2.0;
        lda = k;
        a = BLASData.r8mat_test(transa, lda, m, k);
        ldb = k;
        b = BLASData.r8mat_test(transb, ldb, k, n);
        beta = 3.0;
        ldc = m;
        c = BLASData.r8mat_test(transc, ldc, m, n);

        BLAS3D.dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, ref c, ldc);

        typeMethods.r8mat_print(m, n, c, "  C = alpha * A' * B + beta * C:");

        //
        //  C = alpha * A * B' + beta * C.
        //
        transa = 'N';
        transb = 'T';
        transc = 'N';
        m = 4;
        n = 5;
        k = 3;
        alpha = 2.0;
        lda = m;
        a = BLASData.r8mat_test(transa, lda, m, k);
        ldb = n;
        b = BLASData.r8mat_test(transb, ldb, k, n);
        beta = 3.0;
        ldc = m;
        c = BLASData.r8mat_test(transc, ldc, m, n);

        BLAS3D.dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, ref c, ldc);

        typeMethods.r8mat_print(m, n, c, "  C = alpha * A * B' + beta * C:");

        //
        //  C = alpha * A' * B' + beta * C.
        //
        transa = 'T';
        transb = 'T';
        transc = 'N';
        m = 4;
        n = 5;
        k = 3;
        alpha = 2.0;
        lda = k;
        a = BLASData.r8mat_test(transa, lda, m, k);
        ldb = n;
        b = BLASData.r8mat_test(transb, ldb, k, n);
        beta = 3.0;
        ldc = m;
        c = BLASData.r8mat_test(transc, ldc, m, n);

        BLAS3D.dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, ref c, ldc);

        typeMethods.r8mat_print(m, n, c, "  C = alpha * A' * B' + beta * C:");

    }

    private static void dtrmm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRMM_TEST tests DTRMM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("DTRMM_TEST");
        Console.WriteLine("  DTRMM multiplies a triangular matrix A and a");
        Console.WriteLine("  rectangular matrix B");
        Console.WriteLine("");
        Console.WriteLine("  1: B = alpha * A  * B;");
        Console.WriteLine("  2: B = alpha * A' * B;");
        //
        //  B = alpha * A * B.
        //
        char side = 'L';
        char uplo = 'U';
        char transa = 'N';
        char diag = 'N';
        int m = 4;
        int n = 5;
        double alpha = 2.0;
        int lda = m;
        int ldb = m;

        double[] a = new double[lda * m];
        for (j = 0; j < m; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < m; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        char transb = 'N';
        double[] b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  B = alpha * A * B:");

        //
        //  B = alpha * A' * B.
        //
        side = 'L';
        uplo = 'U';
        transa = 'T';
        diag = 'N';
        m = 4;
        n = 5;
        alpha = 2.0;
        lda = m;
        ldb = m;

        a = new double[lda * m];
        for (j = 0; j < m; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < m; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        transb = 'N';
        b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  B = alpha * A' * B:");

    }

    private static void dtrsm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRSM_TEST tests DTRSM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("DTRSM_TEST");
        Console.WriteLine("  DTRSM solves a linear system involving a triangular");
        Console.WriteLine("  matrix A and a rectangular matrix B.");
        Console.WriteLine("");
        Console.WriteLine("  1: Solve A  * X  = alpha * B;");
        Console.WriteLine("  2: Solve A' * X  = alpha * B;");
        Console.WriteLine("  3: Solve X  * A  = alpha * B;");
        Console.WriteLine("  4: Solve X  * A' = alpha * B;");
        //
        //  Solve A * X = alpha * B.
        //
        char side = 'L';
        char uplo = 'U';
        char transa = 'N';
        char diag = 'N';
        int m = 4;
        int n = 5;
        double alpha = 2.0;
        int lda = m;
        int ldb = m;

        double[] a = new double[lda * m];

        for (j = 0; j < m; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < m; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        char transb = 'N';
        double[] b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  X = inv ( A ) * alpha * B:");

        //
        //  Solve A' * X = alpha * B.
        //
        side = 'L';
        uplo = 'U';
        transa = 'T';
        diag = 'N';
        m = 4;
        n = 5;
        alpha = 2.0;
        lda = m;
        ldb = m;

        a = new double[lda * m];
        for (j = 0; j < m; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < m; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        transb = 'N';
        b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  X = inv ( A' ) * alpha * B:");

        //
        //  Solve X * A = alpha * B.
        //
        side = 'R';
        uplo = 'U';
        transa = 'N';
        diag = 'N';
        m = 4;
        n = 5;
        alpha = 2.0;
        lda = n;
        ldb = m;

        a = new double[lda * n];
        for (j = 0; j < n; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < n; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        transb = 'N';
        b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  X = alpha * B * inv ( A ):");

        //
        //  Solve X * A'' = alpha * B.
        //
        side = 'R';
        uplo = 'U';
        transa = 'T';
        diag = 'N';
        m = 4;
        n = 5;
        alpha = 2.0;
        lda = n;
        ldb = m;

        a = new double[lda * n];
        for (j = 0; j < n; j++)
        {
            for (i = 0; i <= j; i++)
            {
                a[i + j * lda] = i + j + 2;
            }

            for (i = j + 1; i < n; i++)
            {
                a[i + j * lda] = 0.0;
            }
        }

        transb = 'N';
        b = BLASData.r8mat_test(transb, ldb, m, n);

        BLAS3D.dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, ref b, ldb);

        typeMethods.r8mat_print(m, n, b, "  X = alpha * B * inv ( A' ):");
    }
}