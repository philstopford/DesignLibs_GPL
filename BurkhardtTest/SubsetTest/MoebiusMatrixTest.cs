using System;
using Burkardt;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class MoebiusMatrixTest
{
    public static void moebius_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOEBIUS_MATRIX_TEST tests MOEBIUS_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 11;

        int[] ih =
        {
            0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0
        };
        int[] matrix = new int[N * N];

        Console.WriteLine("");
        Console.WriteLine("MOEBIUS_MATRIX_TEST");
        Console.WriteLine("  MOEBIUS_MATRIX computes the Moebius matrix.");

        typeMethods.i4mat_print(N, N, ih, "  The input matrix:");

        MoebiusMatrix.moebius_matrix(N, ih, ref matrix);

        typeMethods.i4mat_print(N, N, matrix, "  The Moebius matrix:");

    }

}