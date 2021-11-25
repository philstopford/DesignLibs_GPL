using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class MatrixTest
{
    public static void matrix_product_opt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MATRIX_PRODUCT_OPT_TEST tests MATRIX_PRODUCT_OPT.
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
        const int N = 6;

        int cost = 0;
        int i;
        int[] order = new int[N - 1];
        int[] rank = { 4, 2, 3, 1, 2, 2, 3 };

        Console.WriteLine("");
        Console.WriteLine("MATRIX_PRODUCT_OPT_TEST");
        Console.WriteLine("  MATRIX_PRODUCT_OPT seeks the optimal order");
        Console.WriteLine("  for a chain of matrix products.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix ranks:");
        Console.WriteLine("");
        Console.WriteLine("   I    R    C");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(5) + "  "
                                                      + rank[i].ToString().PadLeft(5) + "  "
                                                      + rank[i + 1].ToString().PadLeft(5) + "");
        }

        MatrixProduct.matrix_product_opt(N, rank, ref cost, ref order);

        Console.WriteLine("");
        Console.WriteLine("  Optimal cost is " + cost + "");

        typeMethods.i4vec1_print(N - 1, order, "  Ordering:");
    }
}