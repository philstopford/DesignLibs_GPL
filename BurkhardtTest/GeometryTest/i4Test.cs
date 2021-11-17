using System;
using Burkardt.Types;

namespace GeometryTest;

public static class i4Test
{
    public static void test0322()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0322 tests I4COL_FIND_ITEM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 5;
        int N = 4;
        int TEST_NUM = 3;

        int[] a = new int[M * N];
        int col = 0;
        int i;
        int item;
        int[] item_test = {34, 12, 90};
        int j;
        int row = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0322");
        Console.WriteLine("  I4COL_FIND_ITEM finds the first occurrence of");
        Console.WriteLine("  an item in an integer array of columns.");

        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * M] = 10 * (i + 1) + j + 1;
            }
        }

        typeMethods.i4mat_print(M, N, a, "  The matrix of columns:");

        for (test = 0; test < TEST_NUM; test++)
        {
            item = item_test[test];

            typeMethods.i4col_find_item(M, N, a, item, ref row, ref col);

            Console.WriteLine("  Item " + item
                                        + " occurs in row " + row
                                        + " and column " + col + "");
        }
    }

    public static void test0323()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0323 tests I4COL_FIND_PAIR_WRAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 5;
        int N = 4;
        int TEST_NUM = 5;

        int[] a = new int[M * N];
        int col = 0;
        int i;
        int item1;
        int[] item1_test = {22, 32, 22, 54, 54};
        int item2;
        int[] item2_test = {32, 22, 23, 14, 11};
        int j;
        int row = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0323");
        Console.WriteLine("  I4COL_FIND_PAIR_WRAP finds the first occurrence of");
        Console.WriteLine("  a pair of item in an integer array of columns.");
        Console.WriteLine("  Items in the array are ordered by column, and");
        Console.WriteLine("  wraparound is allowed.");

        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * M] = 10 * (i + 1) + j + 1;
            }
        }

        typeMethods.i4mat_print(M, N, a, "  The matrix of columns:");

        for (test = 0; test < TEST_NUM; test++)
        {
            item1 = item1_test[test];
            item2 = item2_test[test];

            typeMethods.i4col_find_pair_wrap(M, N, a, item1, item2, ref row, ref col);

            Console.WriteLine("  Item " + item1
                                        + " followed by item " + item2
                                        + " occurs in row " + row
                                        + " and column " + col + "");
        }
    }

}