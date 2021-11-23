using System;

namespace IndexTest;

using Index = Burkardt.IndexNS.Index;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for INDEX_TEST.
        //
        //  Discussion:
        //
        //    INDEX_TEST tests the INDEX library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("INDEX_TEST:");
        Console.WriteLine("  Test the INDEX library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("INDEX_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests INDEX0 and INDEX1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i_max;
        int i_min;
        int index_max;
        int index_min;
        int value;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  INDEX0 indexes a 1D array with zero base,");
        Console.WriteLine("  INDEX1 indexes a 1D array with  unit base.");
        Console.WriteLine("");
        Console.WriteLine("             Min Index   Max");
        Console.WriteLine("");

        i_min = 1;
        i = 3;
        i_max = 5;
        Console.WriteLine("  1D Index"
                          + "  " + i_min.ToString().PadLeft(4)
                          + "  " + i.ToString().PadLeft(4)
                          + "  " + i_max.ToString().PadLeft(4) + "");

        value = Index.index0(i_min, i, i_max);
        index_min = 0;
        index_max = index_min + i_max - i_min;
        Console.WriteLine("  Index0  "
                          + "  " + index_min.ToString().PadLeft(4)
                          + "  " + value.ToString().PadLeft(4)
                          + "  " + index_max.ToString().PadLeft(4) + "");

        value = Index.index1(i_min, i, i_max);
        index_min = 1;
        index_max = index_min + i_max - i_min;
        Console.WriteLine("  Index1  "
                          + "  " + index_min.ToString().PadLeft(4)
                          + "  " + value.ToString().PadLeft(4)
                          + "  " + index_max.ToString().PadLeft(4) + "");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests INDEX01, INDEX10, INDEX12 and INDEX21.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i_max;
        int i_min;
        int index_max;
        int index_min;
        int j;
        int j_max;
        int j_min;
        int value;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  For a 2D array,");
        Console.WriteLine("  INDEX01 column indexes with zero base,");
        Console.WriteLine("  INDEX10 row indexes with zero base,");
        Console.WriteLine("  INDEX12 column indexes with unit base,");
        Console.WriteLine("  INDEX21 row indexes with unit base.");
        Console.WriteLine("");
        Console.WriteLine("                Min   Index     Max");
        Console.WriteLine("");

        i_min = 1;
        i = 3;
        i_max = 5;
        j_min = 1;
        j = 2;
        j_max = 4;
        Console.WriteLine("  2D Index:"
                          + "  " + i_min.ToString().PadLeft(3) + j_min.ToString().PadLeft(3)
                          + "  " + i.ToString().PadLeft(3) + j.ToString().PadLeft(3)
                          + "  " + i_max.ToString().PadLeft(3) + j_max.ToString().PadLeft(3) + "");

        value = Index.index01(i_min, i, i_max, j_min, j, j_max);
        index_min = 0;
        index_max = index_min + (i_max - i_min + 1) * (j_max - j_min + 1) - 1;
        Console.WriteLine("  INDEX01: "
                          + "  " + index_min.ToString().PadLeft(6)
                          + "  " + value.ToString().PadLeft(6)
                          + "  " + index_max.ToString().PadLeft(6) + "");

        value = Index.index10(i_min, i, i_max, j_min, j, j_max);
        index_min = 0;
        index_max = index_min + (i_max - i_min + 1) * (j_max - j_min + 1) - 1;
        Console.WriteLine("  INDEX10: "
                          + "  " + index_min.ToString().PadLeft(6)
                          + "  " + value.ToString().PadLeft(6)
                          + "  " + index_max.ToString().PadLeft(6) + "");

        value = Index.index12(i_min, i, i_max, j_min, j, j_max);
        index_min = 1;
        index_max = index_min + (i_max - i_min + 1) * (j_max - j_min + 1) - 1;
        Console.WriteLine("  INDEX12: "
                          + "  " + index_min.ToString().PadLeft(6)
                          + "  " + value.ToString().PadLeft(6)
                          + "  " + index_max.ToString().PadLeft(6) + "");

        value = Index.index21(i_min, i, i_max, j_min, j, j_max);
        index_min = 1;
        index_max = index_min + (i_max - i_min + 1) * (j_max - j_min + 1) - 1;
        Console.WriteLine("  INDEX21: "
                          + "  " + index_min.ToString().PadLeft(6)
                          + "  " + value.ToString().PadLeft(6)
                          + "  " + index_max.ToString().PadLeft(6) + "");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests INDEX012, INDEX123, INDEX210, and INDEX321.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i_max;
        int i_min;
        int index_max;
        int index_min;
        int j;
        int j_max;
        int j_min;
        int k;
        int k_max;
        int k_min;
        int m;
        int value;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  For a 3D array,");
        Console.WriteLine("  INDEX012 column indexes with zero base,");
        Console.WriteLine("  INDEX123 column indexes with unit base,");
        Console.WriteLine("  INDEX210 row indexes with zero base.");
        Console.WriteLine("  INDEX321 row indexes with unit base.");
        Console.WriteLine("");
        Console.WriteLine("                   Min      Index        Max");
        Console.WriteLine("");

        i_min = 1;
        i = 3;
        i_max = 5;
        j_min = 1;
        j = 2;
        j_max = 4;
        k_min = 1;
        k = 1;
        k_max = 3;

        m = (i_max - i_min + 1)
            * (j_max - j_min + 1)
            * (k_max - k_min + 1);

        Console.WriteLine("  3D Index:"
                          + "  " + i_min.ToString().PadLeft(3) + j_min.ToString().PadLeft(3) +
                          k_min.ToString().PadLeft(3)
                          + "  " + i.ToString().PadLeft(3) + j.ToString().PadLeft(3) + k.ToString().PadLeft(3)
                          + "  " + i_max.ToString().PadLeft(3) + j_max.ToString().PadLeft(3) +
                          k_max.ToString().PadLeft(3) + "");

        value = Index.index012(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX012:"
                          + "  " + index_min.ToString().PadLeft(9)
                          + "  " + value.ToString().PadLeft(9)
                          + "  " + index_max.ToString().PadLeft(9) + "");

        value = Index.index123(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX123:"
                          + "  " + index_min.ToString().PadLeft(9)
                          + "  " + value.ToString().PadLeft(9)
                          + "  " + index_max.ToString().PadLeft(9) + "");

        value = Index.index210(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX210:"
                          + "  " + index_min.ToString().PadLeft(9)
                          + "  " + value.ToString().PadLeft(9)
                          + "  " + index_max.ToString().PadLeft(9) + "");

        value = Index.index321(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX321:"
                          + "  " + index_min.ToString().PadLeft(9)
                          + "  " + value.ToString().PadLeft(9)
                          + "  " + index_max.ToString().PadLeft(9) + "");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests INDEX0123, INDEX1234, INDEX3210, and INDEX4321.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i_max;
        int i_min;
        int index_max;
        int index_min;
        int j;
        int j_max;
        int j_min;
        int k;
        int k_max;
        int k_min;
        int l;
        int l_max;
        int l_min;
        int m;
        int value;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  For a 4D array,");
        Console.WriteLine("  INDEX0123 column indexes with zero base,");
        Console.WriteLine("  INDEX1234 column indexes with unit base,");
        Console.WriteLine("  INDEX3210 row indexes with zero base,");
        Console.WriteLine("  INDEX4321 row indexes with unit base.");
        Console.WriteLine("");
        Console.WriteLine("                       Min         Index           Max");
        Console.WriteLine("");

        i_min = 1;
        i = 3;
        i_max = 5;
        j_min = 1;
        j = 2;
        j_max = 4;
        k_min = 1;
        k = 1;
        k_max = 3;
        l_min = 1;
        l = 2;
        l_max = 2;

        m = (i_max - i_min + 1)
            * (j_max - j_min + 1)
            * (k_max - k_min + 1)
            * (l_max - l_min + 1);

        Console.WriteLine("  4D Index:  "
                          + "  " + i_min.ToString().PadLeft(3) + j_min.ToString().PadLeft(3) +
                          k_min.ToString().PadLeft(3) + l_min.ToString().PadLeft(3)
                          + "  " + i.ToString().PadLeft(3) + j.ToString().PadLeft(3) + k.ToString().PadLeft(3) +
                          l.ToString().PadLeft(3)
                          + "  " + i_max.ToString().PadLeft(3) + j_max.ToString().PadLeft(3) +
                          k_max.ToString().PadLeft(3) + l_max.ToString().PadLeft(3) + "");

        value = Index.index0123(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX0123: "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.index1234(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX1234: "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.index3210(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX3210: "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.index4321(i_min, i, i_max, j_min, j, j_max, k_min, k, k_max, l_min, l, l_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX4321: "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests INDEX0N, INDEX1N, INDEXN0 and INDEXN1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //
    {
        int[] i = {3, 2, 1, 2};
        int[] i_max = {5, 4, 3, 2};
        int[] i_min = {1, 1, 1, 1};
        int index_max;
        int index_min;
        int j;
        int m;
        int n = 4;
        int value;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  For an N-dimensional array,");
        Console.WriteLine("  INDEX0N column indexes with zero base,");
        Console.WriteLine("  INDEX1N column indexes with unit base,");
        Console.WriteLine("  INDEXN0 row indexes with zero base,");
        Console.WriteLine("  INDEXN1 row indexes with unit base.");
        Console.WriteLine("");
        Console.WriteLine("                       Min         Index           Max");

        m = 1;
        for (j = 0; j < n; j++)
        {
            m *= i_max[j] - i_min[j] + 1;
        }

        Console.WriteLine("  ND Index: "
                          + "  " + i_min[0].ToString().PadLeft(3) + i_min[1].ToString().PadLeft(3) +
                          i_min[2].ToString().PadLeft(3) + i_min[3].ToString().PadLeft(3)
                          + "  " + i[0].ToString().PadLeft(3) + i[1].ToString().PadLeft(3) +
                          i[2].ToString().PadLeft(3) + i[3].ToString().PadLeft(3)
                          + "  " + i_max[0].ToString().PadLeft(3) + i_max[1].ToString().PadLeft(3) +
                          i_max[2].ToString().PadLeft(3) + i_max[3].ToString().PadLeft(3) + "");

        value = Index.index0n(n, i_min, i, i_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX0N:  "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.index1n(n, i_min, i, i_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEX1N:  "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.indexn0(n, i_min, i, i_max);
        index_min = 0;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEXN0:  "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");

        value = Index.indexn1(n, i_min, i, i_max);
        index_min = 1;
        index_max = index_min + m - 1;
        Console.WriteLine("  INDEXN1:  "
                          + "  " + index_min.ToString().PadLeft(12)
                          + "  " + value.ToString().PadLeft(12)
                          + "  " + index_max.ToString().PadLeft(12) + "");
    }
}