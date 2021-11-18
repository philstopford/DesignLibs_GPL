using System;
using Burkardt.Storage;

namespace CompressedColumnToSparseTripletTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ccs_TO_ST_TEST.
        //
        //  Discussion:
        //
        //    ccs_TO_ST_TEST tests the ccs_TO_ST library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ccs_TO_ST_TEST");
        Console.WriteLine("  Test the ccs_TO_ST library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("ccs_TO_ST_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests ccs_TO_ST using a 1-based matrix.
        //
        //  Discussion:
        //
        //    This test uses a trivial matrix whose full representation is:
        //
        //          2  3  0  0  0
        //          3  0  4  0  6
        //      A = 0 -1 -3  2  0
        //          0  0  1  0  0
        //          0  4  2  0  1
        //
        //    The 1-based CCS representation is
        //
        //      #  ICC  CCC  ACC
        //     --  ---  ---  ---
        //      1    1    1    2
        //      2    2         3
        //
        //      3    1    3    3
        //      4    3        -1
        //      5    5         4
        //
        //      6    2    6    4
        //      7    3        -3
        //      8    4         1
        //      9    5         2
        //
        //     10    3   10    2
        //
        //     11    2   11    6
        //     12    5         1
        //
        //     13    *   13
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int NCC = 12;

        double[] acc =  {
                2.0, 3.0,
                3.0, -1.0, 4.0,
                4.0, -3.0, 1.0, 2.0,
                2.0,
                6.0, 1.0
            }
            ;
        double[] ast;
        int[] ccc =  {
                1, 3, 6, 10, 11, 13
            }
            ;
        int[] icc =  {
                1, 2,
                1, 3, 5,
                2, 3, 4, 5,
                3,
                2, 5
            }
            ;
        int[] ist;
        int[] jst;
        int m = 5;
        int n = N;
        int ncc = NCC;
        int nst = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Convert a 1-based CCS matrix to ST format.");
        //
        //  Print the CCS matrix.
        //
        CompressedColumnStorage.ccs_print(m, n, ncc, icc, ccc, acc, "  The CCS matrix:");
        //
        //  Convert it.
        //
        ist = new int[ncc];
        jst = new int[ncc];
        ast = new double[ncc];

        CompressedColumnStorage.ccs_to_st(m, n, ncc, icc, ccc, acc, ref nst, ref ist, ref jst, ref ast);
        //
        //  Print the ST matrix.
        //
        SparseTriplet.st_print(m, n, nst, ist, jst, ast, "  The ST matrix:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests ccs_TO_ST using a 0-based matrix.
        //
        //  Discussion:
        //
        //    This test uses a trivial matrix whose full representation is:
        //
        //          2  3  0  0  0
        //          3  0  4  0  6
        //      A = 0 -1 -3  2  0
        //          0  0  1  0  0
        //          0  4  2  0  1
        //
        //    The 0-based CCS representation is
        //
        //      #  ICC  CCC  ACC
        //     --  ---  ---  ---
        //      0    0    0    2
        //      1    1         3
        //
        //      2    0    2    3
        //      3    2        -1
        //      4    4         4
        //
        //      5    1    5    4
        //      6    2        -3
        //      7    3         1
        //      8    4         2
        //
        //      9    2    9    2
        //
        //     10    1   10    6
        //     11    4         1
        //
        //     12    *   12
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;
        int NCC = 12;

        double[] acc =  {
                2.0, 3.0,
                3.0, -1.0, 4.0,
                4.0, -3.0, 1.0, 2.0,
                2.0,
                6.0, 1.0
            }
            ;
        double[] ast;
        int[] ccc =  {
                0, 2, 5, 9, 10, 12
            }
            ;
        int[] icc =  {
                0, 1,
                0, 2, 4,
                1, 2, 3, 4,
                2,
                1, 4
            }
            ;
        int[] ist;
        int[] jst;
        int m = 5;
        int n = N;
        int ncc = NCC;
        int nst = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Convert a 0-based CCS matrix to ST format.");
        //
        //  Print the CCS matrix.
        //
        CompressedColumnStorage.ccs_print(m, n, ncc, icc, ccc, acc, "  The CCS matrix:");
        //
        //  Convert it.
        //
        ist = new int[ncc];
        jst = new int[ncc];
        ast = new double[ncc];

        CompressedColumnStorage.ccs_to_st(m, n, ncc, icc, ccc, acc, ref nst, ref ist, ref jst, ref ast);
        //
        //  Print the ST matrix.
        //
        SparseTriplet.st_print(m, n, nst, ist, jst, ast, "  The ST matrix:");
    }
}