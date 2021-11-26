﻿using System;
using Burkardt.Types;

namespace SparseTripletTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ST_IO_TEST.
        //
        //  Discussion:
        //
        //    ST_IO_TEST tests the ST_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ST_IO_TEST:");
        Console.WriteLine("  Test the ST_IO library.");
        //
        //  Real double precision tests.
        //
        r8st_write_test();
        r8st_read_test();
        r8st_sort_a_test();

        Console.WriteLine("");
        Console.WriteLine("ST_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void r8st_write_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_WRITE_TEST tests R8ST_WRITE.
        //
        //  Discussion:
        //
        //    The matrix is:
        //
        //      11  12   0   0  15
        //      21  22   0   0   0
        //       0   0  33   0  35
        //       0   0   0  44   0
        //      51   0  53   0  55
        //
        //    The index vectors are 1 based, and so have to be converted to
        //    0-base before being written.
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
        double[] ast =  {
                51.0, 12.0, 11.0, 33.0, 15.0,
                53.0, 55.0, 22.0, 35.0, 44.0,
                21.0
            }
            ;
        int[] ist =  {
                5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2
            }
            ;
        int[] jst =  {
                1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1
            }
            ;
        const int m = 5;
        const int n = 5;
        const int nst = 11;
        const string output_filename = "a5by5_r8.st";

        Console.WriteLine("");
        Console.WriteLine("R8ST_WRITE_TEST");
        Console.WriteLine("  R8ST_WRITE writes an R8ST file.");

        typeMethods.i4vec_dec(nst, ref ist);
        typeMethods.i4vec_dec(nst, ref jst);

        int i_min = typeMethods.i4vec_min(nst, ist);
        int i_max = typeMethods.i4vec_max(nst, ist);
        int j_min = typeMethods.i4vec_min(nst, jst);
        int j_max = typeMethods.i4vec_max(nst, jst);

        typeMethods.r8st_header_print(i_min, i_max, j_min, j_max, m, n, nst);

        typeMethods.r8st_print(m, n, nst, ist, jst, ast, "  Sparse Triple data:");

        typeMethods.r8st_write(output_filename, m, n, nst, ist, jst, ast);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the matrix data to '" + output_filename + "'");

    }

    private static void r8st_read_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_READ_TEST tests R8ST_HEADER_READ, R8ST_DATA_READ.
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
        int i_max = 0;
        int i_min = 0;
        const string input_filename = "kershaw_r8.st";
        int j_max = 0;
        int j_min = 0;
        int m = 0;
        int n = 0;
        int nst = 0;

        Console.WriteLine("");
        Console.WriteLine("R8ST_READ_TEST");
        Console.WriteLine("  R8ST_HEADER_READ reads the header from an R8ST file.");
        Console.WriteLine("  R8ST_DATA_READ reads the data from an R8ST file.");
        Console.WriteLine("");
        Console.WriteLine("  Read the data from '" + input_filename + "'");

        typeMethods.r8st_header_read(input_filename, ref i_min, ref i_max, ref j_min, ref j_max, ref m, ref n, ref nst);

        typeMethods.r8st_header_print(i_min, i_max, j_min, j_max, m, n, nst);

        double[] ast = new double[nst];
        int[] ist = new int[nst];
        int[] jst = new int[nst];

        typeMethods.r8st_data_read(input_filename, m, n, nst, ref ist, ref jst, ref ast);

        typeMethods.r8st_print(m, n, nst, ist, jst, ast,
            "  Sparse Triplet data read from file:");

    }

    private static void r8st_sort_a_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8ST_SORT_A_TEST tests R8ST_SORT_A.
        //
        //  Discussion:
        //
        //    The matrix is:
        //
        //      11  12   0   0  15
        //      21  22   0   0   0
        //       0   0  33   0  35
        //       0   0   0  44   0
        //      51   0  53   0  55
        //
        //    The index vectors are 1 based, and so have to be converted to
        //    0-base before being written.
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
        double[] ast =  {
                51.0, 12.0, 11.0, 33.0, 15.0,
                53.0, 55.0, 22.0, 35.0, 44.0,
                21.0
            }
            ;
        int[] ist =  {
                5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2
            }
            ;
        int[] jst =  {
                1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1
            }
            ;
        const int m = 5;
        const int n = 5;
        const int nst = 11;

        Console.WriteLine("");
        Console.WriteLine("R8ST_SORT_A_TEST");
        Console.WriteLine("  R8ST_SORT_A sorts an R8ST matrix by columns.");

        int i_min = typeMethods.i4vec_min(nst, ist);
        int i_max = typeMethods.i4vec_max(nst, ist);
        int j_min = typeMethods.i4vec_min(nst, jst);
        int j_max = typeMethods.i4vec_max(nst, jst);

        typeMethods.r8st_header_print(i_min, i_max, j_min, j_max, m, n, nst);

        typeMethods.r8st_print(m, n, nst, ist, jst, ast, "  Matrix data before sorting:");

        typeMethods.r8st_sort_a(m, n, nst, ref ist, ref jst, ref ast);

        typeMethods.r8st_print(m, n, nst, ist, jst, ast, "  Matrix data sorted by column:");
    }
}