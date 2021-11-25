﻿using System;
using Burkardt.Table;
using Burkardt.Types;
using Burkardt.Uniform;

namespace TableTest;

internal static class Program
{
    private static void Main()
    {
        //
        //  Purpose:
        //
        //    MAIN is the main program for TABLE_IO_TEST.
        //
        //  Discussion:
        //
        //    TABLE tests the TABLE_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine();
            Console.WriteLine("TABLE_IO_TEST");
            Console.WriteLine("  Test the TABLE_IO library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();

            Console.WriteLine();
            Console.WriteLine("TABLE_IO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }
    }


    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests R8MAT_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 5;
        const int N = 20;

        int i;
        const string output_filename = "r8mat_05_00020.txt";
        double[] table = new double[M * N];

        Console.WriteLine();
        Console.WriteLine("TEST01");
        Console.WriteLine("  R8MAT_WRITE writes an R8MAT file.");

        for (i = 0; i < M; i++)
        {
            int j;
            for (j = 0; j < N; j++)
            {
                table[i + j * M] = (100 * (j + 1) + i + 1) / 10.0;
            }
        }

        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + M);
        Console.WriteLine("  Number of points N  = " + N);

        typeMethods.r8mat_print_some(M, N, table, 1, 1, 5, 5,
            "  5x5 portion of the data written to file:");

        typeMethods.r8mat_transpose_print_some(M, N, table, 1, 1, 5, 5,
            "  5x5 portion of the TRANSPOSED data:");

        typeMethods.r8mat_write(output_filename, M, N, table);

        Console.WriteLine();
        Console.WriteLine("  Wrote the header and data for \""
                          + output_filename + "\"");
    }


    private static void test02 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests R8MAT_HEADER_READ and R8MAT_DATA_READ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string input_filename = "r8mat_05_00020.txt";

        Console.WriteLine();
        Console.WriteLine("TEST02");
        Console.WriteLine("  For data stored in an R8MAT file,");
        Console.WriteLine("  R8MAT_HEADER_READ reads the header");
        Console.WriteLine("  (Information about the dimension of the data)");
        Console.WriteLine("  R8MAT_DATA_READ reads the data.");

        TableHeader h = typeMethods.r8mat_header_read ( input_filename );

        Console.WriteLine();
        Console.WriteLine("  Read the header of \"" + input_filename + "\".");
        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + h.m);
        Console.WriteLine("  Number of points N  = " + h.n);

        double[] table = typeMethods.r8mat_data_read ( input_filename, h.m, h.n );

        Console.WriteLine();
        Console.WriteLine("  Read the data in \"" + input_filename + "\".");

        typeMethods.r8mat_print_some ( h.m, h.n, table, 1, 1, 5, 5, 
            "  5x5 portion of data read from file:" );
    }


    private static void test03 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests I4MAT_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 5;
        const int N = 20;

        const string output_filename = "i4mat_05_00020.txt";
        int[] table = new int[M*N];

        Console.WriteLine();
        Console.WriteLine("TEST03");
        Console.WriteLine("  I4MAT_WRITE writes an I4MAT file.");

        for (int i = 0; i < M; i++ )
        {
            for (int j = 0; j < N; j++ )
            {
                table[i+j*M] = 100 * ( j + 1 ) + i + 1;
            }
        }

        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + M);
        Console.WriteLine("  Number of points N  = " + N);

        typeMethods.i4mat_print_some ( M, N, table, 1, 1, 5, 5, 
            "  5 x 5 portion of data written to file:" );

        typeMethods.i4mat_write ( output_filename, M, N, table );

        Console.WriteLine();
        Console.WriteLine("  Wrote the header and data for \""
                          + output_filename + "\"");
    }


    private static void test04 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests I4MAT_HEADER_READ, I4MAT_DATA_READ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string input_filename = "i4mat_05_00020.txt";

        Console.WriteLine();
        Console.WriteLine("TEST04");
        Console.WriteLine("  For data stored in an I4MAT file,");
        Console.WriteLine("  I4MAT_HEADER_READ reads the header");
        Console.WriteLine("  (Information about the dimension of the data)");
        Console.WriteLine("  I4MAT_HEADER_READ reads the data.");

        TableHeader h = typeMethods.i4mat_header_read ( input_filename );

        Console.WriteLine();
        Console.WriteLine("  Read the header of \"" + input_filename + "\".");
        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + h.m);
        Console.WriteLine("  Number of points N  = " + h.n);

        int[] table = typeMethods.i4mat_data_read (  input_filename, h.m, h.n );

        Console.WriteLine();
        Console.WriteLine("  Read the data in \"" + input_filename + "\".");
            
        typeMethods.i4mat_print_some ( h.m, h.n, table, 1, 1, 5, 5, 
            "  5x5 portion of data read from file:" );
    }


    private static void test05 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests R8MAT_UNIFORM_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 2;
        const int N = 10;

        int seed = 123456789;

        Console.WriteLine();
        Console.WriteLine("TEST05");
        Console.WriteLine("  R8MAT_UNIFORM_01 sets a random R8MAT.");
        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + M);
        Console.WriteLine("  Number of points N  = " + N);

        double[] table = UniformRNG.r8mat_uniform_01 ( M, N, ref seed );

        typeMethods.r8mat_print_some ( M, N, table, 1, 1, 5, 10, 
            "  5x10 portion of random real table dataset:" );

    }


    private static void test06 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests I4MAT_BORDER_CUT and I4MAT_BORDER_ADD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 6;
        const int n = 4;

        Console.WriteLine();
        Console.WriteLine("TEST06");
        Console.WriteLine("  I4MAT_BORDER_CUT cuts off the border;");
        Console.WriteLine("  I4MAT_BORDER_ADD adds a zero border;");
        Console.WriteLine();
        Console.WriteLine("  Spatial dimension M = " + m);
        Console.WriteLine("  Number of points N  = " + n);

        int[] table = typeMethods.i4mat_indicator_new ( m, n );

        typeMethods.i4mat_print ( m, n, table, "  Initial dataset:" );

        int[] table2 = typeMethods.i4mat_border_cut ( m, n, table );

        typeMethods.i4mat_print ( m-2, n-2, table2, "  'Cut' dataset:" );

        int[] table3 = typeMethods.i4mat_border_add ( m - 2, n - 2, table2 );

        typeMethods.i4mat_print ( m, n, table3, "  'Added' dataset:" );
    }
}