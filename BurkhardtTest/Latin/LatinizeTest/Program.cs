using System;
using Burkardt.Latinizer;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.LatinizeTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LATINIZE_TEST.
            //
            //  Discussion:
            //
            //    LATINIZE_TEST tests the LATINIZE library.
            //
            //    The dataset is presumed to be an M by N array of real numbers,
            //    where M is the spatial dimension, and N is the number of sample points.
            //
            //    The dataset is presumed to be stored in a file, with N records,
            //    one per each sample point.  (Comment records may be included, 
            //    which begin with '#'.)
            //
            //    The program reads the data file, "latinizes" the data, and writes
            //    the latinized data to a new file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            Console.WriteLine();
            Console.WriteLine("LATINIZE_TEST");
            Console.WriteLine("  Test the LATINIZE library.");
            Console.WriteLine();
            Console.WriteLine("  Read a dataset of N points in M dimensions,");
            Console.WriteLine("  modify it into a Latin hypercube,");
            Console.WriteLine("  write the modified dataset to a file.");

            test01 ( "cvt_02_00010.txt" );
            test01 ( "cvt_03_00007.txt" );
            test01 ( "cvt_03_00056.txt" );
            test01 ( "cvt_07_00100.txt" );

            Console.WriteLine();
            Console.WriteLine("LATINIZE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }
        
        static void test01 ( string input_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATINIZE_TEST tests the LATINIZE routines.
        //
        //  Discussion:
        //
        //    The dataset is presumed to be an M by N array of real numbers,
        //    where M is the spatial dimension, and N is the number of sample points.
        //
        //    The dataset is presumed to be stored in a file, with N records,
        //    one per each sample point.  (Comment records may be included, 
        //    which begin with '#'.)
        //
        //    The program reads the data file, "latinizes" the data, and writes
        //    the latinized data to a new file.
        //
        //  Modified:
        //
        //    08 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] table;
            //
            //  Need to create the output file name from the input filename.
            //
            string output_filename = TableWriter.file_name_ext_swap ( input_filename, "latin.txt" );

            TableHeader header = TableReader.r8mat_header_read ( input_filename);

            Console.WriteLine();
            Console.WriteLine("  Read the header of \"" + input_filename + "\".");
            Console.WriteLine();
            Console.WriteLine("  Spatial dimension M = " + header.m );
            Console.WriteLine("  Number of points N  = " + header.n);

            table = TableReader.r8mat_data_read ( input_filename, header.m, header.n );

            Console.WriteLine();
            Console.WriteLine("  Read the data in \"" + input_filename + "\".");

            typeMethods.r8mat_transpose_print_some ( header.m, header.n, table, 1, 1, 5, 5, 
                "  Small portion of data read from file:" );

            double[] latTable = Latinize.r8mat_latinize ( header.m, header.n, table );

            Console.WriteLine();
            Console.WriteLine("  Latinized the data.\n");

            typeMethods.r8mat_transpose_print_some ( header.m, header.n, latTable, 1, 1, 5, 5, 
                "  Small portion of Latinized data:" );
            //
            //  Write the data to a file.
            //
            TableWriter.r8mat_write ( output_filename, header.m, header.n, latTable );

            Console.WriteLine();
            Console.WriteLine("  Wrote the latinized data to \"" + output_filename + "\".");
            
            return;
        }        
    }
}