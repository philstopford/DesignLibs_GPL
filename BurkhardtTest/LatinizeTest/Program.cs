using System;
using Latinizer;

namespace LatinizeTest
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
    }
}