using System;
using Burkardt.IO;
using Burkardt.Table;
using Burkardt.Types;

namespace TableLatinizeTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TABLE_LATINIZE.
        //
        //  Discussion:
        //
        //    TABLE_LATINIZE is the main routine of a program to "latinize" a dataset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Usage:
        //
        //    table_latinize input_filename
        //
    {
        int i;
        string input_filename;


        Console.WriteLine("");
        Console.WriteLine("LATINIZE");
        Console.WriteLine("");
        Console.WriteLine("  Read a dataset of N points in M dimensions,");
        Console.WriteLine("  modify it into a Latin hypercube,");
        Console.WriteLine("  write the modified dataset to a file.");

        try
        {
            for (i = 1; i < args.Length; i++)
            {
                input_filename = args[i];
                handle(input_filename);
            }
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("TABLE_LATINIZE:");
            Console.WriteLine("  Please enter the name of a file to be analyzed.");

            input_filename = Console.ReadLine();

            handle(input_filename);

        }

        Console.WriteLine("");
        Console.WriteLine("TABLE_LATINIZE");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void handle(string input_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HANDLE handles a single file.
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
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the input file.
        //
        //  Local parameters:
        //
        //    Local, int DIM_NUM, the spatial dimension of the point set.
        //
        //    Local, int N, the number of points.
        //
        //    Local, double Z[DIM_NUM*N], the point set.
        //
        //    Local, int NS, the number of sample points.
        //
    {
        int dim_num;
        int n;
        string output_filename;
        double[] table;
        //
        //  Need to create the output file name from the input filename.
        //
        output_filename = Files.file_name_ext_swap(input_filename, "latin.txt");

        TableHeader h = typeMethods.r8mat_header_read(input_filename);
        dim_num = h.m;
        n = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + input_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
        Console.WriteLine("  Number of points N  = " + n + "");

        table = typeMethods.r8mat_data_read(input_filename, dim_num, n);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + input_filename + "\".");

        typeMethods.r8mat_print_some(dim_num, n, table, 1, 1, 5, 5,
            "  Small portion of data read from file:");

        typeMethods.r8mat_latinize(dim_num, n, ref table);

        Console.WriteLine("");
        Console.WriteLine("  Latinized the data.");

        typeMethods.r8mat_print_some(dim_num, n, table, 1, 1, 5, 5,
            "  Small portion of Latinized data:");

        typeMethods.r8mat_write(output_filename, dim_num, n, table);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the latinized data to \"" + output_filename + "\".");

    }

}