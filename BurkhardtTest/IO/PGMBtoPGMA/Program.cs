using System;
using Burkardt.IO;

namespace PGMBtoPGMA;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PGMB_TO_PGMA.
        //
        //  Discussion:
        //
        //    PGMB_TO_PGMA converts a binary PGM file to an ASCII PGM file.
        //
        //  Usage:
        //
        //    pgmb_to_pgma file.pgmb file.pgma
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    FILE.PGMB is the name of the input binary PGM file to be read.
        //
        //    FILE.PGMA is the name of the output ASCII PGM file to be created.
        //  
    {
        string file_in_name;
        string file_out_name;
        const bool verbose = false;

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PGMB_TO_PGMA:");
                Console.WriteLine("");
                Console.WriteLine("  Convert a binary PGM file to ASCII PGM format.");
                break;
        }

        //
        //  Get the specification for the input file.
        //
        try
        {
            file_in_name = args[0];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("PGMB_TO_PGMA:");
            Console.WriteLine("  Please enter the input PGMB file name:");

            file_in_name = Console.ReadLine();
        }

        //
        //  Get the specification for the output file.
        //
        try
        {
            file_out_name = args[1];
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("PGMB_TO_PGMA:");
            Console.WriteLine("  Please enter the output PGMA file name:");

            file_out_name = Console.ReadLine();
        }

        bool error = PGMB.pgmb_to_pgma(file_in_name, file_out_name);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PGMB_TO_PGMA - Fatal error!");
                Console.WriteLine("  The conversion was not successful.");
                return;
        }

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PGMB_TO_PGMA:");
                Console.WriteLine("  Normal end of execution.");

                Console.WriteLine("");
                break;
        }
    }
}