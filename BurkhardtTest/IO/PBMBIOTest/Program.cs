using System;
using Burkardt.IO;

namespace PBMBIOTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PBMB_IO_TEST.
        //
        //  Discussion:
        //
        //    PBMB_IO_TEST tests the PBMB_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool error;

        Console.WriteLine("");
        Console.WriteLine("PBMB_IO_TEST:");
        Console.WriteLine("  Test the PBMB_IO library.");

        error = test01();

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PBMB_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST01 terminated with an error.");
                return;
        }

        error = test02();

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PBMB_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST02 terminated with an error.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("PBMB_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static bool test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PBMB_EXAMPLE, PBMB_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b = null;
        bool error;
        string file_out_name = "pbmb_io_prb_01.pbm";
        int xsize = 300;
        int ysize = 300;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  PBMB_EXAMPLE sets up PBMB data.");
        Console.WriteLine("  PBMB_WRITE writes a PBMB file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file \"" + file_out_name + "\".");

        b = new int[xsize * ysize];

        error = PBMB.pbmb_example(xsize, ysize, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PBMB_EXAMPLE failed!");
                return error;
        }

        Console.WriteLine("");
        Console.WriteLine("  PBMB_EXAMPLE has set up the data.");

        error = PBMB.pbmb_write(file_out_name, xsize, ysize, b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PBMB_WRITE failed!");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  PBMB_WRITE was successful.");
                break;
        }

        //
        //  Now have PBMB_READ_TEST look at the file we think we created.
        //
        error = PBMB.pbmb_read_test(file_out_name);

        return error;
    }

    private static bool test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PBMB_READ_HEADER, PBMB_READ_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b = null;
        bool error;
        string file_in_name = "pbmb_io_prb_02.pbm";
        int xsize = 0;
        int ysize = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PBMB_READ reads the header and data of a PBMB file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file \"" + file_in_name + "\".");
        //
        //  Create a data file to read.
        //
        error = PBMB.pbmb_write_test(file_in_name);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  PBMB_WRITE_TEST failed!");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  PBMB_WRITE_TEST created some test data.");
                break;
        }

        //
        //  Now have PBMB_READ try to read it.
        //
        error = PBMB.pbmb_read(file_in_name, ref xsize, ref ysize, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  PBMB_READ failed!");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  PBMB_READ read the test data successfully.");
                break;
        }

        return error;
    }
}