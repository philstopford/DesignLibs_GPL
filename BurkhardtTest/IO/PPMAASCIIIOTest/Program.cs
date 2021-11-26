using System;
using Burkardt.IO;

namespace PPMAASCIIIOTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PPMA_IO_TEST.
        //
        //  Discussion:
        //
        //    PPMA_IO_TEST tests the PPMA_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PPMA_IO_TEST:");
        Console.WriteLine("  Test the PPMA_IO library.");

        bool error = test01();

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST01 terminated with an error.");
                return;
        }

        error = test02();

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PPMA_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST02 terminated with an error.");
                return;
        }

        Console.WriteLine("");
        Console.WriteLine("PPMA_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static bool test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PPMA_EXAMPLE and PPMA_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string file_out_name = "test01.ascii.ppm";
        const int xsize = 300;
        const int ysize = 300;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  PPMA_EXAMPLE sets up PPM data.");
        Console.WriteLine("  PPMA_WRITE writes an ASCII PPM file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file \"" + file_out_name + "\".");

        int[] r = new int[xsize * ysize];
        int[] g = new int[xsize * ysize];
        int[] b = new int[xsize * ysize];

        bool error = PPM_ASCII.ppma_example(xsize, ysize, ref r, ref g, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PPMA_EXAMPLE failed!");
                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("  PPMA_EXAMPLE has set up the data.");

        error = PPM_ASCII.ppma_write(file_out_name, xsize, ysize, r, g, b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PPMA_WRITE failed!");
                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("  PPMA_WRITE was successful.");
        //
        //  Now have PPMA_READ_TEST look at the file we think we created.
        //
        error = PPM_ASCII.ppma_read_test(file_out_name);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PPMA_READ_TEST failed to read the file we wrote!");
                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("  PPMA_READ_TEST was able to read the file we wrote.");

        return false;
    }

    private static bool test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PPMA_READ.
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
        int[] b = new int[1];
        const string file_in_name = "test02.ascii.ppm";
        int[] g = new int[1];
        int k;
        int[] r = new int[1];
        int rgb_max = 0;
        int xsize = 0;
        int ysize = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PPMA_READ reads the header and data of a PPMA file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file \"" + file_in_name + "\".");
        //
        //  Create a data file to read.
        //
        bool error = PPM_ASCII.ppma_write_test(file_in_name);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST02");
                Console.WriteLine("  PPMA_WRITE_TEST failed!");
                return true;
        }

        Console.WriteLine("");
        Console.WriteLine("  PPMA_WRITE_TEST created some test data.");
        //
        //  Now have PPMA_READ try to read it.
        //
        error = PPM_ASCII.ppma_read(file_in_name, ref xsize, ref ysize, ref rgb_max, ref r, ref g, ref b);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("TEST02");
                Console.WriteLine("  PPMA_READ failed!");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  PPMA_READ read the test data successfully.");
        Console.WriteLine("");
        Console.WriteLine("  Ten sample values:");
        Console.WriteLine("");
        for (k = 0; k < 10; k++)
        {
            int i = ((9 - k) * 0 + k * (xsize - 1)) / 9;
            int j = ((9 - k) * 0 + k * (ysize - 1)) / 9;
            Console.WriteLine(i.ToString().PadLeft(4) + "  "
                                                      + j.ToString().PadLeft(4) + "  "
                                                      + r[i * ysize + j].ToString().PadLeft(4) + "  "
                                                      + g[i * ysize + j].ToString().PadLeft(4) + "  "
                                                      + g[i * ysize + j].ToString().PadLeft(4) + "");
        }

        return error;
    }
}