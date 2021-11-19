using System;
using Burkardt.IO;

namespace PBMAIOTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PBMA_IO_TEST.
        //
        //  Discussion:
        //
        //    PBMA_IO_TEST tests the PBMA_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PBMA_IO_TEST:");
        Console.WriteLine("  Test the PBMA_IO library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("PBMA_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PBMA_EXAMPLE, PBMA_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b;
        string file_out_name = "pbma_io_prb_01.ascii.pbm";
        int xsize = 300;
        int ysize = 300;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  PBMA_EXAMPLE sets up ASCII PBM data.");
        Console.WriteLine("  PBMA_WRITE writes an ASCII PBM file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file \"" + file_out_name + "\".");

        b = new int[xsize * ysize];

        PBMA.pbma_example(xsize, ysize, ref b);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_EXAMPLE has set up the data.");

        PBMA.pbma_write(file_out_name, xsize, ysize, b);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_WRITE was successful.");

        //
        //  Now have PBMA_READ_TEST look at the file we think we created.
        //
        PBMA.pbma_read_test(file_out_name);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_READ_TEST was able to read the file.");

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PBMA_READ_HEADER, PBMA_READ_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b = null;
        string file_in_name = "pbma_io_prb_02.ascii.pbm";
        int i;
        int j;
        int k;
        int xsize = 0;
        int ysize = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PBMA_READ reads the header and data of an ASCII PBM file.");
        //
        //  Create a data file to read.
        //
        PBMA.pbma_write_test(file_in_name);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_WRITE_TEST created some test data.");
        //
        //  Now have PBMA_READ try to read it.
        //
        PBMA.pbma_read(file_in_name, ref xsize, ref ysize, ref b);

        Console.WriteLine("");
        Console.WriteLine("  PBMA_READ was able to read the file we created.");
        Console.WriteLine("");
        Console.WriteLine("  Sample data:");
        Console.WriteLine("");

        for (k = 0; k <= 29; k++)
        {
            i = ((29 - k) * 0 + k * (xsize - 1)) / 29;
            j = ((29 - k) * 0 + k * (ysize - 1)) / 29;
            Console.WriteLine(i.ToString().PadLeft(4) + "  "
                                                      + j.ToString().PadLeft(4) + "  "
                                                      + b[i * ysize + j].ToString().PadLeft(6) + "");
        }
    }
}