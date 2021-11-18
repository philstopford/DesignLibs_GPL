using System;
using Burkardt.IO;

namespace PGMAIOTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PGMA_IO_TEST.
        //
        //  Discussion:
        //
        //    PGMA_IO_TEST tests the PGMA_IO library.
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
        Console.WriteLine("PGMA_IO_TEST:");
        Console.WriteLine("  Test the PGMA_IO library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("PGMA_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests PGMA_EXAMPLE, PGMA_WRITE.
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
        string output_name = "pgma_io_test01.ascii.pgm";
        int[] g;
        int xsize = 300;
        int ysize = 300;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  PGMA_EXAMPLE sets up ASCII PGM data.");
        Console.WriteLine("  PGMA_WRITE writes an ASCII PGM file.");
        Console.WriteLine("");
        Console.WriteLine("  Writing the file \"" + output_name + "\".");

        g = new int[xsize * ysize];

        PGMA.pgma_example(xsize, ysize, ref g);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_EXAMPLE has set up the data.");

        PGMA.pgma_write(output_name, xsize, ysize, g);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_WRITE was successful.");

        //
        //  Now have PGMA_READ_TEST look at the file we think we created.
        //
        PGMA.pgma_read_test(output_name);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_READ_TEST was able to read our file.");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests PGMA_READ_HEADER, PGMA_READ_DATA.
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
        string input_name = "pgma_io_test02.ascii.pgm";
        int[] g = new int[1];
        int i;
        int j;
        int k;
        int maxg = 0;
        int xsize = 0;
        int ysize = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  PGMA_READ reads the header and data of an ASCII PGM file.");
        Console.WriteLine("");
        Console.WriteLine("  Reading the file \"" + input_name + "\".");
        //
        //  Create a data file to read.
        //
        PGMA.pgma_write_test(input_name);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_WRITE_TEST created some test data.");
        //
        //  Now have PGMA_READ try to read it.
        //
        PGMA.pgma_read(input_name, ref xsize, ref ysize, ref maxg, ref g);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_READ read the test data successfully.");

        Console.WriteLine("");
        Console.WriteLine("  Sample data:");
        Console.WriteLine("");
        for (k = 0; k <= 9; k++)
        {
            i = ((9 - k) * 0 + k * (xsize - 1)) / 9;
            j = ((9 - k) * 0 + k * (ysize - 1)) / 9;
            Console.WriteLine(i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                      + g[i * ysize + j].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests PGMA_WRITE.
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
        int NGRAY = 11;

        string output_name = "pgma_io_test03.ascii.pgm";
        int[] g;
        double[] gray =
        {
            0.000, 0.291, 0.434, 0.540, 0.629,
            0.706, 0.774, 0.837, 0.895, 0.949,
            1.000
        };
        int i;
        int j;
        int k;
        int xsize = 300;
        int ysize = 300;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  PGMA_WRITE writes an ASCII PGM file.");
        Console.WriteLine("");
        Console.WriteLine("  In this example, we make a sort of grayscale");
        Console.WriteLine("  checkerboard.");

        g = new int[xsize * ysize];

        for (i = 0; i < xsize; i++)
        {
            for (j = 0; j < ysize; j++)
            {
                k = (i + j) * NGRAY / Math.Min(xsize, ysize);
                k %= NGRAY;
                g[i * ysize + j] = (int) (255.0E+00 * gray[k]);
            }
        }

        Console.WriteLine("  Writing the file \"" + output_name + "\".");

        PGMA.pgma_write(output_name, xsize, ysize, g);

        Console.WriteLine("");
        Console.WriteLine("  PGMA_WRITE was successful.");
    }
}