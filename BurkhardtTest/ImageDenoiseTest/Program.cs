using System;
using System.IO;
using Burkardt.Function;
using Burkardt.IO;

namespace ImageDenoiseTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for IMAGE_DENOISE_TEST.
        //
        //  Discussion:
        //
        //    IMAGE_DENOISE_TEST tests the IMAGE_DENOISE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("IMAGE_DENOISE_TEST");
        Console.WriteLine("  Test the IMAGE_DENOISE library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("IMAGE_DENOISE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests GRAY_MEDIAN_NEWS.
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
        int[] g;
        int g_max = 0;
        int[] g2;
        string input_filename = "glassware_noisy.ascii.pgm";
        string[] input_unit;
        int m = 0;
        int n = 0;
        string output_filename = "glassware_median_news.ascii.pgm";

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  GRAY_MEDIAN_NEWS uses a NEWS median filter ");
        Console.WriteLine("  on a noisy grayscale image.");

        Console.WriteLine("");
        Console.WriteLine("  The input file is \"" + input_filename + "\".");

        //
        //  Open the input file and read the data.
        //
        try
        {
            input_unit = File.ReadAllLines(input_filename);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01 - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + input_filename + "\"");
            return;
        }

        int index = 0;

        PGMA.pgma_read_header(input_unit, ref index, ref m, ref n, ref g_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of rows =          " + m + "");
        Console.WriteLine("  Number of columns =       " + n + "");
        Console.WriteLine("  Maximum pixel intensity = " + g_max + "");

        g = new int[m * n];

        PGMA.pgma_read_data(input_unit, ref index, m, n, ref g);

        g2 = NEWS.gray_median_news(m, n, g);
        //
        //  Write the denoised images.
        //
        PGMA.pgma_write(output_filename, m, n, g2);

        Console.WriteLine("");
        Console.WriteLine("  Wrote denoised image to \"" + output_filename + "\".");

    }
}