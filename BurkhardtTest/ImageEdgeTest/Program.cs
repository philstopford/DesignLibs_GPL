using System;
using System.IO;
using Burkardt.Function;
using Burkardt.IO;
using Burkardt.Types;

namespace ImageEdgeTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for IMAGE_EDGE_TEST.
        //
        //  Discussion:
        //
        //    IMAGE_EDGE_TEST tests the IMAGE_EDGE library.
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
        int[] e;
        int[] g;
        int[] g_histo;
        int g_max = 0;
        int i;
        string input_filename = "coins.ascii.pgm";
        string[] input_unit;
        int m = 0;
        int n = 0;
        string output_filename = "coin_edges.ascii.pbm";

        Console.WriteLine("");
        Console.WriteLine("IMAGE_EDGE_TEST");
        Console.WriteLine("  Test the IMAGE_EDGE library.");
        Console.WriteLine("");
        Console.WriteLine("  Demonstrate the NEWS stencil for edge detection");
        Console.WriteLine("  in images.");

        Console.WriteLine("");
        Console.WriteLine("  The input file is \"" + input_filename + "\".");

        try
        {
            input_unit = File.ReadAllLines(input_filename);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("IMAGE_EDGE_TEST - Fatal error!");
            Console.WriteLine("  Could not open the file \"" + input_filename + "\"");
            return;
        }

        int inputIndex = 0;
        PGMA.pgma_read_header(input_unit, ref inputIndex, ref m, ref n, ref g_max);

        Console.WriteLine("");
        Console.WriteLine("  Number of rows =          " + m + "");
        Console.WriteLine("  Number of columns =       " + n + "");
        Console.WriteLine("  Maximum pixel intensity = " + g_max + "");

        g = new int[m * n];

        PGMA.pgma_read_data(input_unit, ref inputIndex, m, n, ref g);


        g_histo = typeMethods.i4mat_histogram(m, n, g, 255);

        Console.WriteLine("");
        Console.WriteLine(" Gray     Count");
        Console.WriteLine("");
        for (i = 0; i <= 255; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                   + "  " + g_histo[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        e = NEWS.news(m, n, g);
        //
        //  Write the edge information as a portable BIT map (0/1).
        //
        PBMA.pbma_write(output_filename, m, n, e);

        Console.WriteLine("");
        Console.WriteLine("  Wrote edge information to \"" + output_filename + "\".");

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("IMAGE_EDGE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}