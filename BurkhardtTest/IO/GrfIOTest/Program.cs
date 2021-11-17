using System;
using Burkardt.IO;

namespace GrfIOTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for GRF_IO_TEST.
        //
        //  Discussion:
        //
        //    GRF_IO_TEST tests the GRF_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("GRF_IO_TEST:");
        Console.WriteLine("  Test the GRF_IO library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("GRF_IO_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests GRF_HEADER_WRITE, GRF_DATA_WRITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int edge_num = 0;
        int[] edge_data;
        int[] edge_pointer;
        int node_num = 0;
        string output_filename = "coxeter.grf";
        double[] xy;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  GRF_HEADER_WRITE writes the header of a GRF file.");
        Console.WriteLine("  GRF_DATA_WRITE writes the data of a GRF file.");

        Grf.grf_example_size(ref node_num, ref edge_num);

        Grf.grf_header_print(node_num, edge_num);

        edge_data = new int[edge_num];
        edge_pointer = new int[node_num + 1];
        xy = new double[2 * node_num];

        Grf.grf_example(node_num, edge_num, ref edge_pointer, ref edge_data, ref xy);

        Grf.grf_write(output_filename, node_num, edge_num, edge_pointer,
            edge_data, xy);

        Console.WriteLine("");
        Console.WriteLine("  Wrote the GRF file \"" + output_filename + "\",");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests GRF_HEADER_READ and GRF_DATA_READ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int edge_num = 0;
        string input_filename = "coxeter.grf";
        int[] edge_data;
        int[] edge_pointer;
        int node_num = 0;
        double[] xy;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  GRF_HEADER_READ reads the header of a GRF file.");
        Console.WriteLine("  GRF_DATA_READ reads the data of a GRF file.");

        Console.WriteLine("");
        Console.WriteLine("  Reading the GRF file \"" + input_filename + "\"");

        Grf.grf_header_read(input_filename, ref node_num, ref edge_num);

        Grf.grf_header_print(node_num, edge_num);

        edge_pointer = new int[node_num + 1];
        edge_data = new int[edge_num];
        xy = new double[2 * node_num];

        Grf.grf_data_read(input_filename, node_num, edge_num, ref edge_pointer,
            ref edge_data, ref xy);

        Grf.grf_data_print(node_num, edge_num, edge_pointer, edge_data, xy);
    }
}