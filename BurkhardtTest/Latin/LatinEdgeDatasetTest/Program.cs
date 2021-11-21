using System;
using Burkardt.Latin;
using Burkardt.Types;

namespace LatinEdgeDatasetTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LATIN_EDGE_DATASET.
        //
        //  Discussion:
        //
        //    LATIN_EDGE_DATASET generates a Latin Edge Square dataset 
        //    and writes it to a file.
        //
        //  Usage:
        //
        //    latin_center_dataset m n seed
        //
        //    where
        //
        //    * M, the spatial dimension,
        //    * N, the number of points to generate,
        //    * SEED, the seed, a positive integer.
        //
        //    creates an M by N Latin Edge Square dataset and writes it to the
        //    file "latin_edge_M_N.txt".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        string output_filename;
        double[] r;
        int seed;


        Console.WriteLine("");
        Console.WriteLine("LATIN_EDGE_DATASET");
        Console.WriteLine("");
        Console.WriteLine("  Generate a Latin Edge Square dataset.");
        //
        //  Get the spatial dimension.
        //
        try
        {
            m = Convert.ToInt32(args[0]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of M");
            m = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");
        //
        //  Get the number of points.
        //
        try
        {
            n = Convert.ToInt32(args[1]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the number of points N");
            n = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("  Number of points N = " + n + "");
        //
        //  Get the seed.
        //
        try
        {
            seed = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of SEED");
            seed = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("  The seed is = " + seed + "");
        //
        //  Compute the data.
        //
        r = new double[m * n];

        r = LatinVariants.latin_edge(m, n, ref seed);
        //
        //  Write it to a file.
        //

        output_filename = "latin_edge_" + m + "_"
                          + n + ".txt";

        typeMethods.r8mat_write(output_filename, m, n, r);

        Console.WriteLine("");
        Console.WriteLine("  The data was written to the file \""
                          + output_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("LATIN_EDGE_DATASET:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}