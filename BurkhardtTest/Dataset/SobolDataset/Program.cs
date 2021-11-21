using System;
using Burkardt.Sobol;
using Burkardt.Types;

namespace SobolDataset;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SOBOL_DATASET.
        //
        //  Discussion:
        //
        //    SOBOL_DATASET generates a Sobol dataset and writes it to a file.
        //
        //  Usage:
        //
        //    sobol_dataset m n skip
        //
        //    where
        //
        //    * M, the spatial dimension,
        //    * N, the number of points to generate,
        //    * SKIP, the number of initial points to skip.
        //
        //    creates an M by N Sobol dataset and writes it to the
        //    file "sobol_M_N.txt".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 December 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        string m_ostring;
        int n;
        string n_ostring;
        string output_filename;
        double[] r;
        int skip;


        Console.WriteLine("");
        Console.WriteLine("SOBOL_DATASET");
        Console.WriteLine("");
        Console.WriteLine("  Generate a Sobol quasirandom dataset.");
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
        //  Get SKIP.
        //
        try
        {
            skip = Convert.ToInt32(args[2]);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("  Enter the value of SKIP");
            skip = Convert.ToInt32(Console.ReadLine());
        }

        Console.WriteLine("  SKIP is = " + skip + "");
        //
        //  Compute the data.
        //
        r = SobolSampler.i8_sobol_generate(m, n, skip);
        //
        //  Write it to a file.
        //
        m_ostring = m.ToString();
        n_ostring = n.ToString();

        output_filename = "sobol_" + m_ostring + "_"
                          + n_ostring + ".txt";

        typeMethods.r8mat_write(output_filename, m, n, r);

        Console.WriteLine("");
        Console.WriteLine("  The data was written to the file \""
                          + output_filename + "\".");

        Console.WriteLine("");
        Console.WriteLine("SOBOL_DATASET:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }
}