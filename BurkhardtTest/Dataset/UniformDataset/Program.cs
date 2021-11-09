using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace UniformDataset
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for UNIFORM_DATASET.
            //
            //  Discussion:
            //
            //    UNIFORM_DATASET generates a uniform pseudorandom dataset 
            //    and writes it to a file.
            //
            //  Usage:
            //
            //    uniform_dataset m n seed
            //
            //    where
            //
            //    * M, the spatial dimension,
            //    * N, the number of points to generate,
            //    * SEED, the seed, a positive integer.
            //
            //    creates an M by N uniform random dataset and writes it to the
            //    file "uniform_M_N.txt".
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 December 2009
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
            int seed;

            Console.WriteLine("");
            Console.WriteLine("UNIFORM_DATASET");
            Console.WriteLine("");
            Console.WriteLine("  Generate a uniform pseudorandom dataset.");
            Console.WriteLine("");
            Console.WriteLine("  The program requests input values from the user:");
            Console.WriteLine("");
            Console.WriteLine("  * M, the spatial dimension,");
            Console.WriteLine("  * N, the number of points to generate,");
            Console.WriteLine("  * SEED, a positive integer.");
            Console.WriteLine("");
            Console.WriteLine("  The program generates the data, writes it to the file");
            Console.WriteLine("");
            Console.WriteLine("    uniform_M_N.txt");
            Console.WriteLine("");
            Console.WriteLine("  where \"M\" and \"N\" are the numeric values specified");
            Console.WriteLine("  by the user.");
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
            r = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);
            //
            //  Write it to a file.
            //
            m_ostring = m.ToString();
            n_ostring = n.ToString();

            output_filename = "uniform_" + m_ostring + "_"
                              + n_ostring + ".txt";

            typeMethods.r8mat_write(output_filename, m, n, r);

            Console.WriteLine("");
            Console.WriteLine("  The data was written to the file \""
                              + output_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("UNIFORM_DATASET:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}