using System;
using Burkardt.Sequence;
using Burkardt.Types;

namespace VanderCorputDataset
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for VAN_DER_CORPUT_DATASET.
            //
            //  Discussion:
            //
            //    VAN_DER_CORPUT_DATASET generates a van der Corput dataset.
            //
            //    This program is meant to be used interactively.  It's also
            //    possible to prepare a simple input file beforehand and use it
            //    in batch mode.
            //
            //  Usage:
            //
            //    van_der_corput_dataset base seed n
            //
            //    where
            //
            //    * BASE, the base of the sequence;
            //    * SEED, the index of the first element to be computed;
            //    * N, the number of points to generate.
            //
            //    The program generates the data and writes it to the file
            //
            //      van_der_corput_BASE_SEED_N.txt
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int base_;
            string output_filename;
            int n;
            double[] r;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("VAN_DER_CORPUT_DATASET");
            Console.WriteLine("");
            Console.WriteLine("  Generate a van der Corput dataset.");
            Console.WriteLine("");
            Console.WriteLine("  The program requests input values from the user:");
            Console.WriteLine("");
            Console.WriteLine("  * BASE, the base,");
            Console.WriteLine("  * SEED, a positive integer.");
            Console.WriteLine("  * N, the number of points to generate,");
            Console.WriteLine("");
            Console.WriteLine("  The program generates the data, writes it to the file");
            Console.WriteLine("");
            Console.WriteLine("    van_der_corput_BASE_SEED_N.txt");
            Console.WriteLine("");
            //
            //  Get the BASE.
            //
            try
            {
                base_ = Convert.ToInt32(args[0]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter BASE, the van der Corput base,");
                Console.WriteLine("  which is often a prime number:");
                Console.WriteLine("  (Try '2' if you don't have a preference.)");
                base_ = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("");
            Console.WriteLine("  BASE = " + base_ + "");
            //
            //  Get SEED.
            //
            try
            {
                seed = Convert.ToInt32(args[1]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter SEED, which is 0 or a positive integer.");
                Console.WriteLine("  (Try '0' or '1' if you don't have a preference.)");
                seed = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("  SEED = " + seed + "");
            //
            //  Get N.
            //
            try
            {
                n = Convert.ToInt32(args[2]);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Enter N, the number of points to generate:");
                Console.WriteLine("  (Try '25' if you don't have a preference.)");
                n = Convert.ToInt32(Console.ReadLine());
            }

            Console.WriteLine("  N = " + n + "");
            //
            //  Compute the data.
            //
            r = new double[n];

            VanDerCorput.i4_to_van_der_corput_sequence(seed, base_, n, ref r);
            //
            //  Write the data to a file.
            //
            output_filename = "van_der_corput_" + base_ + "_" + seed + "_" + n + ".txt";

            typeMethods.r8mat_write(output_filename, 1, n, r);

            Console.WriteLine("");
            Console.WriteLine("  The data was written to the file \""
                              + output_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("VAN_DER_CORPUT_DATASET:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}