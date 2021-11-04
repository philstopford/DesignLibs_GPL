using System;
using Burkardt.IHSNS;
using Burkardt.Types;

namespace IHSDatasetTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for IHS_DATASET.
            //
            //  Discussion:
            //
            //    IHS_DATASET generates an IHS dataset and writes it to a file.
            //
            //    This program is meant to be used interactively.  It's also
            //    possible to prepare a simple input file beforehand and use it
            //    in batch mode.
            //
            //    The program requests input values from the user:
            //
            //    * M, the spatial dimension,
            //    * N, the number of points to generate,
            //    * D, the duplication factor,
            //    * SEED, a seed for the random number generator.
            //
            //    The program generates the data, writes it to the file
            //
            //      ihs_M_N_D_SEED.txt
            //
            //    where "M" and "N" are the numeric values specified by the user,
            //    and then asks the user for more input.   To indicate that no further
            //    computations are desired, it is enough to input a nonsensical
            //    value, such as -1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            int i;
            int[] i4_ihs;
            int j;
            int m;
            int n;
            string output_filename = "";
            int seed;
            int seed_init;
            double[] r8_ihs;
            //
            //  Print introduction and options.
            //

            Console.WriteLine("");
            Console.WriteLine("IHS_DATASET");
            
            Console.WriteLine("  Generate an IHS dataset.");
            Console.WriteLine("");
            Console.WriteLine("  This program is meant to be used interactively.");
            Console.WriteLine("  It is also possible to prepare a simple input");
            Console.WriteLine("  file beforehand and use it in batch mode.");
            Console.WriteLine("");
            Console.WriteLine("  The program requests input values from the user:");
            Console.WriteLine("");
            Console.WriteLine("  * M, the spatial dimension,");
            Console.WriteLine("  * N, the number of points to generate,");
            Console.WriteLine("  * D, the duplication factor,");
            Console.WriteLine("  * SEED, a seed for the random number generator.");
            Console.WriteLine("");
            Console.WriteLine("  The program generates the data, writes it to the file");
            Console.WriteLine("");
            Console.WriteLine("    ihs_M_N_D_SEED.txt");
            Console.WriteLine("");
            Console.WriteLine("  where \"M\" and \"N\" are the numeric values specified");
            Console.WriteLine("  by the user, and then asks the user for more input.");
            Console.WriteLine("");
            Console.WriteLine("  To indicate that no further computations are");
            Console.WriteLine("  desired, it is enough to input a nonsensical value,");
            Console.WriteLine("  such as -1.");

            for (;;)
            {
                Console.WriteLine("  *");
                Console.WriteLine(" *");
                Console.WriteLine("*  Ready to generate a new dataset:");
                Console.WriteLine(" *");
                Console.WriteLine("  *");
                Console.WriteLine("  Enter M, the spatial dimension:");
                Console.WriteLine("  (2 is a small typical value).");
                Console.WriteLine("  (0 or any negative value terminates execution).");

                try
                {
                    m = Convert.ToInt32(Console.ReadLine());

                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET - Fatal error!");
                    Console.WriteLine("  An I/O error occurred while trying to read M.");
                    Console.WriteLine("  Abnormal end of execution.");
                    return;
                }

                Console.WriteLine("  User input M = " + m + "");

                if (m < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET");
                    Console.WriteLine("  The input value of M = " + m + "");
                    Console.WriteLine("  is interpreted as a request for termination.");
                    Console.WriteLine("  Normal end of execution.");
                    return;
                }

                Console.WriteLine("");
                Console.WriteLine("  Enter N, the number of points to generate:");
                Console.WriteLine("  (10 is a small typical value).");
                Console.WriteLine("  (0 or any negative value terminates execution).");

                try
                {
                    n = Convert.ToInt32(Console.ReadLine());

                }
                catch (Exception e)
                {

                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET - Fatal error!");
                    Console.WriteLine("  An I/O error occurred while trying to read N.");
                    Console.WriteLine("  Abnormal end of execution.");
                    break;
                }

                Console.WriteLine("  User input N = " + n + "");

                if (n < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET");
                    Console.WriteLine("  The input value of N = " + n + "");
                    Console.WriteLine("  is interpreted as a request for termination.");
                    Console.WriteLine("  Normal end of execution.");
                    return;
                }


                Console.WriteLine("");
                Console.WriteLine("  Enter D, the duplication factor:");
                Console.WriteLine("  This must be at least 1, but not too large.");
                Console.WriteLine("  (5 is a typical value).");
                Console.WriteLine("  (0 or any negative value terminates execution).");

                try
                {
                    d = Convert.ToInt32(Console.ReadLine());

                }
                catch (Exception e)
                {

                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET - Fatal error!");
                    Console.WriteLine("  An I/O error occurred while trying to read D.");
                    Console.WriteLine("  Abnormal end of execution.");
                    return;
                }

                Console.WriteLine("  User input D = " + d + "");

                if (d < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET");
                    Console.WriteLine("  The input value of D = " + d + "");
                    Console.WriteLine("  is interpreted as a request for termination.");
                    Console.WriteLine("  Normal end of execution.");
                    break;
                }

                Console.WriteLine("");
                Console.WriteLine("  Enter SEED, a seed for the random number generator:");
                Console.WriteLine("  (Try '123456789' if you have no preference.)");
                Console.WriteLine("  (0 indicates you want a seed to be chosen for you.)");
                Console.WriteLine("  (Any negative value terminates execution).");

                try
                {
                    seed = Convert.ToInt32(Console.ReadLine());

                }
                catch (Exception e)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET - Fatal error!");
                    Console.WriteLine("  An I/O error occurred while trying to read SEED.");
                    Console.WriteLine("  Abnormal end of execution.");
                    return;
                }

                Console.WriteLine("  User input SEED = " + seed + "");

                if (seed < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("IHS_DATASET");
                    Console.WriteLine("  The input value of SEED = " + seed + "");
                    Console.WriteLine("  is interpreted as a request for termination.");
                    Console.WriteLine("  Normal end of execution.");
                    return;
                }

                if (seed == 0)
                {
                    seed = entropyRNG.RNG.nextint();
                    Console.WriteLine("");
                    Console.WriteLine("  The actual value of SEED will be = " + seed + "");
                }

                //
                //  Compute the integer data.
                //
                seed_init = seed;
                i4_ihs = IHS.ihs(m, n, d, ref seed);
                //
                //  Convert from integer to real.
                //
                r8_ihs = new double[m * n];
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        r8_ihs[i + j * m] = (double) (2 * i4_ihs[i + j * m] - 1)
                                            / (double) (2 * n);
                    }
                }

                //
                //  Write the data to a file.
                //
                output_filename += "ihs_"
                                   + m + "_"
                                   + n + "_"
                                   + d + "_"
                                   + seed_init + ".txt";

                typeMethods.r8mat_write(output_filename, m, n, r8_ihs);

                Console.WriteLine("");
                Console.WriteLine("  The data was written to the file \""
                                  + output_filename + "\"");
            }

            Console.WriteLine("");
        }
    }
}