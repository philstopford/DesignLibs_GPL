using System;
using Burkardt;
using Burkardt.IO;

namespace PGMBIOTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN tests the PGMB_IO_TEST library.
            //
            //  Discussion:
            //
            //    PGMB_IO_TEST tests the PGMB_IO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool error;

            Console.WriteLine("");
            Console.WriteLine("PGMB_IO_TEST:");
            Console.WriteLine("  Test the PGMB_IO library.");

            error = test01();

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST01 terminated with an error.");
                return;
            }

            error = test02();

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("PGMB_IO_TEST - Fatal error!");
                Console.WriteLine("  TEST02 terminated with an error.");
                return;
            }

            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("PGMB_IO_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static bool test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests PGMB_EXAMPLE, PGMB_WRITE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool error;
            string file_out_name = "pgmb_io_test01.pgm";
            int[] g;
            int i;
            int indexg;
            int j;
            int maxg;
            int xsize = 300;
            int ysize = 300;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  PGMB_EXAMPLE sets up PGMB data.");
            Console.WriteLine("  PGMB_WRITE writes a PGMB file.");
            Console.WriteLine("");
            Console.WriteLine("  Writing the file \"" + file_out_name + "\".");

            g = new int[xsize * ysize];

            error = PGMB.pgmb_example(xsize, ysize, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PGMB_EXAMPLE failed!");
                return error;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_EXAMPLE has set up the data.");
            }

            maxg = 0;
            indexg = 0;

            for (i = 0; i < xsize; i++)
            {
                for (j = 0; j < ysize; j++)
                {
                    if (maxg < g[indexg])
                    {
                        maxg = g[indexg];
                    }

                    indexg = indexg + 1;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Gray scale data has maximum value " + (int)maxg + "");

            error = PGMB.pgmb_write(file_out_name, xsize, ysize, g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("TEST01 - Fatal error!");
                Console.WriteLine("  PGMB_WRITE failed!");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_WRITE was successful.");
            }

            //
            //  Now have PGMB_READ_TEST look at the file we think we created.
            //
            error = PGMB.pgmb_read_test(file_out_name);

            return error;
        }

        static bool test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests PGMB_READ_HEADER, PGMB_READ_DATA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool error;
            string file_in_name = "pgmb_io_test02.pgm";
            int[] g = null;
            int maxg = 0;
            int xsize = 0;
            int ysize = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  PGMB_READ reads the header and data of a PGMB file.");
            Console.WriteLine("");
            Console.WriteLine("  Reading the file \"" + file_in_name + "\".");
            //
            //  Create a data file to read.
            //
            error = PGMB.pgmb_write_test(file_in_name);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_WRITE_TEST failed!");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_WRITE_TEST created some test data.");
            }

            //
            //  Now have PGMB_READ try to read it.
            //
            error = PGMB.pgmb_read(file_in_name, ref xsize, ref ysize, ref maxg, ref g);

            if (error)
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_READ failed!");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  PGMB_READ read the test data successfully.");
            }

            return error;
        }
    }
}