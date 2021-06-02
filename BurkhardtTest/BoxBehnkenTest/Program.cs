using System;
using Burkardt.BoxBehnkenNS;
using Burkardt.Table;
using Burkardt.Types;

namespace Burkardt.BoxBehnkenTest
{
    class Program
    {
        static void Main(string[] args)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for BOX_BEHNKEN_TEST.
            //
            //  Discussion:
            //
            //    BOX_BEHNKEN_TEST tests the BOX_BEHNKEN library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                Console.WriteLine();
                Console.WriteLine("BOX_BEHNKEN_TEST");
                Console.WriteLine("  Test the BOX_BEHNKEN library.");

                test01();
                test02();

                Console.WriteLine();
                Console.WriteLine("BOX_BEHNKEN_TEST");
                Console.WriteLine("  Normal end of execution.");
                Console.WriteLine();
            }
            
            
            static void test01 ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests BOX_BEHNKEN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                int dim_num = 3;
                double[] range = new [] {
                    0.0, 10.0,  5.0,
                    1.0, 11.0, 15.0 };

                Console.WriteLine();
                Console.WriteLine("TEST01");
                Console.WriteLine("  BOX_BEHNKEN computes a Box-Behnken dataset.");

                typeMethods.r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

                int x_num = BoxBehnken.box_behnken_size ( dim_num );

                Console.WriteLine();
                Console.WriteLine("  For dimension DIM_NUM = " + dim_num);
                Console.WriteLine("  the Box-Behnken design is of size " + x_num);

                double[] x = BoxBehnken.box_behnken ( dim_num, x_num, range );

                typeMethods.r8mat_transpose_print ( dim_num, x_num, x, "  The Box-Behnken design:" );
            }
            
            
            static void test02()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests R8MAT_WRITE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            {
                int dim_num = 4;
                string file_out_name = "box_behnken_04_33.txt";
                double[] range = {
                    0.0, 0.0, 0.0, 0.0,
                    1.0, 1.0, 1.0, 1.0 };

                Console.WriteLine();
                Console.WriteLine("TEST02");
                Console.WriteLine("  R8MAT_WRITE writes a Box-Behnken dataset");
                Console.WriteLine("  to a file.");

                typeMethods.r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

                int x_num = BoxBehnken.box_behnken_size ( dim_num );

                Console.WriteLine();
                Console.WriteLine("  For dimension DIM_NUM = " + dim_num);
                Console.WriteLine("  the Box-Behnken design is of size " + x_num);

                double[] x = BoxBehnken.box_behnken ( dim_num, x_num, range );

                TableWriter.r8mat_write ( file_out_name, dim_num, x_num, x );

                Console.WriteLine();
                Console.WriteLine("  The data was written to the file \"" + file_out_name + "\".");
            }
        }
    }
}
