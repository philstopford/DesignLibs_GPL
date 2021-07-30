using System;
using Burkardt.TOMS;
using Burkardt.Types;

namespace TOMS112Test
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TOMS112_TEST.
            //
            //  Discussion:
            //
            //    TOMS112_TEST tests the TOMS112 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool inside;
            int n = 5;
            int test;
            int test_num = 4;
            double[] x = { 0.0,  1.0,  2.0,  1.0,  0.0 };
            double x0;
            double[] x0_test = { 1.0, 3.0, 0.0, 0.5 };
            double[] y = { 0.0, 0.0, 1.0, 2.0, 2.0 };
            double y0;
            double[] y0_test = { 1.0, 4.0, 2.0, -0.25 };

            Console.WriteLine("");
            Console.WriteLine("TOMS112_TEST");
            Console.WriteLine("  POINT_IN_POLYGON determines if a point is in a polygon.");

            typeMethods.r8vec2_print ( n, x, y, "  The polygon vertices:" );

            Console.WriteLine("");
            Console.WriteLine("        Px        Py  Inside");
            Console.WriteLine("");

            for ( test = 0; test < test_num; test++ )
            {
                x0 = x0_test[test];
                y0 = y0_test[test];
 
                inside = TOMS.point_in_polygon ( n, x, y, x0, y0 );

                Console.WriteLine("  " + x0.ToString().PadLeft(8)
                    + "  " + y0.ToString().PadLeft(8)
                    + "       " + inside + "");
            }
            //
            //  Terminate.  
            //
            Console.WriteLine("");
            Console.WriteLine("TOMS112_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}