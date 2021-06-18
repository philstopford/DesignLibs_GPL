using System;
using Burkardt.CircleNS;
using Burkardt.Types;

namespace CircleArcGridTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for CIRCLE_ARC_GRID_TEST.
            //
            //  Discussion:
            //
            //    CIRCLE_ARC_GRID_TEST tests the CIRCLE_ARC_GRID library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("CIRCLE_ARC_GRID_TEST");

            Console.WriteLine("  Test the CIRCLE_ARC_GRID library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("CIRCLE_ARC_GRID_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 demonstrates the use of CIRCLE_ARC_GRID.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a = new double[2];
            double[] c = new double[2];
            string filename = "arc.txt";
            int n;
            double r;
            double[] xy;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Compute points along a 90 degree arc");

            r = 2.0;
            c[0] = 5.0;
            c[1] = 5.0;
            a[0] = 0.0;
            a[1] = 90.0;
            n = 10;
            //
            //  Echo the input.
            //
            Console.WriteLine("");
            Console.WriteLine("  Radius =           " + r + "");
            Console.WriteLine("  Center =           " + c[0] + "  " + c[1] + "");
            Console.WriteLine("  Angle 1 =          " + a[0] + "");
            Console.WriteLine("  Angle 2 =          " + a[1] + "");
            ;
            Console.WriteLine("  Number of points = " + n + "");
            //
            //  Compute the data.
            //
            xy = Circle.circle_arc_grid(r, c, a, n);
            //
            //  Print a little of the data.
            //
            typeMethods.r82vec_print_part(n, xy, 5, "  A few of the points:");
            //
            //  Write the data.
            //
            typeMethods.r8mat_write(filename, 2, n, xy);
            Console.WriteLine("");
            Console.WriteLine("  Data written to \"" + filename + "\".");
        }
    }
}