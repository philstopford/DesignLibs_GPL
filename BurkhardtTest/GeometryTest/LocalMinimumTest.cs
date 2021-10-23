using System;
using Burkardt.Geometry;

namespace GeometryTest
{
    public static class LocalMinimumTest
    {

        public static void test046()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST046 tests MINABS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double x1;
            double x2;
            double x3;
            double xmin = 0;
            double y1;
            double y2;
            double y3;
            double ymin = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST046");
            Console.WriteLine("  MINABS finds the minimum of a function");
            Console.WriteLine("    F(X) = a * ABS ( X ) + B");
            Console.WriteLine("  within an interval, given three data points.");
            //
            //  Case 1: the three points lie on a straight line.
            //  (XMIN=9,YMIN=2).
            //
            x1 = 14.0;
            y1 = 7.0;

            x2 = 9.0;
            y2 = 2.0;

            x3 = 12.0;
            y3 = 5.0;

            LocalMinimum.minabs(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            Console.WriteLine("");
            Console.WriteLine("  The points lie on a straight line.");
            Console.WriteLine("  XMIN = " + xmin + "  YMIN = " + ymin + "");
            //
            //  Case 2: the three points straddle a minimum.
            //  (XMIN=7, YMIN=2).
            //
            x1 = 3.0;
            y1 = 6.0;

            x2 = 12.0;
            y2 = 7.0;

            x3 = 9.0;
            y3 = 4.0;

            LocalMinimum.minabs(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            Console.WriteLine("");
            Console.WriteLine("  The points straddle a minimum.");
            Console.WriteLine("  XMIN = " + xmin + "  YMIN = " + ymin + "");
            //
            //  Case 3: the three points straddle a maximum.
            //  (XMIN=2, YMIN=5).
            //
            x1 = 11.0;
            y1 = 6.0;

            x2 = 6.0;
            y2 = 9.0;

            x3 = 2.0;
            y3 = 5.0;

            LocalMinimum.minabs(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            Console.WriteLine("");
            Console.WriteLine("  The points straddle a maximum.");
            Console.WriteLine("  XMIN = " + xmin + "  YMIN = " + ymin + "");

        }

        public static void test047()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST047 tests MINQUAD.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool result;
            double x1;
            double x2;
            double x3;
            double xmin = 0;
            double y1;
            double y2;
            double y3;
            double ymin = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST047");
            Console.WriteLine("  MINQUAD finds the minimum of a function");
            Console.WriteLine("    F(X) = A * X^2 + B * X + C");
            Console.WriteLine("  within an interval, given three data points.");
            //
            //  Case 1: a minimum is in the interval.
            //  y = ( x - 1 )**2 + 4
            //
            x1 = 0.0;
            y1 = (x1 - 1.0) * (x1 - 1.0) + 4.0;

            x2 = 2.0;
            y2 = (x2 - 1.0) * (x2 - 1.0) + 4.0;

            x3 = 3.0;
            y3 = (x3 - 1.0) * (x3 - 1.0) + 4.0;

            result = LocalMinimum.minquad(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            if (!result)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning");
                Console.WriteLine("  MINQUAD returned an error code.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The minimum lies in the interval.");
                Console.WriteLine("  X1, Y1 = " + x1 + "  " + y1 + "");
                Console.WriteLine("  X2, Y2 = " + x2 + "  " + y2 + "");
                Console.WriteLine("  X3, Y3 = " + x3 + "  " + y3 + "");
                Console.WriteLine("  XMIN = " + xmin + " YMIN = " + ymin + "");
            }

            //
            //  Case 2: the minimum is to the left of the interval.
            //  y = ( x - 1 )**2 + 4
            //
            x1 = 2.0;
            y1 = (x1 - 1.0) * (x1 - 1.0) + 4.0;

            x2 = 4.0;
            y2 = (x2 - 1.0) * (x2 - 1.0) + 4.0;

            x3 = 5.0;
            y3 = (x3 - 1.0) * (x3 - 1.0) + 4.0;

            result = LocalMinimum.minquad(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            if (!result)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning");
                Console.WriteLine("  MINQUAD returned an error code.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The minimum is to the left of the interval");
                Console.WriteLine("  X1, Y1 = " + x1 + "  " + y1 + "");
                Console.WriteLine("  X2, Y2 = " + x2 + "  " + y2 + "");
                Console.WriteLine("  X3, Y3 = " + x3 + "  " + y3 + "");
                Console.WriteLine("  XMIN = " + xmin + " YMIN = " + ymin + "");
            }

            //
            //  Case 3: the function is flat.
            //
            x1 = 11.0;
            y1 = 6.0;

            x2 = 6.0;
            y2 = 6.0;

            x3 = 2.0;
            y3 = 6.0;

            result = LocalMinimum.minquad(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            if (!result)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning");
                Console.WriteLine("  MINQUAD returned an error code.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The function is flat.");
                Console.WriteLine("  X1, Y1 = " + x1 + "  " + y1 + "");
                Console.WriteLine("  X2, Y2 = " + x2 + "  " + y2 + "");
                Console.WriteLine("  X3, Y3 = " + x3 + "  " + y3 + "");
                Console.WriteLine("  XMIN = " + xmin + " YMIN = " + ymin + "");
            }

            //
            //  Case 4: the function has a maximum.
            //  y = - ( x - 1 )**2 + 4
            //
            x1 = 0.0;
            y1 = -(x1 - 1.0) * (x1 - 1.0) + 4.0;

            x2 = 2.0;
            y2 = -(x2 - 1.0) * (x2 - 1.0) + 4.0;

            x3 = 3.0;
            y3 = -(x3 - 1.0) * (x3 - 1.0) + 4.0;

            result = LocalMinimum.minquad(x1, y1, x2, y2, x3, y3, ref xmin, ref ymin);

            if (!result)
            {
                Console.WriteLine("");
                Console.WriteLine("Warning");
                Console.WriteLine("  MINQUAD returned an error code.");
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  The function has a maximum.");
                Console.WriteLine("  X1, Y1 = " + x1 + "  " + y1 + "");
                Console.WriteLine("  X2, Y2 = " + x2 + "  " + y2 + "");
                Console.WriteLine("  X3, Y3 = " + x3 + "  " + y3 + "");
                Console.WriteLine("  XMIN = " + xmin + " YMIN = " + ymin + "");
            }

        }

    }
}