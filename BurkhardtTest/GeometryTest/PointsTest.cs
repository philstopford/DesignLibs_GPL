using System;
using System.Linq;
using Burkardt.PointsNS;
using Burkardt.Types;

namespace GeometryTest
{
    public static class PointsTest
    {
        public static void points_centroid_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_CENTROID_2D_TEST tests POINTS_CENTROID_2D.
            //
            //  !....3&11....
            //  !............
            //  !............
            //  X..9.........
            //  !.....5......
            //  !...........6
            //  !.4.2...10...
            //  !.....8...12.
            //  V............
            //  !..7.........
            //  !......1.....
            //  !............
            //  !............
            //  !----V----X--
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int N = 12;

            int centroid_index;
            double[] p =
            {
                7.0, 3.0,
                4.0, 7.0,
                5.0, 13.0,
                2.0, 7.0,
                6.0, 9.0,
                12.0, 8.0,
                3.0, 4.0,
                6.0, 6.0,
                3.0, 10.0,
                8.0, 7.0,
                5.0, 13.0,
                10.0, 6.0
            };

            Console.WriteLine("");
            Console.WriteLine("POINTS_CENTROID_2D_TEST");
            Console.WriteLine("  POINTS_CENTROID_2D computes the centroid of a");
            Console.WriteLine("  discrete set of points.");

            typeMethods.r8mat_transpose_print(DIM_NUM, N, p, "  The points:");

            centroid_index = Geometry.points_centroid_2d(N, p);

            Console.WriteLine("");
            Console.WriteLine("  The centroid is point #: " + centroid_index + "");
        }

        public static void points_colin_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POINTS_COLIN_2D_TEST tests POINTS_COLIN_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double colin;
            int i;
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0
            };
            double[] p2;
            double[] p2_test =
            {
                10.0, 10.0,
                0.0, 1.0,
                1.0, 0.0
            };
            double[] p3;
            double[] p3_test =
            {
                5.0, 4.99,
                100.0, 0.0,
                0.5, 0.86602539
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("POINTS_COLIN_2D_TEST");
            Console.WriteLine("  POINTS_COLIN_2D estimates the colinearity");
            Console.WriteLine("  of three points.");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+DIM_NUM * test).ToArray();
                p2 = p2_test.Skip(+DIM_NUM * test).ToArray();
                p3 = p3_test.Skip(+DIM_NUM * test).ToArray();

                Console.WriteLine("");

                if (test == 1)
                {
                    Console.WriteLine("  Points almost on a line: Expect tiny COLIN.");
                }
                else if (test == 2)
                {
                    Console.WriteLine("  Two points close, one far: Expect tiny COLIN.");
                }
                else if (test == 3)
                {
                    Console.WriteLine("  Equilateral triangle: Expect COLIN = 1.");
                }

                colin = Geometry.points_colin_2d(p1, p2, p3);

                Console.WriteLine("");
                string cout = "  P1: ";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                cout = "  P2: ";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p2[i].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);
                cout = "  P3: ";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p3[i].ToString().PadLeft(8);
                }

                Console.WriteLine(cout);

                Console.WriteLine("");
                Console.WriteLine("  Colinearity index = " + colin + "");
            }

        }

    }
}