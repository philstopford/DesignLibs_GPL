﻿using System;
using System.Linq;
using Burkardt.Types;

namespace GeometryTest
{
    public static class HalfTest
    {

        public static void halfplane_contains_point_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HALFPLANE_CONTAINS_POINT_2D_TEST tests HALFPLANE_CONTAINS_POINT_2D
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 4;

            bool[] expected = {true, false, true, false};
            double[] p;
            double[] p_test =
            {
                1.0, 1.0,
                1.0, -1.0,
                -1.0, 1.0,
                2.0, 200.0
            };
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                -5.0, -5.0,
                3.0, 150.0
            };
            double[] p2;
            double[] p2_test =
            {
                2.0, 0.0,
                2.0, 0.0,
                10.0, 10.0,
                1.0, 50.0
            };
            bool temp;
            int test;

            Console.WriteLine("");
            Console.WriteLine("HALFPLANE_CONTAINS_POINT_2D_TEST");
            Console.WriteLine("  HALFPLANE_CONTAINS_POINT_2D determines whether a");
            Console.WriteLine("  halfplane bounded by PA:PB contains the");
            Console.WriteLine("  point P.");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                temp = Burkardt.Geometry.Half.halfplane_contains_point_2d(p1, p2, p);

                Console.WriteLine("");
                Console.WriteLine("  P1 = " + p1[0].ToString().PadLeft(12)
                                            + "  " + p1[1].ToString().PadLeft(12) + "");
                Console.WriteLine("  P2 = " + p2[0].ToString().PadLeft(12)
                                            + "  " + p2[1].ToString().PadLeft(12) + "");
                Console.WriteLine("  P = " + p[0].ToString().PadLeft(12)
                                           + "  " + p[1].ToString().PadLeft(12) + "");
                Console.WriteLine("  Contains? = " + temp
                                                   + "  Correct = " + expected[test] + "");
            }

        }

        public static void test029()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST029 tests HALFSPACE_IMP_TRIANGLE_INT_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            double a;
            double b;
            double c;
            double d;
            int int_num;
            int j;
            int test;
            int test_num = 6;
            double[] p = new double[DIM_NUM * 4];
            double[] t = new double[DIM_NUM * 3];

            Console.WriteLine("");
            Console.WriteLine("TEST029");
            Console.WriteLine("  HALFSPACE_IMP_TRIANGLE_INT_3D finds");
            Console.WriteLine("  intersection points of an implicit");
            Console.WriteLine("  halfspace and a triangle.");

            a = 1.0;
            b = -2.0;
            c = -3.0;
            d = 6.0;

            Console.WriteLine("");
            Console.WriteLine("  The implicitly defined bounding plane");
            Console.WriteLine("  has the form: A*X + B*Y + C*Z + D = 0.");
            Console.WriteLine("  A,B,C,D = " + a + "  " + b + "  " + c + "  " + d + "");

            for (test = 0; test < test_num; test++)
            {

                if (test == 0)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 1)
                {
                    t[0 + 0 * 3] = -6.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 2)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 3.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 2.0;
                }
                else if (test == 3)
                {
                    t[0 + 0 * 3] = -6.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 4.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 3.0;
                }
                else if (test == 4)
                {
                    t[0 + 0 * 3] = -8.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 5)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 4.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 4.0;
                }

                Console.WriteLine("");
                Console.WriteLine("  Case " + test + "");

                typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

                int_num = Burkardt.Geometry.Half.halfspace_imp_triangle_int_3d(a, b, c, d, t, ref p);

                Console.WriteLine("");
                Console.WriteLine("  Number of intersection points is " + int_num + "");
                Console.WriteLine("");

                if (0 < int_num)
                {
                    for (j = 0; j < int_num; j++)
                    {
                        Console.WriteLine("  " + j.ToString().PadLeft(4)
                                               + "  " + p[0 + j * DIM_NUM].ToString().PadLeft(10)
                                               + "  " + p[1 + j * DIM_NUM].ToString().PadLeft(10)
                                               + "  " + p[2 + j * DIM_NUM].ToString().PadLeft(10) + "");
                    }
                }
            }

        }

        public static void test030()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST030 tests HALFSPACE_NORM_TRIANGLE_INT_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;

            int int_num;
            int j;
            double[] p = new double[DIM_NUM * 4];
            double[] pn = {2.0, -4.0, -6.0};
            double[] pp = {-6.0, 0.0, 0.0};
            double[] t = new double[DIM_NUM * 3];
            int test;
            int test_num = 6;

            Console.WriteLine("");
            Console.WriteLine("TEST030");
            Console.WriteLine("  HALFSPACE_NORM_TRIANGLE_INT_3D finds");
            Console.WriteLine("  intersection points of a normal form");
            Console.WriteLine("  halfspace and a triangle.");

            typeMethods.r8vec_print(DIM_NUM, pp, "  A point on the plane:");
            typeMethods.r8vec_print(DIM_NUM, pn, "  The normal vector:");

            for (test = 0; test < test_num; test++)
            {
                if (test == 0)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 1)
                {
                    t[0 + 0 * 3] = -6.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 2)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 3.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 2.0;
                }
                else if (test == 3)
                {
                    t[0 + 0 * 3] = -6.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 4.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 3.0;
                }
                else if (test == 4)
                {
                    t[0 + 0 * 3] = -8.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = -1.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = -2.0;
                }
                else if (test == 5)
                {
                    t[0 + 0 * 3] = 0.0;
                    t[1 + 0 * 3] = 0.0;
                    t[2 + 0 * 3] = 0.0;

                    t[0 + 1 * 3] = 0.0;
                    t[1 + 1 * 3] = 4.0;
                    t[2 + 1 * 3] = 0.0;

                    t[0 + 2 * 3] = 0.0;
                    t[1 + 2 * 3] = 0.0;
                    t[2 + 2 * 3] = 4.0;
                }

                Console.WriteLine("");
                Console.WriteLine("  Case " + test + "");
                typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

                int_num = Burkardt.Geometry.Half.halfspace_norm_triangle_int_3d(pp, pn, t, ref p);

                Console.WriteLine("");
                Console.WriteLine("  Number of intersection points is " + int_num + "");
                Console.WriteLine("");

                if (0 < int_num)
                {
                    for (j = 0; j < int_num; j++)
                    {
                        Console.WriteLine("  " + j.ToString().PadLeft(4)
                                               + "  " + p[0 + j * 3].ToString().PadLeft(10)
                                               + "  " + p[1 + j * 3].ToString().PadLeft(10)
                                               + "  " + p[2 + j * 3].ToString().PadLeft(10) + "");
                    }
                }

            }

        }

    }
}