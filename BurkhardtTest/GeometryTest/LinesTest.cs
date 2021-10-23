using System;
using System.Linq;
using Burkardt;
using Burkardt.LineNS;
using Burkardt.Types;

namespace GeometryTest
{
    public static class LinesTest
    {

        public static void test0327()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0327 tests LINE_EXP_NORMAL_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double[] normal;
            double[] p1 = {1.0, 3.0};
            double[] p2 = {4.0, 0.0};

            Console.WriteLine("");
            Console.WriteLine("TEST0327");
            Console.WriteLine("  LINE_EXP_NORMAL_2D determines a unit normal vector");
            Console.WriteLine("  to a given explicit line.");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point 1: ");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point 2: ");

            normal = Geometry.line_exp_normal_2d(p1, p2);

            typeMethods.r8vec_print(DIM_NUM, normal, "  Normal vector N:");
        }

        public static void line_exp_perp_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_EXP_PERP_2D_TEST tests LINE_EXP_PERP_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 July 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            bool flag = false;
            double[] p1 = {1.0, 3.0};
            double[] p2 = {4.0, 0.0};
            double[] p3;
            double[] p3test =
            {
                0.0, 0.0,
                5.0, -1.0,
                5.0, 3.0
            };
            double[] p4;
            int test;

            Console.WriteLine("");
            Console.WriteLine("LINE_EXP_PERP_2D_TEST");
            Console.WriteLine("  LINE_EXP_PERP_2D is given an explicit line (P1,P2),");
            Console.WriteLine("  and another point P3.  It then finds a point");
            Console.WriteLine("  P4 on (P1,P2) so that (P1,P2) is perpendicular");
            Console.WriteLine("  to (P3,P4).");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point P2:");

            for (test = 0; test < TEST_NUM; test++)
            {
                p3 = p3test.Skip(+test * DIM_NUM).ToArray();

                typeMethods.r8vec_print(DIM_NUM, p3, "  Point P3:");

                p4 = Geometry.line_exp_perp_2d(p1, p2, p3, ref flag);

                typeMethods.r8vec_print(DIM_NUM, p4, "  Point P4:");

            }
        }

        public static void line_exp_point_dist_2d()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_EXP_POINT_DIST_2D_TEST tests LINE_EXP_POINT_DIST_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double dist;
            double[] p;
            double[] p1 = {1.0, 3.0};
            double[] p2 = {4.0, 0.0};
            double[] p_test =
            {
                0.0, 0.0,
                5.0, -1.0,
                5.0, 3.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("LINE_EXP_POINT_DIST_2D_TEST");
            Console.WriteLine("  LINE_EXP_POINT_DIST_2D finds the distance from");
            Console.WriteLine("  an explicit line to a point in 2D.");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point 1: ");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point 2: ");

            for (test = 0; test < TEST_NUM; test++)
            {
                p = p_test.Skip(test * DIM_NUM).ToArray();

                typeMethods.r8vec_print(DIM_NUM, p, "  Point: ");

                dist = Geometry.line_exp_point_dist_2d(p1, p2, p);

                Console.WriteLine("  Distance = " + dist + "");
            }

        }

        public static void test0336()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0336 tests LINE_EXP_POINT_DIST_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int TEST_NUM = 3;

            double dist;
            double[] p;
            double[] p1 = {1.0, 3.0, 2.0};
            double[] p2 = {4.0, 0.0, 1.0};
            double[] p_test =
            {
                0.0, 0.0, 2.0,
                5.0, -1.0, 1.0,
                5.0, 3.0, 3.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST0336");
            Console.WriteLine("  LINE_EXP_POINT_DIST_3D finds the distance");
            Console.WriteLine("  from an explicit line to a point in 3D.");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point 1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point 2:");

            for (test = 0; test < TEST_NUM; test++)
            {
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                typeMethods.r8vec_print(DIM_NUM, p, "  Point:");

                dist = Geometry.line_exp_point_dist_3d(p1, p2, p);

                Console.WriteLine("  Distance = " + dist + "");
            }
        }

        public static void test0337()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0337 tests LINE_EXP_POINT_DIST_SIGNED_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double dist;
            double[] p;
            double[] p1 = {1.0, 3.0};
            double[] p2 = {4.0, 0.0};
            double[] p_test =
            {
                0.0, 0.0,
                5.0, -1.0,
                5.0, 3.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST0337");
            Console.WriteLine("  LINE_EXP_POINT_DIST_SIGNED_2D finds the signed");
            Console.WriteLine("  distance to a point from an explicit line.");

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point 1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point 2:");

            for (test = 0; test < TEST_NUM; test++)
            {
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                typeMethods.r8vec_print(DIM_NUM, p, "  Point:");

                dist = Geometry.line_exp_point_dist_signed_2d(p1, p2, p);

                Console.WriteLine("  Signed distance = " + dist + "");
            }

        }

        public static void test034()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST034 tests LINE_EXP_POINT_NEAR_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double dist = 0;
            double[] p;
            double[] p1 = {1.0, 3.0};
            double[] p2 = {4.0, 0.0};
            double[] pn = new double[DIM_NUM];
            double[] p_test =
            {
                0.0, 0.0,
                5.0, -1.0,
                5.0, 3.0
            };
            double t = 0;
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST034");
            Console.WriteLine("  LINE_EXP_POINT_NEAR_2D finds the point on");
            Console.WriteLine("  a line nearest in point in 2D.");

            typeMethods.r8vec_print(DIM_NUM, p1, "  The point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  The point P2:");

            for (test = 0; test < TEST_NUM; test++)
            {
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                typeMethods.r8vec_print(DIM_NUM, p, "  The point P:");

                Geometry.line_exp_point_near_2d(p1, p2, p, ref pn, ref dist, ref t);

                typeMethods.r8vec_print(DIM_NUM, pn, "  Nearest point PN:");

                Console.WriteLine("  Distance = " + dist + "");
                Console.WriteLine("  Relative line position T = " + t + "");
            }
        }

        public static void test0345()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0345 tests LINE_EXP2IMP_2D and LINE_IMP2EXP_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double a1 = 1.0;
            double a2 = 0;
            double b1 = 2.0;
            double b2 = 0;
            double c1 = 3.0;
            double c2 = 0;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TEST0345");
            Console.WriteLine("  LINE_EXP2IMP_2D converts explicit to implicit lines.");
            Console.WriteLine("  LINE_IMP2EXP_2D converts implicit to explicit lines.");

            Console.WriteLine("");
            Console.WriteLine("  Implicit line A = " + a1
                                                     + " B = " + b1
                                                     + " C = " + c1 + "");

            Geometry.line_imp2exp_2d(a1, b1, c1, ref p1, ref p2);

            typeMethods.r8vec_print(DIM_NUM, p1, "  The point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  The point P2:");

            Geometry.line_exp2imp_2d(p1, p2, ref a2, ref b2, ref c2);

            Console.WriteLine("");
            Console.WriteLine("  Recovered implicit line A = " + a2
                                                               + " B = " + b2
                                                               + " C = " + c2 + "");

        }

        public static void test0346()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0346 tests LINE_EXP2PAR_2D and LINE_PAR2EXP_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double f1 = 1.0;
            double f2 = 0;
            double g1 = 2.0;
            double g2 = 0;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double x1 = 3.0;
            double x2 = 0;
            double y1 = 4.0;
            double y2 = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST0346");
            Console.WriteLine("  LINE_EXP2PAR_2D converts explicit to parametric lines.");
            Console.WriteLine("  LINE_PAR2EXP_2D converts parametric to explicit lines.");

            Console.WriteLine("");
            Console.WriteLine("  Parametric line:");
            Console.WriteLine("    F  = " + f1 + " G  = " + g1 + "");
            Console.WriteLine("    X0 = " + x1 + " Y0 = " + y1 + "");

            Geometry.line_par2exp_2d(f1, g1, x1, y1, ref p1, ref p2);

            typeMethods.r8vec_print(DIM_NUM, p1, "  The point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  The point P2:");

            Geometry.line_exp2par_2d(p1, p2, ref f2, ref g2, ref x2, ref y2);

            Console.WriteLine("");
            Console.WriteLine("  Recovered parametric line:");
            Console.WriteLine("    F  = " + f2 + " G  = " + g2 + "");
            Console.WriteLine("    X0 = " + x2 + " Y0 = " + y2 + "");

        }

        public static void line_imp_point_dist_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_IMP_POINT_DIST_2D_TEST tests LINE_IMP_POINT_DIST_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double a = 0;
            double[] atest = {2.0, 2.0, 2.0};
            double b = 0;
            double[] btest = {5.0, 5.0, 5.0};
            double c = 0;
            double[] ctest = {3.0, 3.0, 3.0};
            double dist = 0;
            double[] p;
            double[] p_test =
            {
                0.0, 6.0,
                0.0, 5.0,
                0.0, 4.0
            };
            int test = 0;

            Console.WriteLine("");
            Console.WriteLine("LINE_IMP_POINT_DIST_2D_TEST");
            Console.WriteLine("  LINE_IMP_POINT_DIST_2D finds the distance from");
            Console.WriteLine("  a point P to a line A * X + B * Y + C = 0.");
            Console.WriteLine("");
            Console.WriteLine("   X       Y       A       B       C       DIST");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                a = atest[test];
                b = btest[test];
                c = ctest[test];
                p = p_test.Skip(+test * DIM_NUM).ToArray();

                dist = Geometry.line_imp_point_dist_2d(a, b, c, p);

                Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                       + "  " + p[1].ToString().PadLeft(8)
                                       + "  " + a.ToString().PadLeft(8)
                                       + "  " + b.ToString().PadLeft(8)
                                       + "  " + c.ToString().PadLeft(8)
                                       + "  " + dist.ToString().PadLeft(8) + "");
            }
        }

        public static void test0351()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0351 tests LINE_PAR_POINT_NEAR_2D and LINE_PAR_POINT_DIST_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double dist;
            double f;
            double g;
            int i;
            double[] p = new double [2];
            double[] pn;
            double[] p_test =
            {
                0.0, 0.0,
                5.0, -1.0,
                5.0, 3.0
            };
            int test;
            int test_num = TEST_NUM;
            double x0;
            double y0;

            Console.WriteLine("");
            Console.WriteLine("TEST0351");
            Console.WriteLine("  LINE_PAR_POINT_NEAR_2D finds the point on");
            Console.WriteLine("  a parametric line (X0,Y0,F,G) nearest a point P in 2D.");

            x0 = 1.0;
            y0 = 3.0;
            f = +1.0;
            g = -1.0;

            Console.WriteLine("");
            Console.WriteLine("  Parametric line:");
            Console.WriteLine("  X(t) = " + x0 + " + " + f + " * t");
            Console.WriteLine("  Y(t) = " + y0 + " + " + g + " * t");

            for (test = 0; test < test_num; test++)
            {
                for (i = 0; i < 2; i++)
                {
                    p[i] = p_test[i + test * 2];
                }

                typeMethods.r8vec_print(2, p, "  The point P:");

                dist = Geometry.line_par_point_dist_2d(f, g, x0, y0, p);

                Console.WriteLine("  Distance = " + dist + "");

                pn = Geometry.line_par_point_near_2d(f, g, x0, y0, p);

                typeMethods.r8vec_print(2, pn, "  Nearest point PN:");

                dist = typeMethods.r8vec_norm_affine(2, p, pn);

                Console.WriteLine("  Distance recomputed = " + dist + "");


            }

            return;

        }

        public static void test0352()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0352 tests LINE_PAR_POINT_DIST_3D and LINE_PAR_POINT_NEAR_3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double dist;
            double f;
            double g;
            double h;
            int i;
            double[] p = new double[3];
            double[] pn;
            double[] p_test =
            {
                0.0, 0.0, 2.0,
                5.0, -1.0, 1.0,
                5.0, 3.0, 3.0
            };
            int test;
            int test_num = TEST_NUM;
            double x0;
            double y0;
            double z0;

            Console.WriteLine("");
            Console.WriteLine("TEST0352");
            Console.WriteLine("  LINE_PAR_POINT_DIST_3D finds the distance");
            Console.WriteLine("  from a parametric line to a point in 3D.");

            x0 = 1.0;
            y0 = 3.0;
            z0 = 2.0;

            f = +3.0;
            g = -3.0;
            h = -1.0;

            Console.WriteLine("");
            Console.WriteLine("  Parametric line:");
            Console.WriteLine("  X(t) = " + x0 + " + " + f + " * t");
            Console.WriteLine("  Y(t) = " + y0 + " + " + g + " * t");
            Console.WriteLine("  Z(t) = " + z0 + " + " + h + " * t");

            for (test = 0; test < test_num; test++)
            {
                for (i = 0; i < 3; i++)
                {
                    p[i] = p_test[i + test * 3];
                }

                typeMethods.r8vec_print(3, p, "  The point P:");

                dist = Geometry.line_par_point_dist_3d(f, g, h, x0, y0, z0, p);

                Console.WriteLine("  Distance = " + dist + "");

                pn = Geometry.line_par_point_near_3d(f, g, h, x0, y0, z0, p);

                typeMethods.r8vec_print(3, pn, "  Nearest point PN:");

                dist = typeMethods.r8vec_norm_affine(3, p, pn);

                Console.WriteLine("  Distance recomputed = " + dist + "");


            }

            return;

        }

        public static void test038()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST038 tests LINES_EXP_ANGLE_3D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int TEST_NUM = 2;

            double angle;
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0, 0.0,
                1.0, 2.0, 0.0
            };
            double[] p2;
            double[] p2_test =
            {
                1.0, 2.0, 0.0,
                1.0, 2.0, 0.0
            };
            double[] q1;
            double[] q1_test =
            {
                0.0, 3.0, 3.0,
                1.0, 2.0, -1.0
            };
            double[] q2;
            double[] q2_test =
            {
                3.0, 0.0, 3.0,
                1.0, 2.0, 3.0
            };
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST038");
            Console.WriteLine("  LINES_EXP_ANGLE_3D finds the angle between");
            Console.WriteLine("  two explicit lines in 3D;");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
                q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

                angle = Geometry.lines_exp_angle_3d(p1, p2, q1, q2);

                Console.WriteLine("");
                Console.WriteLine("  Angle between lines is " + angle + "");
            }

        }

        public static void test0385()

            //****************************************************************************80*
            //
            //  Purpose:
            //
            //    TEST0385 tests LINES_EXP_DIST_3D and LINES_EXP_DIST_3D_2;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int TEST_NUM = 2;

            double dist;
            double dist2;
            int i;
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0, 0.0,
                4.0, -3.0, 0.0
            };
            double[] p2;
            double[] p2_test =
            {
                1.0, 2.0, 0.0,
                -8.0, 6.0, 0.0
            };
            double[] q1;
            double[] q1_test =
            {
                0.0, 3.0, 3.0,
                3.0, 4.0, -1.0
            };
            double[] q2;
            double[] q2_test =
            {
                3.0, 0.0, 3.0,
                3.0, 4.0, 3.0
            };
            int test;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("TEST0385");
            Console.WriteLine("  LINES_EXP_DIST_3D finds the distance between");
            Console.WriteLine("  two explicit lines in 3D;");
            Console.WriteLine("  LINES_EXP_DIST_3D_2 finds the distance between");
            Console.WriteLine("  two explicit lines in 3D;");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
                q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

                dist = Geometry.lines_exp_dist_3d(p1, p2, q1, q2);
                dist2 = Geometry.lines_exp_dist_3d_2(p1, p2, q1, q2);

                Console.WriteLine("");
                Console.WriteLine("");
                cout = "  P1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  P2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                Console.WriteLine("  LINES_EXP_DIST_3D =   " + dist + "");
                Console.WriteLine("  LINES_EXP_DIST_3D_2 = " + dist2 + "");
            }

        }

        public static void test03855()

            //****************************************************************************80*
            //
            //  Purpose:
            //
            //    TEST03855 tests LINES_EXP_DIST_3D and LINES_EXP_DIST_3D_2;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int TEST_NUM = 2;

            int i;
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0, 0.0,
                4.0, -3.0, 0.0
            };
            double[] p2;
            double[] p2_test =
            {
                1.0, 2.0, 0.0,
                -8.0, 6.0, 0.0
            };
            double[] pn = new double[DIM_NUM];
            double[] q1;
            double[] q1_test =
            {
                0.0, 3.0, 3.0,
                3.0, 4.0, -1.0
            };
            double[] q2;
            double[] q2_test =
            {
                3.0, 0.0, 3.0,
                3.0, 4.0, 3.0
            };
            double[] qn = new double[DIM_NUM];
            int test;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("TEST03855");
            Console.WriteLine("  LINES_EXP_NEAR_3D finds nearest points on");
            Console.WriteLine("  two explicit lines in 3D;");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
                q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

                Console.WriteLine("");
                Console.WriteLine("");
                cout = "  P1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  P2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);

                Geometry.lines_exp_near_3d(p1, p2, q1, q2, ref pn, ref qn);

                Console.WriteLine("");
                cout = "  PN:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + pn[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  QN:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + qn[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
            }

        }

        public static void test0386()

            //****************************************************************************80*
            //
            //  Purpose:
            //
            //    TEST0386 tests LINES_EXP_EQUAL_2D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 6;

            bool equal;
            int i;
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0
            };
            double[] p2;
            double[] p2_test =
            {
                1.0, 2.0,
                1.0, 2.0,
                1.0, 2.0,
                1.0, 2.0,
                1.0, 2.0,
                1.0, 2.0
            };
            double[] q1;
            double[] q1_test =
            {
                0.0, 0.0,
                1.0, 2.0,
                0.0, 0.0,
                7.0, 14.0,
                1.0, 2.0,
                0.0, 10.0
            };
            double[] q2;
            double[] q2_test =
            {
                1.0, 2.0,
                0.0, 0.0,
                2.0, 4.0,
                5.5, 11.0,
                3.0, 5.0,
                1.0, 12.0
            };
            int test;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("TEST0386");
            Console.WriteLine("  LINES_EXP_EQUAL_2D tries to determine if two");
            Console.WriteLine("  explicit lines in 2D are equal;");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
                q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

                Console.WriteLine("");
                cout = "  P1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  P2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);

                equal = Geometry.lines_exp_equal_2d(p1, p2, q1, q2);

                if (equal)
                {
                    Console.WriteLine("  The lines are equal.");
                }
                else
                {
                    Console.WriteLine("  The lines are distinct.");
                }
            }

        }

        public static void lines_exp_int_2d_test()

            //****************************************************************************80*
            //
            //  Purpose:
            //
            //    LINES_EXP_INT_2D_TEST tests LINES_EXP_INT_2D;
            //
            //  Discussion:
            //
            //    Test #1:
            //
            //      x + 2y -  4 = 0
            //      x -  y -  1 = 0
            //
            //    Test #2:
            //
            //      x + 2y -  4 = 0
            //     2x + 4y -  1 = 0
            //
            //    Test #3:
            //
            //      x + 2y -  4 = 0
            //    -3x - 6y + 12 = 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            int i;
            int ival = 0;
            double[] p = new double[DIM_NUM];
            double[] p1;
            double[] p1_test =
            {
                0.0, 0.0,
                0.0, 2.0,
                0.0, 2.0
            };
            double[] p2;
            double[] p2_test =
            {
                4.0, 0.0,
                4.0, 0.0,
                4.0, 0.0
            };
            double[] q1;
            double[] q1_test =
            {
                0.0, -1.0,
                0.0, 0.25,
                0.0, 2.0
            };
            double[] q2;
            double[] q2_test =
            {
                1.0, 0.0,
                0.5, 0.0,
                4.0, 0.0
            };
            int test;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine("LINES_EXP_INT_2D_TEST");
            Console.WriteLine("  LINES_EXP_INT_2D finds intersections of");
            Console.WriteLine("  two explicit lines in 2D;");

            for (test = 0; test < TEST_NUM; test++)
            {
                p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
                p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
                q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
                q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

                Console.WriteLine("");
                cout = "  P1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  P2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q1:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q1[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);
                cout = "  Q2:";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + q2[i].ToString().PadLeft(12);
                }

                Console.WriteLine(cout);

                Geometry.lines_exp_int_2d(p1, p2, q1, q2, ref ival, ref p);

                if (ival == 1)
                {
                    cout = "  Intersection at";
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + p[i].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                }
                else if (ival == 0)
                {
                    Console.WriteLine("  Lines are parallel, no intersection.");
                }
                else if (ival == 2)
                {
                    Console.WriteLine("  Lines are coincident.");
                }
                else
                {
                    Console.WriteLine("  Unknown return value of IVAL = " + ival + "");
                }
            }

        }

        public static void test040()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST040 tests LINES_IMP_ANGLE_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double a1;
            double a2;
            double angle;
            double b1;
            double b2;
            double c1;
            double c2;

            Console.WriteLine("");
            Console.WriteLine("TEST040");
            Console.WriteLine("  For two lines written in implicit form:");
            Console.WriteLine("  LINES_IMP_ANGLE_2D finds the angle;");
            Console.WriteLine("");
            //
            //  x + 2y - 4 = 0
            //
            a1 = 1.0;
            b1 = 2.0;
            c1 = -4.0;

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  x - y - 1 = 0
            //
            a2 = 1.0;
            b2 = -1.0;
            c2 = -1.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            angle = Geometry.lines_imp_angle_2d(a1, b1, c1, a2, b2, c2);

            Console.WriteLine("");
            Console.WriteLine("  Angle between lines is " + Helpers.radians_to_degrees(angle) + "");

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  2x + 4y - 1 = 0
            //
            a2 = 2.0;
            b2 = +4.0;
            c2 = -1.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            angle = Geometry.lines_imp_angle_2d(a1, b1, c1, a2, b2, c2);

            Console.WriteLine("");
            Console.WriteLine("  Angle between lines is " + Helpers.radians_to_degrees(angle) + "");

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  -3x - 6y +12 = 0
            //
            a2 = -3.0;
            b2 = -6.0;
            c2 = +12.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            angle = Geometry.lines_imp_angle_2d(a1, b1, c1, a2, b2, c2);

            Console.WriteLine("");
            Console.WriteLine("  Angle between lines is " + Helpers.radians_to_degrees(angle) + "");

        }

        public static void test041()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST041 tests LINES_IMP_DIST_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 3;

            double a1;
            double[] a1_test = {4.0, 2.0, 1.0};
            double a2;
            double[] a2_test = {4.0, 4.0, 2.0};
            double b1;
            double[] b1_test = {-1.0, -1.0, 2.0};
            double b2;
            double[] b2_test = {-1.0, -2.0, 3.0};
            double c1;
            double[] c1_test = {3.0, 0.0, 2.0};
            double c2;
            double[] c2_test = {12.0, 6.0, 1.0};
            double dist;
            int test;

            Console.WriteLine("");
            Console.WriteLine("TEST041");
            Console.WriteLine("  LINES_IMP_DIST_3D finds the distance between");
            Console.WriteLine("  two implicit lines in 2D.");
            Console.WriteLine("");
            Console.WriteLine("   A1      B1      C1      A2      B2      C2   DIST");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                a1 = a1_test[test];
                b1 = b1_test[test];
                c1 = c1_test[test];
                a2 = a2_test[test];
                b2 = b2_test[test];
                c2 = c2_test[test];

                dist = Geometry.lines_imp_dist_2d(a1, b1, c1, a2, b2, c2);

                Console.WriteLine("  " + a1.ToString().PadLeft(8)
                                       + "  " + b1.ToString().PadLeft(8)
                                       + "  " + c1.ToString().PadLeft(8)
                                       + "  " + a2.ToString().PadLeft(8)
                                       + "  " + b2.ToString().PadLeft(8)
                                       + "  " + c2.ToString().PadLeft(8)
                                       + "  " + dist.ToString().PadLeft(8) + "");
            }

        }

        public static void test0415()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0415 tests LINES_IMP_INT_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double a1;
            double a2;
            double b1;
            double b2;
            double c1;
            double c2;
            int ival = 0;
            double[] p = new double[2];

            Console.WriteLine("");
            Console.WriteLine("TEST0415");
            Console.WriteLine("  For two lines written in implicit form:");
            Console.WriteLine("  LINES_IMP_INT_2D finds the intersection.");
            Console.WriteLine("");
            //
            //  x + 2y - 4 = 0
            //
            a1 = 1.0;
            b1 = 2.0;
            c1 = -4.0;

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  x - y - 1 = 0
            //
            a2 = 1.0;
            b2 = -1.0;
            c2 = -1.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            Geometry.lines_imp_int_2d(a1, b1, c1, a2, b2, c2, ref ival, ref p);

            if (ival == 1)
            {
                Console.WriteLine("  Intersection at " + p[0] + "  " + p[1] + "");
            }
            else if (ival == 0)
            {
                Console.WriteLine("  Lines are parallel, no intersection.");
            }
            else if (ival == 2)
            {
                Console.WriteLine("  Lines are coincident.");
            }
            else
            {
                Console.WriteLine("  Unknown return value of ival = " + ival + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  2x + 4y - 1 = 0
            //
            a2 = 2.0;
            b2 = +4.0;
            c2 = -1.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            Geometry.lines_imp_int_2d(a1, b1, c1, a2, b2, c2, ref ival, ref p);

            if (ival == 1)
            {
                Console.WriteLine("  Intersection at " + p[0] + "  " + p[1] + "");
            }
            else if (ival == 0)
            {
                Console.WriteLine("  Lines are parallel, no intersection.");
            }
            else if (ival == 2)
            {
                Console.WriteLine("  Lines are coincident.");
            }
            else
            {
                Console.WriteLine("  Unknown return value of IVAL = " + ival + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Line 1 coefficients: " + a1 + "  " + b1 + "  " + c1 + "");
            //
            //  -3x - 6y +12 = 0
            //
            a2 = -3.0;
            b2 = -6.0;
            c2 = +12.0;
            Console.WriteLine("  Line 2 coefficients: " + a2 + "  " + b2 + "  " + c2 + "");

            Geometry.lines_imp_int_2d(a1, b1, c1, a2, b2, c2, ref ival, ref p);

            if (ival == 1)
            {
                Console.WriteLine("  Intersection at " + p[0] + "  " + p[1] + "");
            }
            else if (ival == 0)
            {
                Console.WriteLine("  Lines are parallel, no intersection.");
            }
            else if (ival == 2)
            {
                Console.WriteLine("  Lines are coincident.");
            }
            else
            {
                Console.WriteLine("  Unknown return value of IVAL = " + ival + "");
            }

        }

        public static void test0416()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0416 tests LINES_PAR_INT_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double f1 = 0;
            double f2 = 0;
            double g1 = 0;
            double g2 = 0;
            double[] pint = new double[DIM_NUM];
            double t1 = 0;
            double t2 = 0;
            double x1 = 0;
            double x2 = 0;
            double y1 = 0;
            double y2 = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST0416");
            Console.WriteLine("  LINES_PAR_INT_2D finds the intersection of");
            Console.WriteLine("  two lines written in parametric form.");
            Console.WriteLine("");
            //
            //  x - 2y = -1
            //
            x1 = 0.0;
            y1 = 1.0;
            f1 = 2.0;
            g1 = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Line 1 parameters:" + x1.ToString().PadLeft(8)
                                                     + "  " + y1.ToString().PadLeft(8)
                                                     + "  " + f1.ToString().PadLeft(8)
                                                     + "  " + g1.ToString().PadLeft(8) + "");
            //
            //  x + y - 8 = 0
            //
            x2 = 10.0;
            y2 = -2.0;
            f2 = 1.0;
            g2 = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Line 2 parameters:" + x2.ToString().PadLeft(8)
                                                     + "  " + y2.ToString().PadLeft(8)
                                                     + "  " + f2.ToString().PadLeft(8)
                                                     + "  " + g2.ToString().PadLeft(8) + "");

            Geometry.lines_par_int_2d(f1, g1, x1, y1, f2, g2, x2, y2, ref t1, ref t2, ref pint);

            Console.WriteLine("");
            Console.WriteLine("  Line 1 evaluated at T1:");
            Console.WriteLine("");
            Console.WriteLine("    T1 =   " + t1 + "");
            Console.WriteLine("    X(T1)= " + x1 + f1 * t1 + "");
            Console.WriteLine("    Y(T1)= " + y1 + g1 * t1 + "");
            Console.WriteLine("");
            Console.WriteLine("  Line 2 evaluated at T2:");
            Console.WriteLine("");
            Console.WriteLine("    T2 =   " + t2 + "");
            Console.WriteLine("    X(T2)= " + x2 + f2 * t2 + "");
            Console.WriteLine("    Y(T2)= " + y2 + g2 * t2 + "");
            Console.WriteLine("");
            Console.WriteLine("  Reported intersection PINT:");
            Console.WriteLine("");
            Console.WriteLine("    " + pint[0] + "");
            Console.WriteLine("    " + pint[1] + "");

        }

    }
}