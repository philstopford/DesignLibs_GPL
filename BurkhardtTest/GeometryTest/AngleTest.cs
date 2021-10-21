using System;
using Burkardt;
using Burkardt.Geometry;
using Burkardt.Types;

namespace GeometryTest
{
    public static class AngleTest
    {
        public static void angle_box_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_BOX_2D_TEST tests ANGLE_BOX_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double dist;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];
            double[] p4 = new double[DIM_NUM];
            double[] p5 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("ANGLE_BOX_2D_TEST");
            Console.WriteLine("  ANGLE_BOX_2D");
            Console.WriteLine("");
            Console.WriteLine("  Compute P4 and P5, normal to");
            Console.WriteLine("  line through P1 and P2, and");
            Console.WriteLine("  line through P2 and P3,");
            Console.WriteLine("  and DIST units from P2.");
            //
            //  These points define the lines
            //    y = 0
            //  and
            //    y = 2x-6
            //
            p1[0] = 0.0;
            p1[1] = 0.0;
            p2[0] = 3.0;
            p2[1] = 0.0;
            p3[0] = 4.0;
            p3[1] = 2.0;
            dist = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  DIST " + dist.ToString().PadLeft(14) + "");
            Console.WriteLine("  P1:  " + p1[0].ToString().PadLeft(14) + "  " + p1[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P2:  " + p2[0].ToString().PadLeft(14) + "  " + p2[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P3:  " + p3[0].ToString().PadLeft(14) + "  " + p3[1].ToString().PadLeft(14) + "");

            Angle.angle_box_2d(dist, p1, p2, p3, ref p4, ref p5);

            Console.WriteLine("  P4:  " + p4[0].ToString().PadLeft(14) + "  " + p4[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P5:  " + p5[0].ToString().PadLeft(14) + "  " + p5[1].ToString().PadLeft(14) + "");
            //
            //  These points define the lines
            //    y = 0
            //  and
            //    y = 2x-6
            //
            p1[0] = 0.0;
            p1[1] = 0.0;
            p2[0] = 3.0;
            p2[1] = 0.0;
            p3[0] = 2.0;
            p3[1] = -2.0;
            dist = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  DIST " + dist.ToString().PadLeft(14) + "");
            Console.WriteLine("  P1:  " + p1[0].ToString().PadLeft(14) + "  " + p1[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P2:  " + p2[0].ToString().PadLeft(14) + "  " + p2[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P3:  " + p3[0].ToString().PadLeft(14) + "  " + p3[1].ToString().PadLeft(14) + "");

            Angle.angle_box_2d(dist, p1, p2, p3, ref p4, ref p5);

            Console.WriteLine("  P4:  " + p4[0].ToString().PadLeft(14) + "  " + p4[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P5:  " + p5[0].ToString().PadLeft(14) + "  " + p5[1].ToString().PadLeft(14) + "");
            //
            //  By setting P1 = P2, we are asking to extend the line
            //    y = 2x-6
            //  from P3 to P2 through to the other side.
            //
            p1[0] = 3.0;
            p1[1] = 0.0;
            p2[0] = 3.0;
            p2[1] = 0.0;
            p3[0] = 2.0;
            p3[1] = -2.0;
            dist = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  DIST " + dist.ToString().PadLeft(14) + "");
            Console.WriteLine("  P1:  " + p1[0].ToString().PadLeft(14) + "  " + p1[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P2:  " + p2[0].ToString().PadLeft(14) + "  " + p2[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P3:  " + p3[0].ToString().PadLeft(14) + "  " + p3[1].ToString().PadLeft(14) + "");

            Angle.angle_box_2d(dist, p1, p2, p3, ref p4, ref p5);

            Console.WriteLine("  P4:  " + p4[0].ToString().PadLeft(14) + "  " + p4[1].ToString().PadLeft(14) + "");
            Console.WriteLine("  P5:  " + p5[0].ToString().PadLeft(14) + "  " + p5[1].ToString().PadLeft(14) + "");

        }

        public static void angle_contains_ray_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_CONTAINS_RAY_2D_TEST tests ANGLE_CONTAINS_RAY_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 6;

            int angle;
            int angle_num = 12;
            double angle_rad;
            bool inside;
            double[] p = new double[2];
            double[] p1 = new double[2];
            double[] p2 = new double[2];
            double[] p3 = new double[2];
            int test;

            Console.WriteLine("");
            Console.WriteLine("ANGLE_CONTAINS_RAY_2D_TEST");
            Console.WriteLine("  ANGLE_CONTAINS_RAY_2D sees if a ray lies within an angle.");

            for (test = 0; test < TEST_NUM; test++)
            {
                if (test == 0)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 1.0;
                    p3[1] = 1.0;
                }
                else if (test == 1)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = 1.0;
                }
                else if (test == 2)
                {
                    p1[0] = 1.0;
                    p1[1] = -1.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = 1.0;
                }
                else if (test == 3)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = -1.0;
                    p3[1] = 0.0;
                }
                else if (test == 4)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = -1.0;
                }
                else if (test == 5)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 1.0;
                    p3[1] = -0.01;
                }

                typeMethods.r8vec_print(2, p1, "  Vertex A");
                typeMethods.r8vec_print(2, p2, "  Vertex B");
                typeMethods.r8vec_print(2, p3, "  Vertex C");

                Console.WriteLine("");
                Console.WriteLine("       X            Y       Inside?");
                Console.WriteLine("");

                for (angle = 0; angle <= angle_num; angle++)
                {
                    angle_rad = (double) (angle) * 2.0 * Math.PI / (double) angle_num;

                    p[0] = Math.Cos(angle_rad);
                    p[1] = Math.Sin(angle_rad);

                    inside = Angle.angle_contains_ray_2d(p1, p2, p3, p);

                    Console.WriteLine("  " + p[0].ToString().PadLeft(12)
                                           + "  " + p[1].ToString().PadLeft(12)
                                           + "  " + inside + "");
                }

            }
        }

        public static void angle_deg_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_DEG_2D_TEST tests ANGLE_DEG_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            int i;
            int angle_num = 12;
            double temp1;
            double temp2;
            double thetad;
            double thetar;
            double[] v1 = new double[DIM_NUM];
            double[] v2 = new double[DIM_NUM];
            double[] v3 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("TANGLE_DEG_2D_TEST");
            Console.WriteLine("  ANGLE_DEG_2D computes an angle;");
            Console.WriteLine("");
            Console.WriteLine("      X           Y          Theta         atan2  ANGLE_DEG_2D");
            Console.WriteLine("");

            v1[0] = 1.0;
            v1[1] = 0.0;
            v3[0] = 0.0;
            v3[1] = 0.0;

            for (i = 0; i <= angle_num; i++)
            {
                thetad = (double) (i) * 360.0 / (double) (angle_num);
                thetar = Helpers.degrees_to_radians(thetad);

                v2[0] = Math.Cos(thetar);
                v2[1] = Math.Sin(thetar);

                temp1 = Helpers.radians_to_degrees(Math.Atan2(v2[1], v2[0]));

                temp2 = Angle.angle_deg_2d(v1, v3, v2);

                Console.WriteLine("  " + v2[0].ToString().PadLeft(10)
                                       + "  " + v2[1].ToString().PadLeft(10)
                                       + "  " + thetad.ToString().PadLeft(10)
                                       + "  " + temp1.ToString().PadLeft(10)
                                       + "  " + temp2.ToString().PadLeft(10) + "");
            }

        }

        public static void angle_half_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_HALF_2D_TEST tests ANGLE_HALF_2D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            double angle_deg;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];
            double[] p4;
            double r;

            Console.WriteLine("");
            Console.WriteLine("ANGLE_HALF_2D_TEST");
            Console.WriteLine("  ANGLE_HALF_2D computes the half angle between two rays;");
            Console.WriteLine("  The angle is defined by the points (P1,P2,P3)");
            Console.WriteLine("  or by the rays P2-->P3, P2-->P1");

            p2[0] = 5.0;
            p2[1] = 3.0;

            angle_deg = 75.0;
            r = 3.0;
            p1[0] = p2[0] + r * typeMethods.r8_cosd(angle_deg);
            p1[1] = p2[1] + r * typeMethods.r8_sind(angle_deg);

            angle_deg = 15.0;
            r = 2.0;
            p3[0] = p2[0] + r * typeMethods.r8_cosd(angle_deg);
            p3[1] = p2[1] + r * typeMethods.r8_sind(angle_deg);

            typeMethods.r8vec_print(DIM_NUM, p1, "  Point P1:");
            typeMethods.r8vec_print(DIM_NUM, p2, "  Point P2:");
            typeMethods.r8vec_print(DIM_NUM, p3, "  Point P3:");

            p4 = Angle.angle_half_2d(p1, p2, p3);

            typeMethods.r8vec_print(DIM_NUM, p4,
                "  End point of unit ray from P2, defining half angle, P4:");

            angle_deg = 45.0;
            r = 1.0;
            p4[0] = p2[0] + r * typeMethods.r8_cosd(angle_deg);
            p4[1] = p2[1] + r * typeMethods.r8_sind(angle_deg);

            typeMethods.r8vec_print(DIM_NUM, p4, "  Expected value of P4:");


        }

        public static void angle_rad_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_2D_TEST tests ANGLE_RAD_2D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 6;

            double angle;
            double[] p1 = new double[DIM_NUM];
            double[] p2 = new double[DIM_NUM];
            double[] p3 = new double[DIM_NUM];
            int test;

            Console.WriteLine("");
            Console.WriteLine("ANGLE_RAD_2D_TEST");
            Console.WriteLine("  ANGLE_RAD_2D computes the angle between two rays;");
            Console.WriteLine("");

            for (test = 1; test <= TEST_NUM; test++)
            {
                if (test == 1)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 1.0;
                    p3[1] = 1.0;
                }
                else if (test == 2)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = 1.0;
                }
                else if (test == 3)
                {
                    p1[0] = 1.0;
                    p1[1] = -1.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = 1.0;
                }
                else if (test == 4)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = -1.0;
                    p3[1] = 0.0;
                }
                else if (test == 5)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 0.0;
                    p3[1] = -1.0;
                }
                else if (test == 6)
                {
                    p1[0] = 1.0;
                    p1[1] = 0.0;

                    p2[0] = 0.0;
                    p2[1] = 0.0;

                    p3[0] = 1.0;
                    p3[1] = -0.01;
                }

                typeMethods.r8vec_print(DIM_NUM, p1, "  Vertex A");
                typeMethods.r8vec_print(DIM_NUM, p2, "  Vertex B");
                typeMethods.r8vec_print(DIM_NUM, p3, "  Vertex C");

                angle = Angle.angle_rad_2d(p1, p2, p3);

                Console.WriteLine("");
                Console.WriteLine("  Angle = " + angle + "");
                Console.WriteLine("");
            }
        }

        public static void angle_rad_3d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_3D_TEST tests ANGLE_RAD_3D;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 3;
            int TEST_NUM = 3;

            int i;
            double[] p1 = new double[DIM_NUM];
            double[] p1_test =
            {
                1.0, 0.0, 0.0,
                1.0, 2.0, 3.0,
                0.0, 0.0, 1.0
            };
            double[] p2 = {0.0, 0.0, 0.0};
            double[] p3 = {0.0, 0.0, 1.0};
            double temp1;
            double temp2;
            int test;

            Console.WriteLine("");
            Console.WriteLine("ANGLE_RAD_3D_TEST");
            Console.WriteLine("  ANGLE_RAD_3D computes an angle;");
            Console.WriteLine("");
            Console.WriteLine("         X           Y           Z   ANGLE_RAD_3D  (Degrees)");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    p1[i] = p1_test[i + test * DIM_NUM];
                }

                temp1 = Angle.angle_rad_3d(p1, p2, p3);
                temp2 = Helpers.radians_to_degrees(temp1);

                string cout = "";
                for (i = 0; i < DIM_NUM; i++)
                {
                    cout += "  " + p1[i].ToString().PadLeft(10);
                }

                Console.WriteLine(cout + "  " + temp1.ToString().PadLeft(10)
                                  + "  " + temp2.ToString().PadLeft(10) + "");
            }
        }

        public static void angle_rad_nd_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_RAD_ND_TEST tests ANGLE_RAD_ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;

            int i;
            int angle_num = 12;
            double temp1;
            double temp2;
            double thetad;
            double thetar;
            double[] v1 = new double[DIM_NUM];
            double[] v2 = new double[DIM_NUM];

            Console.WriteLine("");
            Console.WriteLine("ANGLE_RAD_ND_TEST");
            Console.WriteLine("  ANGLE_RAD_ND computes an angle.");
            Console.WriteLine("");
            Console.WriteLine("      X           Y          Theta         atan2  ANGLE_RAD_ND");
            Console.WriteLine("");

            v1[0] = 1.0;
            v1[1] = 0.0;

            for (i = 0; i <= angle_num; i++)
            {
                thetad = (double) (i) * 360.0 / (double) (angle_num);
                thetar = Helpers.degrees_to_radians(thetad);

                v2[0] = Math.Cos(thetar);
                v2[1] = Math.Sin(thetar);

                temp1 = Helpers.radians_to_degrees(Math.Atan2(v2[1], v2[0]));

                temp2 = Angle.angle_rad_nd(DIM_NUM, v1, v2);

                Console.WriteLine("  " + v2[0].ToString().PadLeft(10)
                                       + "  " + v2[1].ToString().PadLeft(10)
                                       + "  " + thetad.ToString().PadLeft(10)
                                       + "  " + temp1.ToString().PadLeft(10)
                                       + "  " + temp2.ToString().PadLeft(10) + "");
            }

        }

        public static void angle_turn_2d_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ANGLE_TURN_2D_TEST tests ANGLE_TURN_2D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int DIM_NUM = 2;
            int TEST_NUM = 13;

            double[] p1 = new double[DIM_NUM];
            double[] p2 = {0.0, 0.0};
            double[] p3 = {1.0, 0.0};
            double pi = 3.141592653589793;
            int test;
            double theta;
            double theta_degrees;
            double turn;

            Console.WriteLine("");
            Console.WriteLine("ANGLE_TURN_2D_TEST");
            Console.WriteLine("  ANGLE_TURN_2D computes the turning angle");
            Console.WriteLine("  defined by the line segments [P1,P2] and [P2,P3].");

            Console.WriteLine("");
            Console.WriteLine("  Our three points are:");
            Console.WriteLine("");
            Console.WriteLine("    P1 = (C,S)");
            Console.WriteLine("    P2 = (0,0)");
            Console.WriteLine("    P3 = (1,0)");
            Console.WriteLine("");
            Console.WriteLine("  C = cosine ( theta ), S = sine ( theta ).");
            Console.WriteLine("");
            Console.WriteLine("  Test  Theta    Turn");
            Console.WriteLine("");

            for (test = 1; test <= TEST_NUM; test++)
            {
                theta = 2.0 * pi * (double) (test - 1)
                        / (double) (TEST_NUM - 1);

                theta_degrees = 360.0 * (double) (test - 1)
                                / (double) (TEST_NUM - 1);

                p1[0] = Math.Cos(theta);
                p1[1] = Math.Sin(theta);

                turn = Angle.angle_turn_2d(p1, p2, p3);

                Console.WriteLine("  " + test.ToString().PadLeft(4)
                                       + "  " + theta_degrees.ToString().PadLeft(5)
                                       + "  " + turn.ToString().PadLeft(14) + "");
            }

        }

    }
}