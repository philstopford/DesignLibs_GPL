using System;
using System.Globalization;
using System.Linq;
using Burkardt.Geometry;
using Burkardt.Uniform;

namespace GeometryTest;

public static class SegmentTest
{

    public static void test0418()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0418 tests SEGMENTS_CURVATURE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 13;

        double curvature;
        double[] p1 = {0.0, 0.0};
        double[] p2 = {1.0, 0.0};
        double[] p3 = new double[DIM_NUM];
        double theta;
        double theta_degrees;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0418");
        Console.WriteLine("  SEGMENTS_CURVATURE_2D computes the local curvature");
        Console.WriteLine("  defined by the line segments [P1,P2] and [P2,P3].");

        Console.WriteLine("");
        Console.WriteLine("  Our three points are:");
        Console.WriteLine("");
        Console.WriteLine("    P1 = (0,0)");
        Console.WriteLine("    P2 = (1,0)");
        Console.WriteLine("    P3 = (C,S)");
        Console.WriteLine("");
        Console.WriteLine("  C = cosine ( theta), S = sine ( theta ).");
        Console.WriteLine("");
        Console.WriteLine("  Test  Theta  Curvature");
        Console.WriteLine("");

        for (test = 1; test <= TEST_NUM; test++)
        {
            theta = 2.0 * Math.PI * (test - 1)
                    / (TEST_NUM - 1);

            theta_degrees = 360.0 * (test - 1)
                            / (TEST_NUM - 1);

            p3[0] = Math.Cos(theta);
            p3[1] = Math.Sin(theta);

            curvature = Segments.segments_curvature_2d(p1, p2, p3);

            Console.WriteLine("  " + test.ToString().PadLeft(4)
                                   + "  " + theta_degrees.ToString().PadLeft(5)
                                   + "  " + curvature.ToString().PadLeft(14) + "");
        }
    }

    public static void test042()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST042 tests SEGMENTS_DIST_2D.
        //
        //  Discussion:
        //
        //    Case 1, parallel, not coincident.
        //    Case 2, parallel, coincident, overlapping.
        //    Case 3, parallel, coincident, disjoint.
        //    Case 4, nonparallel, intersecting.
        //    Case 5, nonparallel, disjoint.
        //    Case 6 and 7, should be same, because simply a translation by 50;
        //    Case 8 and 9, answers should be same.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 9;

        double dist;
        int i;
        double[] p1;
        double[] p1_test =
        {
            2.0, 3.0,
            2.0, 3.0,
            2.0, 3.0,
            2.0, 3.0,
            2.0, 3.0,
            57.0, 53.0,
            7.0, 3.0,
            0.0, 0.0,
            -10.0, -10.0
        };
        double[] p2;
        double[] p2_test =
        {
            8.0, 6.0,
            8.0, 6.0,
            8.0, 6.0,
            8.0, 6.0,
            8.0, 6.0,
            58.0, 53.0,
            8.0, 3.0,
            100.0, 100.0,
            100.0, 100.0
        };
        double[] q1;
        double[] q1_test =
        {
            8.0, 3.0,
            4.0, 4.0,
            14.0, 9.0,
            0.0, 8.0,
            7.0, 3.0,
            65.0, 45.0,
            15.0, -5.0,
            50.0, 0.0,
            50.0, 0.0
        };
        double[] q2;
        double[] q2_test =
        {
            14.0, 6.0,
            14.0, 9.0,
            16.0, 10.0,
            5.0, 3.0,
            9.0, -1.0,
            57.0, 53.0,
            7.0, 3.0,
            60.0, 0.0,
            60.0, 0.0
        };
        int test;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST042");
        Console.WriteLine("  SEGMENTS_DIST_2D computes the distance between");
        Console.WriteLine("  line segments in 2D.");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = p1_test.Skip(+test * DIM_NUM).ToArray();
            p2 = p2_test.Skip(+test * DIM_NUM).ToArray();
            q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
            q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

            dist = Segments.segments_dist_2d(p1, p2, q1, q2);

            Console.WriteLine("");

            switch (test)
            {
                case 0:
                    Console.WriteLine("  Same slope, different intercepts.");
                    break;
                case 1:
                    Console.WriteLine("  Same slope, same intercepts, overlapping.");
                    Console.WriteLine("  Distance should be 0.");
                    break;
                case 2:
                    Console.WriteLine("  Same slope, same intercepts, disjoint.");
                    Console.WriteLine("  Distance should be sqrt(45)=6.7082038");
                    break;
                case 3:
                    Console.WriteLine("  Different slopes, intersecting.");
                    Console.WriteLine("  Distance should be 0.");
                    break;
                case 4:
                    Console.WriteLine("  Different slopes, not intersecting.");
                    break;
                case 5:
                    Console.WriteLine("  Simple problem.");
                    Console.WriteLine("  Distance should be 0.");
                    break;
                case 6:
                    Console.WriteLine("  Same data, translated by 50.");
                    Console.WriteLine("  Distance should be 0.");
                    break;
                case 7:
                    Console.WriteLine("  Diagonal and horizontal.");
                    Console.WriteLine("  Distance should be sqrt(2500/2)=35.355339");
                    break;
                case 8:
                    Console.WriteLine("  Same data, except first segment extended.");
                    Console.WriteLine("  Distance should be sqrt(2500/2)=35.355339");
                    break;
            }

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
            Console.WriteLine("  Distance = " + dist + "");
        }
    }

    public static void test043()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST043 tests SEGMENTS_DIST_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double dist;
        double[] p1 = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        double[] p3 = new double[DIM_NUM];
        double[] p4 = new double[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST043");
        Console.WriteLine("  SEGMENTS_DIST_3D computes the distance between two");
        Console.WriteLine("  line segments in 3D.");
        Console.WriteLine(" ");
        Console.WriteLine("  Case   Computed    True");
        Console.WriteLine(" ");
        //
        //  Case 0, colinear, not coincident.
        //
        //  LS1: (2,3,0) + t * (2,1,0) for t = 0 to 3.
        //  LS2: (11,6,4) + t * (2,1,0) for t = 0 to 3.
        //  Distance is 5.
        //
        p1[0] = 0.0;
        p1[1] = 0.0;
        p1[2] = 0.0;

        p2[0] = 1.0;
        p2[1] = 0.0;
        p2[2] = 0.0;

        p3[0] = 2.0;
        p3[1] = 0.0;
        p3[2] = 0.0;

        p4[0] = 4.0;
        p4[1] = 0.0;
        p4[2] = 0.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 0 + "  " + dist + "  " + 1.0 + "");
        //
        //  Case 1, parallel, not coincident.
        //
        //  LS1: (2,3,0) + t * (2,1,0) for t = 0 to 3.
        //  LS2: (11,6,4) + t * (2,1,0) for t = 0 to 3.
        //  Distance is 5.
        //
        p1[0] = 2.0;
        p1[1] = 3.0;
        p1[2] = 0.0;

        p2[0] = 8.0;
        p2[1] = 6.0;
        p2[2] = 0.0;

        p3[0] = 11.0;
        p3[1] = 6.0;
        p3[2] = 4.0;

        p4[0] = 17.0;
        p4[1] = 9.0;
        p4[2] = 4.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 1.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 5.0 + "");
        //
        //  Case 2, parallel, coincident, overlapping.
        //
        //  (1,2,3) + t * ( 1,-1,2)
        //  LS1: t = 0 to t = 3.
        //  Distance is 0.
        //
        p1[0] = 1.0;
        p1[1] = 2.0;
        p1[2] = 3.0;

        p2[0] = 4.0;
        p2[1] = -1.0;
        p2[2] = 9.0;

        p3[0] = 3.0;
        p3[1] = 0.0;
        p3[2] = 7.0;

        p4[0] = 6.0;
        p4[1] = -3.0;
        p4[2] = 13.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 2.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 0.0 + "");
        //
        //  Case 3, parallel, coincident, disjoint.
        //
        //  LS1: (3,4,5) + t * ( 2,2,1) for 0 <= t <= 2.
        //  LS2: (3,4,5) + t * ( 2,2,1) for 3 <= t <= 5.
        //  Distance = 3.
        //
        p1[0] = 3.0;
        p1[1] = 4.0;
        p1[2] = 5.0;

        p2[0] = 7.0;
        p2[1] = 8.0;
        p2[2] = 7.0;

        p3[0] = 9.0;
        p3[1] = 10.0;
        p3[2] = 8.0;

        p4[0] = 13.0;
        p4[1] = 14.0;
        p4[2] = 10.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 3.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 3.0 + "");
        //
        //  Case 4, nonparallel, could intersect, and does intersect.
        //
        //  L1: (1,1,1) + t * (0,1,2)
        //  L2: (0,2,3) + t * (1,0,0)
        //  intersect at (1,2,3)
        //  Distance is 0.
        //
        p1[0] = 1.0;
        p1[1] = 1.0;
        p1[2] = 1.0;

        p2[0] = 1.0;
        p2[1] = 4.0;
        p2[2] = 7.0;

        p3[0] = 0.0;
        p3[1] = 2.0;
        p3[2] = 3.0;

        p4[0] = 5.0;
        p4[1] = 2.0;
        p4[2] = 3.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 4.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 0.0 + "");
        //
        //  Case 5, nonparallel, could intersect, and does not intersect.
        //
        //  L1: (1,1,1) + t * (0,1,2)
        //  L2: (0,2,3) + t * (1,0,0)
        //  lines intersect at (1,2,3), line segments do not.
        //  Distance is 1.0;
        //
        p1[0] = 1.0;
        p1[1] = 1.0;
        p1[2] = 1.0;

        p2[0] = 1.0;
        p2[1] = 4.0;
        p2[2] = 7.0;

        p3[0] = 0.0;
        p3[1] = 2.0;
        p3[2] = 3.0;

        p4[0] = -5.0;
        p4[1] = 2.0;
        p4[2] = 3.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 5.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 1.0 + "");
        //
        //  Case 6, nonparallel, can not intersect, "end-to-end".
        //
        //  L1: (2,2,1) + t * (0,1,2)  0 <= t <= 5
        //  L2: (0,0,0) + t * (-1,-1,-1) 0 <= t <= 5
        //  Distance is 3.
        //
        p1[0] = 2.0;
        p1[1] = 2.0;
        p1[2] = 1.0;

        p2[0] = 2.0;
        p2[1] = 7.0;
        p2[2] = 11.0;

        p3[0] = 0.0;
        p3[1] = 0.0;
        p3[2] = 0.0;

        p4[0] = -5.0;
        p4[1] = -5.0;
        p4[2] = -5.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 6.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 3.0 + "");
        //
        //  Case 7, nonparallel, can not intersect, "end-to-mid".
        //
        //  L1: (1,1,1) + t * (0,1,2) 0 <= t <= 5
        //  L2: (0,4,7) + t * (-1,0,0) 0 <= t <= 5
        //  Distance is 1.
        //
        p1[0] = 1.0;
        p1[1] = 1.0;
        p1[2] = 1.0;

        p2[0] = 1.0;
        p2[1] = 6.0;
        p2[2] = 11.0;

        p3[0] = 0.0;
        p3[1] = 4.0;
        p3[2] = 7.0;

        p4[0] = -5.0;
        p4[1] = 4.0;
        p4[2] = 7.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 7.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 1.0 + "");
        //
        //  Case 8, nonparallel, can not intersect, "mid-to-mid".
        //
        //  L1: (0,5,10) + t * (1,-1,0) 0 <= t <= 5
        //  L2: (0,0,0) + t * (1,1,0) 0 <= t <= 6
        //  Distance = 10.
        //
        p1[0] = 0.0;
        p1[1] = 5.0;
        p1[2] = 10.0;

        p2[0] = 5.0;
        p2[1] = 0.0;
        p2[2] = 10.0;

        p3[0] = 0.0;
        p3[1] = 0.0;
        p3[2] = 0.0;

        p4[0] = 6.0;
        p4[1] = 6.0;
        p4[2] = 0.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 8.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 10.0 + "");
        //
        //  Case 9, nonparallel, can not intersect, "mid-to-end".
        //
        //  L1: (-2,0,0) + t * (1,0,0) 0 <= t <= 12
        //  L2: (-2,8,1) + t * (9,-4,-1) 0 <= t <= 1
        //  Distance = 4.
        //
        p1[0] = -2.0;
        p1[1] = 0.0;
        p1[2] = 0.0;

        p2[0] = 10.0;
        p2[1] = 0.0;
        p2[2] = 0.0;

        p3[0] = -2.0;
        p3[1] = 8.0;
        p3[2] = 1.0;

        p4[0] = 7.0;
        p4[1] = 4.0;
        p4[2] = 0.0;

        dist = Segments.segments_dist_3d(p1, p2, p3, p4);

        Console.WriteLine("  " + 9.ToString().PadLeft(3)
                               + "  " + dist.ToString().PadLeft(12) + "  " + 4.0 + "");

    }

    public static void test044()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST044 tests SEGMENTS_INT_1D.
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
        int TEST_NUM = 7;

        double dist;
        int test;
        double p1;
        double p2;
        double q1;
        double[] q1_test = {-1.0, 3.0, 1.0, 0.5, 0.25, 0.5, 2.0};
        double q2;
        double[] q2_test = {1.0, 2.0, 2.0, -3.0, 0.50, 0.5, 2.0};
        double r1 = 0;
        double r2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST044");
        Console.WriteLine("  SEGMENTS_INT_1D determines the intersection [R1,R2]");
        Console.WriteLine("  of line segments [P1,P2] and [Q1,Q2] in 1D.");
        Console.WriteLine("");
        Console.WriteLine("  DIST is negative for overlap,");
        Console.WriteLine("  0 for point intersection,");
        Console.WriteLine("  positive if there is no overlap.");
        Console.WriteLine("");
        Console.WriteLine("  Test        P1        P2        Q1        Q2        R1        R2       DIST");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = -1.0;
            p2 = 1.0;
            q1 = q1_test[test];
            q2 = q2_test[test];

            dist = Segments.segments_int_1d(p1, p2, q1, q2, ref r1, ref r2);

            Console.WriteLine("  " + test.ToString().PadLeft(2)
                                   + "  " + p1.ToString().PadLeft(8)
                                   + "  " + p2.ToString().PadLeft(8)
                                   + "  " + q1.ToString().PadLeft(8)
                                   + "  " + q2.ToString().PadLeft(8)
                                   + "  " + r1.ToString().PadLeft(8)
                                   + "  " + r2.ToString().PadLeft(8)
                                   + "  " + dist.ToString().PadLeft(8) + "");
        }
    }

    public static void test045()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST045 tests SEGMENTS_INT_2D.
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
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        int flag = 0;
        int i;
        double[] p1 = {-1.0, 3.0};
        double[] p2 = {1.0, 1.0};
        double[] q1;
        double[] q1_test =
        {
            -1.0, 1.0,
            3.0, -1.0,
            0.0, 0.0,
            1.0, 2.0
        };
        double[] q2;
        double[] q2_test =
        {
            1.0, -1.0,
            2.0, 0.0,
            0.0, 9.0,
            3.0, 2.0
        };
        double[] r = new double[DIM_NUM];
        int test;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST045");
        Console.WriteLine("  SEGMENTS_INT_2D searches for an intersection of two");
        Console.WriteLine("  line segments in 2D.");
        Console.WriteLine("");
        Console.WriteLine("  All tests use the same line segment 1:");
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

        for (test = 0; test < TEST_NUM; test++)
        {
            q1 = q1_test.Skip(+test * DIM_NUM).ToArray();
            q2 = q2_test.Skip(+test * DIM_NUM).ToArray();

            Console.WriteLine("");
            Console.WriteLine("");
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

            Segments.segments_int_2d(p1, p2, q1, q2, ref flag, ref r);

            switch (flag)
            {
                case 0:
                    Console.WriteLine("");
                    Console.WriteLine("  The line segments do not intersect.");
                    break;
                case 1:
                {
                    Console.WriteLine("");
                    cout = "  The line segments intersect at:";
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + r[i].ToString().PadLeft(12);
                    }

                    Console.WriteLine(cout);
                    break;
                }
            }
        }
    }

    public static void test036()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST036 tests SEGMENT_CONTAINS_POINT_1D.
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
        int TEST_NUM = 4;

        double p;
        double[] p_test = {3.0, 7.5, 20.0, 5.0};
        double p1;
        double[] p1_test = {2.0, 10.0, 8.0, 88.0};
        double p2;
        double[] p2_test = {6.0, -10.0, 10.0, 88.0};
        double t = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST036");
        Console.WriteLine("  SEGMENT_CONTAINS_POINT_1D determines if a point");
        Console.WriteLine("  lies within a line segment in 1D.");
        Console.WriteLine("");
        Console.WriteLine("       P1     P       T");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = p1_test[test];
            p2 = p2_test[test];
            p = p_test[test];

            Segments.segment_contains_point_1d(p1, p2, p, ref t);

            Console.WriteLine("  " + p1.ToString().PadLeft(7)
                                   + "  " + p2.ToString().PadLeft(7)
                                   + "  " + p.ToString().PadLeft(7)
                                   + "  " + t.ToString().PadLeft(12) + "");
        }

    }

    public static void test0365()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0365 tests SEGMENT_POINT_DIST_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double dist;
        int i;
        double[] p;
        double[] p1;
        double[] p2;
        int seed = 123456789;
        int test;
        int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST0365");
        Console.WriteLine("  SEGMENT_POINT_DIST_2D computes the distance");
        Console.WriteLine("  between a line segment and point in 2D");

        for (test = 1; test <= test_num; test++)
        {
            p1 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p2 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

            dist = Segments.segment_point_dist_2d(p1, p2, p);

            Console.WriteLine("");
            Console.WriteLine("  TEST = " + test.ToString().PadLeft(2) + "");
            string cout = "  P1 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p1[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P2 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p2[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P =    ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  DIST = " + dist.ToString().PadLeft(12) + "");
        }

    }

    public static void segment_point_dist_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_DIST_3D_TEST tests SEGMENT_POINT_DIST_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double dist;
        int i;
        double[] p;
        double[] p1;
        double[] p2;
        int seed = 123456789;
        int test;
        int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("SEGMENT_POINT_DIST_3D_TEST");
        Console.WriteLine("  SEGMENT_POINT_DIST_3D computes the distance");
        Console.WriteLine("  between a line segment and point in 3D");

        for (test = 1; test <= test_num; test++)
        {
            p1 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p2 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

            dist = Segments.segment_point_dist_3d(p1, p2, p);

            Console.WriteLine("");
            Console.WriteLine("  TEST = " + test.ToString().PadLeft(12) + "");
            string cout = "  P1 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p1[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P2 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p2[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P =    ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  DIST = " + dist.ToString().PadLeft(12) + "");
        }

    }

    public static void segment_point_near_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR_2D_TEST tests SEGMENT_POINT_NEAR_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double dist = 0;
        int i;
        double[] p;
        double[] p1;
        double[] p2;
        double[] pn = new double[DIM_NUM];
        int seed = 123456789;
        double t = 0;
        int test;
        int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("SEGMENT_POINT_NEAR_2D_TEST");
        Console.WriteLine("  SEGMENT_POINT_NEAR_2D computes the nearest point");
        Console.WriteLine("  on a line segment to a point in 2D");

        for (test = 1; test <= test_num; test++)
        {
            p1 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p2 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

            Segments.segment_point_near_2d(p1, p2, p, ref pn, ref dist, ref t);

            Console.WriteLine("");
            Console.WriteLine("  TEST = " + test.ToString().PadLeft(2) + "");
            string cout = "  P1 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p1[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P2 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p2[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P =    ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  PN =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += pn[i].ToString().PadLeft(12);
            }

            Console.WriteLine("");
            Console.WriteLine("  DIST = " + dist.ToString().PadLeft(12) + "");
            Console.WriteLine("  T =    " + t.ToString().PadLeft(12) + "");
        }

    }

    public static void segment_point_near_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR_3D_TEST tests SEGMENT_POINT_NEAR_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double dist = 0;
        int i;
        double[] p;
        double[] p1;
        double[] p2;
        double[] pn = new double[DIM_NUM];
        int seed = 123456789;
        double t = 0;
        int test;
        int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("SEGMENT_POINT_NEAR_3D_TEST");
        Console.WriteLine("  SEGMENT_POINT_NEAR_3D computes the nearest point");
        Console.WriteLine("  on a line segment to a point in 3D");

        for (test = 1; test <= test_num; test++)
        {
            p1 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p2 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            p = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);

            Segments.segment_point_near_3d(p1, p2, p, ref pn, ref dist, ref t);

            Console.WriteLine("");
            Console.WriteLine("  TEST = " + test.ToString().PadLeft(2) + "");
            string cout = "  P1 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p1[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P2 =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p2[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  P =    ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
            cout = "  PN =   ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += pn[i].ToString().PadLeft(12);
            }

            Console.WriteLine("");
            Console.WriteLine("  DIST = " + dist.ToString().PadLeft(12) + "");
            Console.WriteLine("  T    = " + t.ToString().PadLeft(12) + "");
        }

    }

    public static void segment_point_near_3d_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEGMENT_POINT_NEAR_3D_TEST2 tests SEGMENT_POINT_NEAR_3D.
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
        int DIM_NUM = 3;

        double dist = 0;
        double[] p = new double[DIM_NUM];
        double[] p1 = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        double[] pn = new double[DIM_NUM];
        double t = 0;

        Console.WriteLine("");
        Console.WriteLine("SEGMENT_POINT_NEAR_3D_TEST2");
        Console.WriteLine("  SEGMENT_POINT_NEAR_3D computes the nearest point on a");
        Console.WriteLine("  line segment, to a given point, in 3 space.");
        Console.WriteLine("");
        Console.WriteLine("  Case  T  Distance    PN");
        Console.WriteLine("");
        //
        //  Case 1, point is nearest end of segment.
        //
        //  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
        //  P (11,6,4)
        //  Distance is 5.
        //
        p1[0] = 2.0;
        p1[1] = 3.0;
        p1[2] = 0.0;

        p2[0] = 8.0;
        p2[1] = 6.0;
        p2[2] = 0.0;

        p[0] = 11.0;
        p[1] = 6.0;
        p[2] = 4.0;

        Segments.segment_point_near_3d(p1, p2, p, ref pn, ref dist, ref t);

        Console.WriteLine("  " + 1.ToString().PadLeft(6)
                               + "  " + t.ToString().PadLeft(10)
                               + "  " + dist.ToString().PadLeft(10)
                               + "  " + pn[0].ToString().PadLeft(10)
                               + "  " + pn[1].ToString().PadLeft(10)
                               + "  " + pn[2].ToString().PadLeft(10) + "");
        //
        //  Case 2, point is nearest interior point of segment.
        //
        //  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
        //  P (4,4,1)
        //  Distance is 1.
        //
        p1[0] = 2.0;
        p1[1] = 3.0;
        p1[2] = 0.0;

        p2[0] = 8.0;
        p2[1] = 6.0;
        p2[2] = 0.0;

        p[0] = 4.0;
        p[1] = 4.0;
        p[2] = 1.0;

        Segments.segment_point_near_3d(p1, p2, p, ref pn, ref dist, ref t);

        Console.WriteLine("  " + 2.ToString().PadLeft(6)
                               + "  " + t.ToString().PadLeft(10)
                               + "  " + dist.ToString().PadLeft(10)
                               + "  " + pn[0].ToString().PadLeft(10)
                               + "  " + pn[1].ToString().PadLeft(10)
                               + "  " + pn[2].ToString().PadLeft(10) + "");
        //
        //  Case 3, point is on the line.
        //
        //  LS: (2,3,0) + t * (2,1,0) for t = 0 to 3.
        //  P (6,5,0)
        //  Distance is 0.
        //
        p1[0] = 2.0;
        p1[1] = 3.0;
        p1[2] = 0.0;

        p2[0] = 8.0;
        p2[1] = 6.0;
        p2[2] = 0.0;

        p[0] = 6.0;
        p[1] = 5.0;
        p[2] = 0.0;

        Segments.segment_point_near_3d(p1, p2, p, ref pn, ref dist, ref t);

        Console.WriteLine("  " + 3.ToString().PadLeft(6)
                               + "  " + t.ToString().PadLeft(10)
                               + "  " + dist.ToString().PadLeft(10)
                               + "  " + pn[0].ToString().PadLeft(10)
                               + "  " + pn[1].ToString().PadLeft(10)
                               + "  " + pn[2].ToString().PadLeft(10) + "");

    }

    public static void test201 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST201 tests STRING_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int VEC_NUM = 15;

        int i;
        int jstrng;
        int[] order = new int[VEC_NUM];
        double[] p1 = {
            0.0, 0.0,
            3.0, 4.0,
            2.0, 2.0,
            3.0, 2.0,
            2.0, 1.0,
            1.0, 1.0,
            0.0, 5.0,
            1.0, 2.0,
            3.0, 2.0,
            0.0, 0.0,
            5.0, 5.0,
            3.0, 3.0,
            2.0, 4.0,
            7.0, 4.0,
            1.0, 0.0 };
        double[] p2 = {
            1.0, 1.0,
            2.0, 4.0,
            1.0, 3.0,
            2.0, 3.0,
            2.0, 2.0,
            1.0, 2.0,
            1.0, 6.0,
            1.0, 3.0,
            3.0, 3.0,
            1.0, 0.0,
            6.0, 6.0,
            3.0, 4.0,
            2.0, 3.0,
            5.0, 5.0,
            2.0, 1.0 };
        int[] string_ = new int[VEC_NUM];
        int string_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST201");
        Console.WriteLine("  STRING_2D takes a set of line segments, and");
        Console.WriteLine("  strings them together.");
        Console.WriteLine("");
        Console.WriteLine("  I     P1     P2");
        Console.WriteLine("");
        for ( i = 0; i < VEC_NUM; i++ )
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + p1[0+i*2].ToString().PadLeft(10)
                                   + "  " + p1[1+i*2].ToString().PadLeft(10)
                                   + "  " + p2[0+i*2].ToString().PadLeft(10)
                                   + "  " + p1[1+i*2].ToString().PadLeft(10) + "");
        }

        Segments.string_2d ( VEC_NUM, p1, p2, ref string_num, ref order, ref string_ );

        Console.WriteLine("");
        Console.WriteLine("  Found " + string_num + " groups of segments.");
        Console.WriteLine("");
        Console.WriteLine("  STRING, ORDER, P1, P2");
        Console.WriteLine("");

        jstrng = 1;

        for ( i = 0; i < VEC_NUM; i++ )
        {
            if ( jstrng < string_[i] )
            {
                Console.WriteLine("");
                jstrng += 1;
            }
            Console.WriteLine("  " + string_[i].ToString().PadLeft(3)
                                   + "  " + order[i].ToString().PadLeft(3)
                                   + "  " + p1[0+i*2].ToString().PadLeft(10)
                                   + "  " + p1[1+i*2].ToString().PadLeft(10)
                                   + "  " + p2[0+i*2].ToString().PadLeft(10)
                                   + "  " + p2[1+i*2].ToString().PadLeft(10) + "");
        }

    }
}