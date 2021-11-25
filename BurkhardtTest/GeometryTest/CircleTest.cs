using System;
using System.Globalization;
using System.Linq;
using Burkardt.CircleNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class CircleTest
{
    public static void circle_dia2imp_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_DIA2IMP_2D_TEST tests CIRCLE_DIA2IMP_2D.
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
        const int DIM_NUM = 2;

        double[] pc = new double[DIM_NUM];
        double[] p1 = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        double r = 0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_DIA2IMP_2D_TEST");
        Console.WriteLine("  CIRCLE_DIA2IMP_2D converts a diameter to an");
        Console.WriteLine("  implicit circle in 2D.");

        const double theta = 2.0;

        p1[0] = 2.0 + 5.0 * Math.Cos(theta);
        p1[1] = 2.0 + 5.0 * Math.Sin(theta);
        p2[0] = 2.0 - 5.0 * Math.Cos(theta);
        p2[1] = 2.0 - 5.0 * Math.Sin(theta);

        typeMethods.r8vec_print(DIM_NUM, p1, "  P1:");
        typeMethods.r8vec_print(DIM_NUM, p2, "  P2:");

        Geometry.circle_dia2imp_2d(p1, p2, ref r, ref pc);

        Geometry.circle_imp_print_2d(r, pc, "  The implicit circle:");

    }

    public static void circle_exp_contains_point_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_EXP_CONTAINS_POINT_2D_TEST tests CIRCLE_EXP_CONTAINS_POINT_2D.
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
        double[] p1 = new double[2];
        double[] p2 = new double[2];
        double[] p3 = new double[2];
        double[] p4 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_EXP_CONTAINS_POINT_2D_TEST");
        Console.WriteLine("  CIRCLE_EXP_CONTAINS_POINT_2D determines if a");
        Console.WriteLine("    point lies inside a circle.");
        Console.WriteLine("");
        Console.WriteLine("  Possible return values are:");
        Console.WriteLine("");
        Console.WriteLine("  -1: The point is inside the circle.");
        Console.WriteLine("   0: The point is on the circle.");
        Console.WriteLine("   1: The point is outside the circle");
        Console.WriteLine("   2: Colinear data, the point is on the line.");
        Console.WriteLine("   3: Colinear data, the point is not on the line.");
        Console.WriteLine("   4: Two equal data points, the point is on the line.");
        Console.WriteLine("   5: Two equal data points, the point is not on the line.");
        Console.WriteLine("   6: All data points equal, the point is equal.");
        Console.WriteLine("   7: All data points equal, the point is not equal.");
        //
        //  This point is inside.
        //
        p1[0] = 4.0;
        p1[1] = 2.0;

        p2[0] = 1.0;
        p2[1] = 5.0;

        p3[0] = -2.0;
        p3[1] = 2.0;

        p4[0] = 2.0;
        p4[1] = 3.0;

        int inside = Geometry.circle_exp_contains_point_2d(p1, p2, p3, p4);

        Console.WriteLine("");
        Console.WriteLine("  P1 = " + p1[0] + "  " + p1[1] + "");
        Console.WriteLine("  P2 = " + p2[0] + "  " + p2[1] + "");
        Console.WriteLine("  P3 = " + p3[0] + "  " + p3[1] + "");
        Console.WriteLine("  P4 = " + p4[0] + "  " + p4[1] + "");
        Console.WriteLine("  INSIDE = " + inside + "");
        //
        //  This point is actually right on the circle.
        //
        p1[0] = 4.0;
        p1[1] = 2.0;

        p2[0] = 1.0;
        p2[1] = 5.0;

        p3[0] = -2.0;
        p3[1] = 2.0;

        p4[0] = 1.0;
        p4[1] = -1.0;

        inside = Geometry.circle_exp_contains_point_2d(p1, p2, p3, p4);

        Console.WriteLine("");
        Console.WriteLine("  P1 = " + p1[0] + "  " + p1[1] + "");
        Console.WriteLine("  P2 = " + p2[0] + "  " + p2[1] + "");
        Console.WriteLine("  P3 = " + p3[0] + "  " + p3[1] + "");
        Console.WriteLine("  P4 = " + p4[0] + "  " + p4[1] + "");
        Console.WriteLine("  INSIDE = " + inside + "");
        //
        //  This point is outside.
        //
        p1[0] = 4.0;
        p1[1] = 2.0;

        p2[0] = 1.0;
        p2[1] = 5.0;

        p3[0] = -2.0;
        p3[1] = 2.0;

        p4[0] = 4.0;
        p4[1] = 6.0;

        inside = Geometry.circle_exp_contains_point_2d(p1, p2, p3, p4);

        Console.WriteLine("");
        Console.WriteLine("  P1 = " + p1[0] + "  " + p1[1] + "");
        Console.WriteLine("  P2 = " + p2[0] + "  " + p2[1] + "");
        Console.WriteLine("  P3 = " + p3[0] + "  " + p3[1] + "");
        Console.WriteLine("  P4 = " + p4[0] + "  " + p4[1] + "");
        Console.WriteLine("  INSIDE = " + inside + "");

    }

    public static void circle_exp2imp_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_EXP2IMP_2D_TEST tests CIRCLE_EXP2IMP_2D.
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
        const int DIM_NUM = 2;
        const int TEST_NUM = 3;

        double[] pc = new double[DIM_NUM];
        double r = 0;
        int test;
        double[] test_t =
        {
            4.0, 2.0,
            1.0, 5.0,
            -2.0, 2.0,
            4.0, 2.0,
            5.0, 4.0,
            6.0, 6.0,
            4.0, 2.0,
            1.0, 5.0,
            4.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_EXP2IMP_2D_TEST");
        Console.WriteLine("  CIRCLE_EXP2IMP_2D computes the radius and");
        Console.WriteLine("  center of the circle through three points.");

        for (test = 0; test < TEST_NUM; test++)
        {
            double[] t = test_t.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  The points::");

            double[] p1 = t;
            double[] p2 = p1.Skip(DIM_NUM).ToArray();
            double[] p3 = p2.Skip(DIM_NUM).ToArray();

            Geometry.circle_exp2imp_2d(p1, p2, p3, ref r, ref pc);

            Geometry.circle_imp_print_2d(r, pc, "  The implicit circle:");
        }
    }

    public static void circle_imp_point_dist_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINT_DIST_2D_TEST tests CIRCLE_IMP_POINT_DIST_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center = {0.0, 0.0};
        int i;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_IMP_POINT_DIST_2D_TEST");
        Console.WriteLine("  CIRCLE_IMP_POINT_DIST_2D finds the");
        Console.WriteLine("  distance from a point to a circle.");

        double r = 5.0;
        Geometry.circle_imp_print_2d(r, center, "  The circle:");

        Console.WriteLine("");
        Console.WriteLine("       X       Y       D");
        Console.WriteLine("");

        int seed = 123456789;

        for (i = 1; i <= 10; i++)
        {
            double[] p = UniformRNG.r8vec_uniform_ab_new(2, -10.0, +10.0, ref seed);
            double d = Geometry.circle_imp_point_dist_2d(r, center, p);
            Console.WriteLine("  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + d.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    public static void circle_imp_points_arc_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_IMP_POINTS_ARC_2D_TEST tests CIRCLE_IMP_POINTS_ARC_2D.
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
        const int N = 13;
        const int DIM_NUM = 2;

        const double r = 2.0;
        double[] p = new double[DIM_NUM * N];
        double[] pc = {5.0, -2.0};

        const double theta1 = Math.PI / 2.0;
        const double theta2 = 3.0 * Math.PI / 2.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_IMP_POINTS_ARC_2D_TEST");
        Console.WriteLine("  CIRCLE_IMP_POINTS_ARC_2D returns points on a");
        Console.WriteLine("  circular arc.");
        Console.WriteLine("");
        Console.WriteLine("  The circle will have center " + pc[0] + ", " + pc[1] + "");
        Console.WriteLine("  and radius R = " + r + "");
        Console.WriteLine("");
        Console.WriteLine("  The arc extends from THETA1 = " + theta1 + "");
        Console.WriteLine("  to THETA2 = " + theta2 + "");

        Geometry.circle_imp_points_arc_2d(r, pc, theta1, theta2, N, ref p);

        typeMethods.r8mat_transpose_print(DIM_NUM, N, p, "  Sample results:");
    }

    public static void circle_llr2imp_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LLR2IMP_2D_TEST tests CIRCLE_LLR2IMP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        const double p_hi = 10.0;
        const double p_lo = -10.0;
        int test;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LLR2IMP_2D_TEST");
        Console.WriteLine("  CIRCLE_LLR2IMP_2D is given:");
        Console.WriteLine("  a line through P1 and P2,");
        Console.WriteLine("  a line through Q1 and Q2,");
        Console.WriteLine("  and a radius R,");
        Console.WriteLine("  and determines the centers C of 4 circles");
        Console.WriteLine("  of the given radius, tangent to both lines.");

        int seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            double[] p1 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] p2 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] q1 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] q2 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);

            double r_lo = 1.0;
            double r_hi = 5.0;
            double r = UniformRNG.r8_uniform_ab(r_lo, r_hi, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Radius R = " + r + "");

            Console.WriteLine("  Point #P1: ( "
                              + p1[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p1[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #P2: ( "
                              + p2[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p2[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #Q1: ( "
                              + q1[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + q1[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #Q2: ( "
                              + q2[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + q2[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");

            double[] pc = Geometry.circle_llr2imp_2d(p1, p2, q1, q2, r);

            Console.WriteLine("  Center #1: ( "
                              + pc[0 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Center #2: ( "
                              + pc[0 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Center #3: ( "
                              + pc[0 + 2 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 2 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Center #4: ( "
                              + pc[0 + 3 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 3 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            //
            //  Check that the lines are the right distance from the center.
            //
            double d1 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(p1, p2, pc);
            double d2 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(p1, p2, pc, +1 * 2);
            double d3 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(p1, p2, pc, +2 * 2);
            double d4 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(p1, p2, pc, +3 * 2);

            Console.WriteLine("  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d3.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d4.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

            d1 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(q1, q2, pc);
            d2 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(q1, q2, pc, +1 * 2);
            d3 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(q1, q2, pc, +2 * 2);
            d4 = Burkardt.LineNS.Geometry.line_exp_point_dist_2d(q1, q2, pc, +3 * 2);

            Console.WriteLine("  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d3.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d4.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }

    }

    public static void circle_lune_angle_by_height_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D_TEST tests CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const int n_test = 6;

        const double r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D_TEST");
        Console.WriteLine("  CIRCLE_LUNE_ANGLE_BY_HEIGHT_2D computes the angle of a");
        Console.WriteLine("  circular lune based on the 'height' of the circular triangle.");
        Console.WriteLine("");
        Console.WriteLine("      R            H        Angle");
        Console.WriteLine("");

        for (i = -n_test; i <= n_test; i++)
        {
            double h = i * r / n_test;

            double angle = Geometry.circle_lune_angle_by_height_2d(r, h);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + angle.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void circle_lune_area_by_angle_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_AREA_BY_ANGLE_2D_TEST tests CIRCLE_LUNE_AREA_BY_ANGLE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] pc = {0.0, 0.0};
        const double r = 1.0;
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LUNE_AREA_BY_ANGLE_2D_TEST");
        Console.WriteLine("  CIRCLE_LUNE_AREA_BY_ANGLE_2D computes the area of a");
        Console.WriteLine("  circular lune, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc.");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        Area");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            const double theta1 = 0.0;
            double theta2 = test * 2.0 * Math.PI / test_num;

            double area = Geometry.circle_lune_area_by_angle_2d(r, pc, theta1, theta2);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void circle_lune_area_by_height_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_AREA_BY_HEIGHT_2D_TEST tests CIRCLE_LUNE_AREA_BY_HEIGHT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n_test = 6;
        const double r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LUNE_AREA_BY_HEIGHT_2D_TEST");
        Console.WriteLine("  CIRCLE_LUNE_AREA_BY_HEIGHT_2D computes the area of a");
        Console.WriteLine("  circular lune, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc.");
        Console.WriteLine("");
        Console.WriteLine("      R            Height        Area");
        Console.WriteLine("");

        for (i = -n_test; i <= n_test; i++)
        {
            double h = i * r / n_test;

            double area = Geometry.circle_lune_area_by_height_2d(r, h);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void circle_lune_centroid_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_CENTROID_2D_TEST tests CIRCLE_LUNE_CENTROID_2D.
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
        double[] pc = {5.0, 3.0};
        const double r = 2.0;
        int test;
        const int test_num = 12;
        const double theta1 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LUNE_CENTROID_2D_TEST");
        Console.WriteLine("  CIRCLE_LUNE_CENTROID_2D computes the centroid of a");
        Console.WriteLine("  circular lune, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc.");

        Geometry.circle_imp_print_2d(r, pc, "  The implicit circle:");

        Console.WriteLine("");
        Console.WriteLine("  The first angle of our lune is always 0.");
        Console.WriteLine("");
        Console.WriteLine("  THETA2           X             Y");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            double theta2 = test * 2.0 * Math.PI / test_num;

            double[] centroid = Geometry.circle_lune_centroid_2d(r, pc, theta1, theta2);

            Console.WriteLine("  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + centroid[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + centroid[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }

    }

    public static void circle_lune_height_by_angle_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D_TEST tests CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        const int n_test = 12;

        const double r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D_TEST");
        Console.WriteLine("  CIRCLE_LUNE_HEIGHT_BY_ANGLE_2D computes the height of");
        Console.WriteLine("  the triangle of a circular lune, given the subtended angle.");
        Console.WriteLine("");
        Console.WriteLine("      R            Angle        Height");
        Console.WriteLine("");

        for (i = 0; i <= n_test; i++)
        {
            double angle = i * 2.0 * Math.PI / n_test;

            double height = Geometry.circle_lune_height_by_angle_2d(r, angle);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + angle.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + height.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void circle_pppr2imp_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_PPPR2IMP_3D_TEST tests CIRCLE_PPPR2IMP_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 3;

        double[] normal = new double[DIM_NUM];
        const double p_hi = 10.0;
        const double p_lo = -10.0;
        double[] pc = new double[DIM_NUM * 2];
        int test;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_PPPR2IMP_3D_TEST");
        Console.WriteLine("  CIRCLE_PPPR2IMP_3D is given 3D points P1, P2, P3,");
        Console.WriteLine("  and a radius R,");
        Console.WriteLine("  and determines the centers C of two circles");
        Console.WriteLine("  of the given radius, passing through P1 and P2");
        Console.WriteLine("  and lying in the plane of P1, P2 and P3.");

        int seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            double[] p1 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] p2 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] p3 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);

            double r_lo = typeMethods.r8vec_distance(DIM_NUM, p1, p2);
            double r_hi = r_lo + 5.0;
            double r = UniformRNG.r8_uniform_ab(r_lo, r_hi, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Radius R = " + r + "");

            Console.WriteLine("  Point #1: ( "
                              + p1[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p1[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #2: ( "
                              + p2[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p2[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #3: ( "
                              + p3[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p3[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");

            Geometry.circle_pppr2imp_3d(p1, p2, p3, r, ref pc, ref normal);

            Console.WriteLine("  Center #1: ( "
                              + pc[0 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Center #2: ( "
                              + pc[0 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            //
            //  Check that the points are the right distance from the center.
            //
            double d11 = typeMethods.r8vec_distance(DIM_NUM, p1, pc);
            double d21 = typeMethods.r8vec_distance(DIM_NUM, p2, pc);
            double d12 = typeMethods.r8vec_distance(DIM_NUM, p1, pc, +1 * DIM_NUM);
            double d22 = typeMethods.r8vec_distance(DIM_NUM, p2, pc, +1 * DIM_NUM);

            Console.WriteLine("  " + d11.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d21.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d12.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d22.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
            //
            //  Check that the radial vector to the point is perpendicular to NORMAL.
            //
            d11 = 0.0;
            d21 = 0.0;
            d12 = 0.0;
            d22 = 0.0;
            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                d11 += normal[i] * (p1[i] - pc[i + 0 * DIM_NUM]);
                d21 += normal[i] * (p2[i] - pc[i + 0 * DIM_NUM]);
                d12 += normal[i] * (p1[i] - pc[i + 1 * DIM_NUM]);
                d22 += normal[i] * (p2[i] - pc[i + 1 * DIM_NUM]);
            }

            Console.WriteLine("  " + d11.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d21.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d12.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d22.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    public static void circle_ppr2imp_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_PPR2IMP_2D_TEST tests CIRCLE_PPR2IMP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        const double p_hi = 10.0;
        const double p_lo = -10.0;
        int test;
        const int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_PPR2IMP_2D_TEST");
        Console.WriteLine("  CIRCLE_PPR2IMP_2D is given points P1 and P2,");
        Console.WriteLine("  and a radius R,");
        Console.WriteLine("  and determines the centers C of two circles");
        Console.WriteLine("  of the given radius, passing through P1 and P2.");

        int seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            double[] p1 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);
            double[] p2 = UniformRNG.r8vec_uniform_ab_new(DIM_NUM, p_lo, p_hi, ref seed);

            double r_lo = typeMethods.r8vec_distance(DIM_NUM, p1, p2);
            double r_hi = r_lo + 5.0;
            double r = UniformRNG.r8_uniform_ab(r_lo, r_hi, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Radius R = " + r + "");

            Console.WriteLine("  Point #1: ( "
                              + p1[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p1[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Point #2: ( "
                              + p2[0].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + p2[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");

            double[] pc = Geometry.circle_ppr2imp_2d(p1, p2, r);

            Console.WriteLine("  Center #1: ( "
                              + pc[0 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");
            Console.WriteLine("  Center #2: ( "
                              + pc[0 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ","
                              + pc[1 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12) + ")");

            double d11 = typeMethods.r8vec_distance(DIM_NUM, p1, pc);
            double d21 = typeMethods.r8vec_distance(DIM_NUM, p2, pc);
            double d12 = typeMethods.r8vec_distance(DIM_NUM, p1, pc, +1 * DIM_NUM);
            double d22 = typeMethods.r8vec_distance(DIM_NUM, p2, pc, +1 * DIM_NUM);

            Console.WriteLine("  " + d11.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d21.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d12.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d22.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    public static void circle_sector_area_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_AREA_2D_TEST tests CIRCLE_SECTOR_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] pc = {0.0, 0.0};
        const double r = 1.0;
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SECTOR_AREA_2D_TEST");
        Console.WriteLine("  CIRCLE_SECTOR_AREA_2D computes the area of a");
        Console.WriteLine("  circular sector, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc to the center.");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        Area");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            double theta1 = 0.0;
            double theta2 = test * 2.0 * Math.PI / test_num;

            double area = Geometry.circle_sector_area_2d(r, pc, theta1, theta2);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void circle_sector_centroid_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_SECTOR_CENTROID_2D_TEST tests CIRCLE_SECTOR_CENTROID_2D.
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
        double[] pc = {5.0, 3.0};
        const double r = 2.0;
        int test;
        const int test_num = 12;
        const double theta1 = 0.0;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_SECTOR_CENTROID_2D_TEST");
        Console.WriteLine("  CIRCLE_SECTOR_CENTROID_2D computes the centroid of a");
        Console.WriteLine("  circular sector, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc to the center.");

        Geometry.circle_imp_print_2d(r, pc, "  The implicit circle:");

        Console.WriteLine("");
        Console.WriteLine("  The first angle of our sector is always 0.");
        Console.WriteLine("");
        Console.WriteLine("  THETA2           X             Y");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            double theta2 = test * 2.0 * Math.PI / test_num;

            double[] centroid = Geometry.circle_sector_centroid_2d(r, pc, theta1, theta2);

            Console.WriteLine("  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + centroid[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + centroid[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }

    }

    public static void circle_triangle_area_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_TRIANGLE_AREA_2D_TEST tests CIRCLE_TRIANGLE_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] pc = {0.0, 0.0};
        const double r = 1.0;
        int test;
        const int test_num = 12;

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_TRIANGLE_AREA_2D_TEST");
        Console.WriteLine("  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a");
        Console.WriteLine("  triangle, defined by joining the endpoints");
        Console.WriteLine("  of a circular arc and the center.");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        Area");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            double theta1 = 0.0;
            double theta2 = test * 2.0 * Math.PI / test_num;

            double area = Geometry.circle_triangle_area_2d(r, pc, theta1, theta2);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }


    public static void test0155()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0155 tests CIRCLE_EXP2IMP_2D.
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
        const int DIM_NUM = 2;
        const int TEST_NUM = 13;

        double[] pc = new double[DIM_NUM];
        double[] p1 = {0.0, 0.0};
        double[] p2 = {1.0, 0.0};
        double[] p3 = new double[DIM_NUM];
        double r = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0155");
        Console.WriteLine("  CIRCLE_EXP2IMP_2D computes the radius and");
        Console.WriteLine("  center of the circle through three points.");
        Console.WriteLine("");
        Console.WriteLine("  We can use this routine to compute, for three");
        Console.WriteLine("  points in space, the circle incident to those");
        Console.WriteLine("  points, and hence the radius of that circle,");
        Console.WriteLine("  and hence the \"curvature\" of those points.");

        Console.WriteLine("");
        Console.WriteLine("  Our three points are:");
        Console.WriteLine("");
        Console.WriteLine("    (0,0)");
        Console.WriteLine("    (1,0)");
        Console.WriteLine("    (C,S)");
        Console.WriteLine("");
        Console.WriteLine("  C = cosine ( theta), S = sine ( theta ).");
        Console.WriteLine("");
        Console.WriteLine("  Test  Theta  Curvature");
        Console.WriteLine("");

        for (test = 1; test <= TEST_NUM; test++)
        {
            double theta = 2.0 * Math.PI * (test - 1)
                           / (TEST_NUM - 1);

            double theta_degrees = 360.0 * (test - 1)
                                   / (TEST_NUM - 1);

            p3[0] = Math.Cos(theta);
            p3[1] = Math.Sin(theta);

            Geometry.circle_exp2imp_2d(p1, p2, p3, ref r, ref pc);

            double curvature = r switch
            {
                > 0.0 => 1.0 / r,
                _ => 0.0
            };

            Console.WriteLine("  " + test.ToString().PadLeft(4)
                                   + "  " + theta_degrees.ToString(CultureInfo.InvariantCulture).PadLeft(5)
                                   + "  " + curvature.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    public static void test0156()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0156 tests CIRCLE_EXP2IMP_2D and CIRCLE_IMP2EXP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        double[] p1 = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        double[] p3 = new double[DIM_NUM];
        double[] pc1 = new double[DIM_NUM];
        double[] pc2 = new double[DIM_NUM];
        double r2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST0156");
        Console.WriteLine("  CIRCLE_EXP2IMP_2D converts an explicit circle");
        Console.WriteLine("  to an implicit circle.");
        Console.WriteLine("  CIRCLE_IMP2EXP_2D converts an implicit circle");
        Console.WriteLine("  to an explicit circle.");

        pc1[0] = 10.0;
        pc1[1] = 5.0;
        const double r1 = 3.0;

        Geometry.circle_imp_print_2d(r1, pc1, "  The implicit circle:");

        Geometry.circle_imp2exp_2d(r1, pc1, p1, p2, p3);

        typeMethods.r8vec_print(DIM_NUM, p1, "  P1:");
        typeMethods.r8vec_print(DIM_NUM, p2, "  P2:");
        typeMethods.r8vec_print(DIM_NUM, p3, "  P3:");

        Geometry.circle_exp2imp_2d(p1, p2, p3, ref r2, ref pc2);

        Geometry.circle_imp_print_2d(r2, pc2, "  The recovered implicit circle:");

    }

    public static void test016()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST016 tests CIRCLE_IMP_POINTS_2D and POLYGON_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        double[] pc = {5.0, -2.0};
        const double r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST016");
        Console.WriteLine("  CIRCLE_IMP_POINTS_2D gets points on a circle;");
        Console.WriteLine("  POLYGON_AREA_2D finds the area of a polygon.");

        Geometry.circle_imp_print_2d(r, pc, "  The implicit circle:");

        Console.WriteLine("");
        Console.WriteLine("  The area = " + Math.PI * r * r + "");

        int n = 8;

        double[] p = Geometry.circle_imp_points_2d(r, pc, n);

        typeMethods.r8mat_transpose_print(DIM_NUM, n, p, "  Sample results:");

        Console.WriteLine("");
        Console.WriteLine("  For any N, the sampled points define a polygon");
        Console.WriteLine("  whose area approximates the circle area.");
        Console.WriteLine("");
        Console.WriteLine("  N      Area");
        Console.WriteLine("");

        for (n = 3; n <= 24; n++)
        {
            p = Geometry.circle_imp_points_2d(r, pc, n);
            double result = Burkardt.Polygon.Geometry.polygon_area_2d(n, p);
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    public static void test0165()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0165 tests CIRCLE_IMP_POINTS_3D.
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
        const int DIM_NUM = 3;
        const int N = 12;

        double[] nc = {1.0, 1.0, 1.0};
        double[] pc = {5.0, -2.0, 1.0};
        const double r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST0165");
        Console.WriteLine("  CIRCLE_IMP_POINTS_3D gets points on a circle in 3D;");

        Geometry.circle_imp_print_3d(r, pc, nc, "  The implicit circle:");

        double[] p = Geometry.circle_imp_points_3d(r, pc, nc, N);

        typeMethods.r8mat_transpose_print(DIM_NUM, N, p, "  Points on the circle:");

    }

    public static void circles_intersect_area_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLES_INTERSECT_AREA_2D_TEST tests CIRCLES_INTERSECT_AREA_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] d_test = {1.5, 1.0, 0.5, 1.5, 1.0, 0.0};
        int i;
        const int ntest = 6;
        double[] r1_test = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        double[] r2_test = {0.5, 0.5, 0.5, 1.0, 1.0, 1.0};

        Console.WriteLine("");
        Console.WriteLine("CIRCLES_INTERSECT_AREA_2D_TEST");
        Console.WriteLine("  CIRCLES_INTERSECT_AREA_2D determines the area of the");
        Console.WriteLine("  intersection of two circes of radius R1 and R2,");
        Console.WriteLine("  with a distance D between the centers.");
        Console.WriteLine("");
        Console.WriteLine("      R1      R2       D    Area");
        Console.WriteLine("");

        for (i = 0; i < ntest; i++)
        {
            double r1 = r1_test[i];
            double r2 = r2_test[i];
            double d = d_test[i];
            double area = Geometry.circles_intersect_area_2d(r1, r2, d);

            Console.WriteLine("  " + r1.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + r2.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + d.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    public static void circles_intersect_points_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLES_INTERSECT_POINTS_2D_TEST tests CIRCLES_INTERSECT_POINTS_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;
        const int TEST_NUM = 5;

        int int_num = 0;
        double[] pint = new double[DIM_NUM * 2];
        double[] pc1 = {0.0, 0.0};
        double[] pc2_test =
        {
            5.0, 5.0,
            7.0710678, 7.0710678,
            4.0, 0.0,
            6.0, 0.0,
            0.0, 0.0
        };
        const double r1 = 5.0;
        double[] r2_test = {0.5, 5.0, 3.0, 3.0, 5.0};
        int test;

        Console.WriteLine("");
        Console.WriteLine("CIRCLES_INTERSECT_POINTS_2D_TEST");
        Console.WriteLine("  CIRCLES_IMP_INT_2D determines the intersections of");
        Console.WriteLine("  two circles in 2D.");

        Geometry.circle_imp_print_2d(r1, pc1, "  The first circle:");

        for (test = 0; test < TEST_NUM; test++)
        {
            double r2 = r2_test[test];
            double[] pc2 = pc2_test.Skip(test * DIM_NUM).ToArray();

            Geometry.circle_imp_print_2d(r2, pc2, "  The second circle:");

            Geometry.circles_intersect_points_2d(r1, pc1, r2, pc2, ref int_num, ref pint);

            switch (int_num)
            {
                case 0:
                    Console.WriteLine("");
                    Console.WriteLine("  The circles do not intersect.");
                    break;
                case 1:
                    Console.WriteLine("");
                    Console.WriteLine("  The circles intersect at one point:");
                    Console.WriteLine("");
                    Console.WriteLine("        P");
                    Console.WriteLine("");
                    Console.WriteLine("  " + pint[0 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + pint[1 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                    break;
                case 2:
                    Console.WriteLine("");
                    Console.WriteLine("  The circles intersect at two points:");
                    Console.WriteLine("");
                    Console.WriteLine("        P");
                    Console.WriteLine("");
                    Console.WriteLine("  " + pint[0 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + pint[1 + 0 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                    Console.WriteLine("  " + pint[0 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + pint[1 + 1 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                    break;
                case 3:
                    Console.WriteLine("");
                    Console.WriteLine("  The circles coincide (infinite intersection).");
                    break;
            }
        }
    }

}