using System;
using System.Globalization;
using Burkardt.Ellipse;
using Burkardt.Types;

namespace GeometryTest;

public static class EllipseTest
{
    public static void ellipse_area1_test()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA1_TEST tests ELLIPSE_AREA1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = {5.0, 1.0, 1.0, 2.0};

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA1_TEST");
            
        Console.WriteLine("  ELLIPSE_AREA1 computes the area of an ellipse.");

        double r = 10.0;

        double area = Geometry.ellipse_area1(a, r);

        Console.WriteLine("");
        Console.WriteLine("  R = " + r + "");
        typeMethods.r8mat_print(2, 2, a, "  Matrix A in ellipse definition x*A*x=r^2");
        Console.WriteLine("  Area = " + area + "");

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA1_TEST");
        Console.WriteLine("  Normal end of execution.");

    }

    public static void ellipse_area2_test()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA2_TEST tests ELLIPSE_AREA2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA2_TEST");
            
        Console.WriteLine("  ELLIPSE_AREA2 computes the area of an ellipse.");

        const double a = 5.0;
        const double b = 2.0;
        const double c = 2.0;
        const double d = 10.0;

        double area = Geometry.ellipse_area2(a, b, c, d);

        Console.WriteLine("");
        Console.WriteLine("  Ellipse: " + a
                                        + " * x^2 + " + b
                                        + " * xy + " + c
                                        + " * y^2 = " + d + "");
        Console.WriteLine("  Area = " + area + "");

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA2_TEST");
        Console.WriteLine("  Normal end of execution.");
    }

    public static void ellipse_area3_test()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA3_TEST tests ELLIPSE_AREA3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA3_TEST");
            
        Console.WriteLine("  ELLIPSE_AREA3 computes the area of an ellipse.");

        const double r1 = 10.0;
        const double r2 = 10.0 / 3.0;

        double area = Geometry.ellipse_area3(r1, r2);

        Console.WriteLine("");
        Console.WriteLine("  Ellipse:  (x/" + r1 + ")^2 + (y/" + r2 + ")^2 = 1");
        Console.WriteLine("  Area = " + area + "");

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_AREA3_TEST");
        Console.WriteLine("  Normal end of execution.");

    }

    public static void test025()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST025 tests ELLIPSE_POINT_DIST_2D.
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

        int i;
        const int n = 10;
        double[] p = new double[DIM_NUM];
        const double r1 = 3.0;
        const double r2 = 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST025:");
        Console.WriteLine("  ELLIPSE_POINT_DIST_2D is given a point P, and");
        Console.WriteLine("  finds the distance to an ellipse in 2D.");
        Console.WriteLine("");
        Console.WriteLine("  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1");
        Console.WriteLine("");
        Console.WriteLine("  R1 = " + r1 + "");
        Console.WriteLine("  R2 = " + r2 + "");
        Console.WriteLine("");
        Console.WriteLine("           P            DIST");
        Console.WriteLine("");

        for (i = -3; i <= n + 3; i++)
        {
            p[0] = ((n - i) * 0.0
                    + i * 4.0)
                   / n;

            p[1] = ((n - i) * 3.0
                    + i * 0.0)
                   / n;

            double dist = Geometry.ellipse_point_dist_2d(r1, r2, p);

            Console.WriteLine("  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + dist.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    public static void ellipse_point_near_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINT_NEAR_2D_TEST tests ELLIPSE_POINT_NEAR_2D.
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

        int i;
        const int n = 10;
        double[] p = new double[DIM_NUM];
        const double r1 = 3.0;
        const double r2 = 2.0;

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_POINT_NEAR_2D_TEST:");
        Console.WriteLine("  ELLIPSE_POINT_NEAR_2D is given a point P, and");
        Console.WriteLine("  finds the nearest point on an ellipse in 2D.");
        Console.WriteLine("");
        Console.WriteLine("  The ellipse is (X/R1)^2 + (Y/R2)^2 = 1");
        Console.WriteLine("");
        Console.WriteLine("  R1 = " + r1 + "");
        Console.WriteLine("  R2 = " + r2 + "");
        Console.WriteLine("");
        Console.WriteLine("           P            PN");
        Console.WriteLine("");

        for (i = -3; i <= n + 3; i++)
        {
            p[0] = ((n - i) * 0.0
                    + i * 4.0)
                   / n;

            p[1] = ((n - i) * 3.0
                    + i * 0.0)
                   / n;

            double[] pn = Geometry.ellipse_point_near_2d(r1, r2, p);

            Console.WriteLine("  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + pn[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + pn[1].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }

    }

    public static void test026()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST026 tests ELLIPSE_POINTS_2D, POLYGON_AREA_2D.
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
        const int N_MAX = 24;

        double[] pc = {5.0, -2.0};
        const double psi = 3.141592653589793 / 6.0;
        const double r1 = 3.0;
        const double r2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST026");
        Console.WriteLine("  ELLIPSE_POINTS_2D returns points on an ellipse;");
        Console.WriteLine("  POLYGON_AREA_2D finds the area of a polygon.");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Ellipse center:");

        Console.WriteLine("");
        Console.WriteLine("  radii R1 = " + r1 + " R2 = " + r2 + "");
        Console.WriteLine("  and angle PSI = " + psi + "");

        double area = Geometry.ellipse_area3(r1, r2);

        Console.WriteLine("  and area = " + area + "");

        int n = 16;
        double[] v = new double[DIM_NUM * n];

        Geometry.ellipse_points_2d(pc, r1, r2, psi, n, ref v);

        typeMethods.r8mat_transpose_print(DIM_NUM, n, v, "  Sample points:");

        Console.WriteLine("");
        Console.WriteLine("  For any N, the sampled points define a polygon");
        Console.WriteLine("  whose area approximates the ellipse area.");
        Console.WriteLine("");
        Console.WriteLine("       N         Area");
        Console.WriteLine("");

        for (n = 3; n <= N_MAX; n++)
        {
            v = new double[DIM_NUM * n];
            Geometry.ellipse_points_2d(pc, r1, r2, psi, n, ref v);
            double result = Burkardt.Polygon.Geometry.polygon_area_2d(n, v);
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void ellipse_points_arc_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_POINTS_ARC_2D_TEST tests ELLIPSE_POINTS_ARC_2D.
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
        int N = 13;

        double[] p = new double[DIM_NUM * N];
        double[] pc = {5.0, -2.0};
        double psi = 3.141592653589793 / 6.0;
        double r1 = 3.0;
        double r2 = 1.0;
        double theta1 = 3.141592653589793 / 2.0;
        double theta2 = 2.0 * 3.141592653589793;

        Console.WriteLine("");
        Console.WriteLine("ELLIPSE_POINTS_ARC_2D_TEST");
        Console.WriteLine("  ELLIPSE_POINTS_ARC_2D returns points on an");
        Console.WriteLine("  elliptical arc.");
        Console.WriteLine("");
        Console.WriteLine("  The ellipse has center " + pc[0] + "  " + pc[1] + "");
        Console.WriteLine("  radii R1 = " + r1 + " R2 = " + r2 + "");
        Console.WriteLine("  and angle PSI = " + psi + "");
        Console.WriteLine("");
        Console.WriteLine("  The arc extends from THETA1 = " + theta1 + "");
        Console.WriteLine("  to THETA2 = " + theta2 + "");

        Geometry.ellipse_points_arc_2d(pc, r1, r2, psi, theta1, theta2, N, ref p);

        typeMethods.r8mat_transpose_print(DIM_NUM, N, p, "  Sample points:");

    }

    public static void test202 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST202 tests SUPER_ELLIPSE_POINTS_2D;
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
        const int DIM_NUM = 2;
        const int N = 24;

        double[] pc = { 5.0, -2.0 };
        double[] p = new double[DIM_NUM*N];

        const double r1 = 3.0;
        const double r2 = 1.0;
        const double expo = 1.5;
        const double psi = Math.PI / 6.0;

        Console.WriteLine("");
        Console.WriteLine("TEST202");
        Console.WriteLine("  SUPER_ELLIPSE_POINTS_2D returns points on a super ellipse;");

        typeMethods.r8vec_print ( DIM_NUM, pc, "  Superellipse center:" );

        Console.WriteLine("");
        Console.WriteLine("  radii R1 = " + r1 + " R2 = " + r2 + "");
        Console.WriteLine("  exponent EXPO = " + expo + "");
        Console.WriteLine("  and angle PSI = " + psi + "");

        Geometry.super_ellipse_points_2d ( pc, r1, r2, expo, psi, N, ref p );

        typeMethods.r8mat_transpose_print ( DIM_NUM, N, p, "  Sample points:" );
    }

}