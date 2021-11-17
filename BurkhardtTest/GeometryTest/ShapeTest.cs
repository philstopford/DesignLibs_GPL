using System;
using System.Linq;
using Burkardt.Geometry;
using Burkardt.Types;

namespace GeometryTest;

public static class ShapeTest
{
    public static void test196()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST196 tests SHAPE_POINT_DIST_2D.
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
        int DIM_NUM = 2;
        int SIDE_NUM = 4;
        int TEST_NUM = 9;

        double dist;
        int i;
        double[] p;
        double[] p1 = {5.0, 0.0};
        double[] pc = {3.0, 0.0};
        double[] p_test =
        {
            3.0, 0.0,
            5.0, 0.0,
            4.0, 0.0,
            10.0, 0.0,
            8.0, 5.0,
            6.0, 6.0,
            1.0, 2.0,
            2.5, -0.5,
            4.0, -1.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST196");
        Console.WriteLine("  For a shape in 2D,");
        Console.WriteLine("  SHAPE_POINT_DIST_2D computes the distance");
        Console.WriteLine("  to a point;");
        Console.WriteLine("");
        Console.WriteLine("  Number of sides:  " + SIDE_NUM + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Center of square:");

        typeMethods.r8vec_print(DIM_NUM, p1, "  Square vertex #1");

        Console.WriteLine("");
        Console.WriteLine("  TEST       X            Y            DIST");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist = Shape.shape_point_dist_2d(pc, p1, SIDE_NUM, p);

            string cout = "  " + test.ToString().PadLeft(6);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout + "  " + dist.ToString().PadLeft(12) + "");
        }

    }

    public static void test197()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST197 tests SHAPE_POINT_DIST_2D.
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
        int DIM_NUM = 2;
        int SIDE_NUM = 6;
        int TEST_NUM = 8;

        double dist;
        double[] p;
        double[] p1 = {5.0, 0.0};
        double[] pc = {3.0, 0.0};
        double[] p_test =
        {
            3.0, 0.0,
            5.0, 0.0,
            4.0, 0.0,
            10.0, 0.0,
            4.0, 1.7320508,
            5.0, 3.4641016,
            3.0, 1.7320508,
            3.0, 0.86602539
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST197");
        Console.WriteLine("  For a shape in 2D,");
        Console.WriteLine("  SHAPE_POINT_DIST_2D computes the distance");
        Console.WriteLine("  to a point;");
        Console.WriteLine("");
        Console.WriteLine("  Number of sides: " + SIDE_NUM + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Center of hexagon:");

        typeMethods.r8vec_print(DIM_NUM, p1, "  Hexagon vertex #1");

        Console.WriteLine("");
        Console.WriteLine("  TEST         X    Y     DIST");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist = Shape.shape_point_dist_2d(pc, p1, SIDE_NUM, p);

            Console.WriteLine("  " + test.ToString().PadLeft(6)
                                   + "  " + p[0].ToString().PadLeft(10)
                                   + "  " + p[1].ToString().PadLeft(10)
                                   + "  " + dist.ToString().PadLeft(10) + "");
        }
    }

    public static void test198()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST198 tests SHAPE_POINT_NEAR_2D.
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
        int DIM_NUM = 2;
        int SIDE_NUM = 6;
        int TEST_NUM = 8;

        double dist = 0;
        double[] p;
        double[] p_test =
        {
            3.0, 0.0,
            5.0, 0.0,
            4.0, 0.0,
            10.0, 0.0,
            4.0, 1.7320508,
            5.0, 3.4641016,
            3.0, 1.7320508,
            3.0, 0.86602539
        };
        double[] p1 = {5.0, 0.0};
        double[] pc = {3.0, 0.0};
        double[] pn = new double[DIM_NUM];
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST198");
        Console.WriteLine("  For a shape in 2D,");
        Console.WriteLine("  SHAPE_POINT_NEAR_2D computes the nearest");
        Console.WriteLine("  point to a point;");
        Console.WriteLine("");
        Console.WriteLine("  Number of sides: " + SIDE_NUM + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Hexagon center:");

        typeMethods.r8vec_print(DIM_NUM, p1, "  Hexagon vertex #1");

        Console.WriteLine("");
        Console.WriteLine("  TEST       X            Y              PN     Dist");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Shape.shape_point_near_2d(pc, p1, SIDE_NUM, p, ref pn, ref dist);

            Console.WriteLine("  " + test.ToString().PadLeft(6)
                                   + "  " + p[0].ToString().PadLeft(10)
                                   + "  " + p[1].ToString().PadLeft(10)
                                   + "  " + pn[0].ToString().PadLeft(10)
                                   + "  " + pn[1].ToString().PadLeft(10)
                                   + "  " + dist.ToString().PadLeft(10) + "");
        }

    }

    public static void test199()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST199 tests SHAPE_RAY_INT_2D.
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
        int DIM_NUM = 2;
        int SIDE_NUM = 6;
        int TEST_NUM = 4;

        double[] p1 = {5.0, 0.0};
        double[] pa;
        double[] pa_test =
        {
            3.0, 0.0,
            3.0, 0.0,
            3.0, -1.0,
            3.0, -1.0
        };
        double[] pb;
        double[] pb_test =
        {
            4.0, 0.0,
            3.0, 1.0,
            3.0, 1.0,
            7.0, 5.0
        };
        double[] pc = {3.0, 0.0};
        double[] pint = new double[DIM_NUM];
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST199");
        Console.WriteLine("  For a shape in 2D,");
        Console.WriteLine("  SHAPE_RAY_INT_2D computes the intersection of");
        Console.WriteLine("  a shape and a ray whose origin is within");
        Console.WriteLine("  the shape.");
        Console.WriteLine("");
        Console.WriteLine("  Number of sides = " + SIDE_NUM + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Hexagon center:");

        Console.WriteLine("");
        Console.WriteLine("  Hexagon vertex #1:");
        Console.WriteLine("");
        Console.WriteLine("  " + p1[0].ToString().PadLeft(10)
                               + "  " + p1[1].ToString().PadLeft(10) + "");

        Console.WriteLine("");
        Console.WriteLine("  TEST       XA          YA          XB"
                          + "          YB          XI          YI");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            pa = pa_test.Skip(+test * DIM_NUM).ToArray();
            pb = pb_test.Skip(+test * DIM_NUM).ToArray();

            Shape.shape_ray_int_2d(pc, p1, SIDE_NUM, pa, pb, ref pint);

            Console.WriteLine("  " + test.ToString().PadLeft(6)
                                   + "  " + pa[0].ToString().PadLeft(10)
                                   + "  " + pa[1].ToString().PadLeft(10)
                                   + "  " + pb[0].ToString().PadLeft(10)
                                   + "  " + pb[1].ToString().PadLeft(10)
                                   + "  " + pint[0].ToString().PadLeft(10)
                                   + "  " + pint[1].ToString().PadLeft(10) + "");
        }

    }

}