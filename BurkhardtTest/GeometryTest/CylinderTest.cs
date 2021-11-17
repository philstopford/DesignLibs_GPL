using System;
using System.Linq;
using Burkardt.Cylinder;
using Burkardt.Types;

namespace GeometryTest;

public static class CylinderTest
{
    public static void cylinder_point_dist_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_DIST_3D_TEST tests CYLINDER_POINT_DIST_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        double dist;
        double[] dist_test =
        {
            3.0, 0.5, 5.0, 8.0, 1.0, 0.25
        };
        double[] p;
        double[] p_test =
        {
            4.0, 0.5, 0.0,
            -0.5, -1.0, 0.0,
            4.0, 6.0, 0.0,
            0.75, -10.0, 0.0,
            0.0, 0.0, 0.0,
            0.25, 1.75, 0.0
        };
        double[] p1 = {0.0, -2.0, 0.0};
        double[] p2 = {0.0, 2.0, 0.0};
        double r = 1.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("CYLINDER_POINT_DIST_3D_TEST");
        Console.WriteLine("  CYLINDER_POINT_DIST_3D computes the distance");
        Console.WriteLine("  to a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Console.WriteLine("");
            Console.WriteLine("  P ="
                              + "  " + p[0].ToString().PadLeft(10)
                              + "  " + p[1].ToString().PadLeft(10)
                              + "  " + p[2].ToString().PadLeft(10) + "");

            dist = Geometry.cylinder_point_dist_3d(p1, p2, r, p);

            Console.WriteLine("  Distance (computed) ="
                              + "  " + dist.ToString().PadLeft(10) + "");
            Console.WriteLine("  Distance (exact)    ="
                              + "  " + dist_test[test].ToString().PadLeft(10) + "");
        }
    }

    public static void cylinder_point_dist_signed_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_DIST_SIGNED_3D_TEST tests CYLINDER_POINT_DIST_SIGNED_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        double dist;
        double[] dist_test =
        {
            3.0, -0.5, 5.0, 8.0, -1.0, -0.25
        };
        double[] p;
        double[] p_test =
        {
            4.0, 0.5, 0.0,
            -0.5, -1.0, 0.0,
            4.0, 6.0, 0.0,
            0.75, -10.0, 0.0,
            0.0, 0.0, 0.0,
            0.25, 1.75, 0.0
        };
        double[] p1 = {0.0, -2.0, 0.0};
        double[] p2 = {0.0, 2.0, 0.0};
        double r = 1.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("CYLINDER_POINT_DIST_SIGNED_3D_TEST");
        Console.WriteLine("  CYLINDER_POINT_DIST_SIGNED_3D computes the signed");
        Console.WriteLine("  distance to a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(test * DIM_NUM).ToArray();

            Console.WriteLine("");
            Console.WriteLine("  P ="
                              + "  " + p[0].ToString().PadLeft(10)
                              + "  " + p[1].ToString().PadLeft(10)
                              + "  " + p[2].ToString().PadLeft(10) + "");

            dist = Geometry.cylinder_point_dist_signed_3d(p1, p2, r, p);

            Console.WriteLine("  Signed distance (computed) ="
                              + "  " + dist.ToString().PadLeft(10) + "");
            Console.WriteLine("  Signed distance (exact)    ="
                              + "  " + dist_test[test].ToString().PadLeft(10) + "");
        }
    }

    public static void test0202()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0202 tests CYLINDER_POINT_INSIDE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        bool inside;
        bool[] inside_test =
        {
            false, true, false, false, true, true
        };
        double[] p;
        double[] p_test =
        {
            4.0, 0.5, 0.0,
            -0.5, -1.0, 0.0,
            4.0, 6.0, 0.0,
            0.75, -10.0, 0.0,
            0.0, 0.0, 0.0,
            0.25, 1.75, 0.0
        };
        double[] p1 = {0.0, -2.0, 0.0};
        double[] p2 = {0.0, 2.0, 0.0};
        double r = 1.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0202");
        Console.WriteLine("  CYLINDER_POINT_INSIDE_3D determines if a point");
        Console.WriteLine("  is inside a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(test * DIM_NUM).ToArray();

            Console.WriteLine("");
            Console.WriteLine("  P ="
                              + "  " + p[0].ToString().PadLeft(10)
                              + "  " + p[1].ToString().PadLeft(10)
                              + "  " + p[2].ToString().PadLeft(10) + "");

            inside = Geometry.cylinder_point_inside_3d(p1, p2, r, p);

            Console.WriteLine("  INSIDE (computed) ="
                              + "  " + inside.ToString().PadLeft(1) + "");
            Console.WriteLine("  INSIDE (exact)    ="
                              + "  " + inside_test[test].ToString().PadLeft(1) + "");
        }
    }

    public static void test0203()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0203 tests CYLINDER_POINT_NEAR_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 6;

        double[] p;
        double[] p_test =
        {
            4.0, 0.5, 0.0,
            -0.5, -1.0, 0.0,
            4.0, 6.0, 0.0,
            0.75, -10.0, 0.0,
            0.0, 0.0, 0.0,
            0.25, 1.75, 0.0
        };
        double[] p1 = {0.0, -2.0, 0.0};
        double[] p2 = {0.0, 2.0, 0.0};
        double[] pn;
        double[] pn_test =
        {
            1.0, 0.5, 0.0,
            -1.0, -1.0, 0.0,
            1.0, 2.0, 0.0,
            0.75, -2.0, 0.0,
            1.0, 0.0, 0.0,
            0.25, 2.0, 0.0
        };
        double r = 1.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0203");
        Console.WriteLine("  CYLINDER_POINT_NEAR_3D computes the nearest point");
        Console.WriteLine("  on a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Console.WriteLine("");
            Console.WriteLine("  P ="
                              + "  " + p[0].ToString().PadLeft(10)
                              + "  " + p[1].ToString().PadLeft(10)
                              + "  " + p[2].ToString().PadLeft(10) + "");

            pn = Geometry.cylinder_point_near_3d(p1, p2, r, p);

            Console.WriteLine("  PN (computed) ="
                              + "  " + pn[0].ToString().PadLeft(10)
                              + "  " + pn[1].ToString().PadLeft(10)
                              + "  " + pn[2].ToString().PadLeft(10) + "");
            Console.WriteLine("  PN (exact)    ="
                              + "  " + pn_test[0 + test * DIM_NUM].ToString().PadLeft(10)
                              + "  " + pn_test[1 + test * DIM_NUM].ToString().PadLeft(10)
                              + "  " + pn_test[2 + test * DIM_NUM].ToString().PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  (Note that case 5 is ambiguous.  The set of nearest");
        Console.WriteLine("  points forms a circle, any of which will do.)");

    }

    public static void test02035()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02035 tests CYLINDER_SAMPLE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 20;

        double[] p;
        double[] p1 = {0.0, -2.0, 0.0};
        double[] p2 = {0.0, 2.0, 0.0};
        double r = 1.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST02035");
        Console.WriteLine("  CYLINDER_SAMPLE_3D samples points in a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        p = Geometry.cylinder_sample_3d(p1, p2, r, N, ref seed);

        typeMethods.r8mat_transpose_print(DIM_NUM, N, p, "  Sample points:");

    }

    public static void test0204()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0204 tests CYLINDER_VOLUME_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] p1 = {1.0, 2.0, 3.0};
        double[] p2 = {5.0, 6.0, 5.0};
        double pi = 3.141592653589793;
        double r = 5.0;
        double volume;

        Console.WriteLine("");
        Console.WriteLine("TEST0204");
        Console.WriteLine("  CYLINDER_VOLUME_3D computes the volume of a cylinder.");

        Console.WriteLine("");
        Console.WriteLine("  Radius R = " + r + "");
        Console.WriteLine("  Center of bottom disk ="
                          + "  " + p1[0].ToString().PadLeft(10)
                          + "  " + p1[1].ToString().PadLeft(10)
                          + "  " + p1[2].ToString().PadLeft(10) + "");
        Console.WriteLine("  Center of top disk =   "
                          + "  " + p2[0].ToString().PadLeft(10)
                          + "  " + p2[1].ToString().PadLeft(10)
                          + "  " + p2[2].ToString().PadLeft(10) + "");

        volume = Geometry.cylinder_volume_3d(p1, p2, r);

        Console.WriteLine("");
        Console.WriteLine("  Volume (computed) = " + volume + "");
        Console.WriteLine("  Volume (exact)    = " + pi * 150.0 + "");

    }

}