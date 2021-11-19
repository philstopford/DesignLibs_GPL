using System;
using System.Globalization;
using System.Linq;
using Burkardt;
using Burkardt.Vector;

namespace GeometryTest;

public static class VectorTest
{
    public static void vector_directions_nd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_DIRECTIONS_ND_TEST tests VECTOR_DIRECTIONS_ND;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 3;

        double[] angle = new double[DIM_NUM];
        double[] angle_degrees = new double[DIM_NUM];
        int j;
        int test;
        double[] v;
        double[] v_test =
        {
            1.0, 0.0, 0.0,
            1.0, 2.0, 3.0,
            0.0, 0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("VECTOR_DIRECTIONS_ND_TEST");
        Console.WriteLine("  VECTOR_DIRECTIONS_ND computes the angles");
        Console.WriteLine("  that a vector makes with the axes.");
        Console.WriteLine("");
        Console.WriteLine("    X       Y       Z      AX      AY      AZ " +
                          "     AX      AY      AZ   ");
        Console.WriteLine("                         (_____Radians_______)" +
                          " (_______Degrees_______)");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            v = v_test.Skip(+test * DIM_NUM).ToArray();

            Geometry.vector_directions_nd(DIM_NUM, v, ref angle);

            for (j = 0; j < DIM_NUM; j++)
            {
                angle_degrees[j] = Helpers.radians_to_degrees(angle[j]);
            }

            Console.WriteLine("  " + v[0].ToString().PadLeft(7)
                                   + "  " + v[1].ToString().PadLeft(7)
                                   + "  " + v[2].ToString().PadLeft(7)
                                   + "  " + angle[0].ToString().PadLeft(7)
                                   + "  " + angle[1].ToString().PadLeft(7)
                                   + "  " + angle[2].ToString().PadLeft(7)
                                   + "  " + angle_degrees[0].ToString().PadLeft(7)
                                   + "  " + angle_degrees[1].ToString().PadLeft(7)
                                   + "  " + angle_degrees[2].ToString().PadLeft(7) + "");
        }
    }

    public static void vector_rotate_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_ROTATE_2D_TEST tests VECTOR_ROTATE_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 3;

        double angle;
        double[] a_test = {30.0, -45.0, 270.0};
        int test;
        double[] v;
        double[] v_test =
        {
            1.0, 0.0,
            0.0, 2.0,
            1.0, 1.0
        };
        double[] w = new double[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("VECTOR_ROTATE_2D_TEST");
        Console.WriteLine("  VECTOR_ROTATE_2D rotates a vector through");
        Console.WriteLine("  a given angle around the origin.");
        Console.WriteLine("");
        Console.WriteLine("    X1      Y1   Angle      X2      Y2");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            v = v_test.Skip(+test * DIM_NUM).ToArray();

            angle = Helpers.degrees_to_radians(a_test[test]);

            Geometry.vector_rotate_2d(v, angle, ref w);

            Console.WriteLine("  " + v[0].ToString().PadLeft(7)
                                   + "  " + v[1].ToString().PadLeft(7)
                                   + "  " + a_test[test].ToString().PadLeft(7)
                                   + "  " + w[0].ToString().PadLeft(7)
                                   + "  " + w[1].ToString().PadLeft(7) + "");
        }

    }

    public static void vector_rotate_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_ROTATE_3D_TEST tests VECTOR_ROTATE_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 5;

        double[] axis = {1.0, 1.0, 1.0};
        double angle;
        double[] a_test = {30.0, -45.0, 90.0, 270.0, 30.0};
        int test;
        double[] v1 = new double[DIM_NUM];
        double[] v1_test =
        {
            1.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 3.0,
            1.0, 1.0, 1.0,
            1.0, 1.0, -2.0
        };
        double[] v2 = new double[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("VECTOR_ROTATE_3D_TEST");
        Console.WriteLine("  VECTOR_ROTATE_3D rotates a vector through");
        Console.WriteLine("  a given angle around the origin.");
        Console.WriteLine("");
        Console.WriteLine("  Rotations will be about the following axis:");
        Console.WriteLine("");
        Console.WriteLine("  " + axis[0].ToString().PadLeft(8)
                               + "  " + axis[1].ToString().PadLeft(8)
                               + "  " + axis[2].ToString().PadLeft(8) + "");
        Console.WriteLine("");
        Console.WriteLine("              V1             Angle             V2");
        Console.WriteLine("    ----------------------  ------  ----------------------");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            v1 = v1_test.Skip(+test * DIM_NUM).ToArray();

            angle = Helpers.degrees_to_radians(a_test[test]);

            Geometry.vector_rotate_3d(v1, axis, angle, ref v2);

            Console.WriteLine("  " + v1[0].ToString("0.####").PadLeft(7)
                                   + "  " + v1[1].ToString("0.####").PadLeft(7)
                                   + "  " + v1[2].ToString("0.####").PadLeft(7)
                                   + "  " + angle.ToString("0.####").PadLeft(7)
                                   + "  " + v2[0].ToString("0.####").PadLeft(7)
                                   + "  " + v2[1].ToString("0.####").PadLeft(7)
                                   + "  " + v2[2].ToString("0.####").PadLeft(7) + "");
        }

        //
        //  Test using an axis that is not of unit length!
        //
        axis[0] = 0.0;
        axis[1] = 0.0;
        axis[2] = 2.0;
        Console.WriteLine("");
        Console.WriteLine("  Rotations will be about the following axis:");
        Console.WriteLine("");
        Console.WriteLine("  " + axis[0].ToString().PadLeft(8)
                               + "  " + axis[1].ToString().PadLeft(8)
                               + "  " + axis[2].ToString().PadLeft(8) + "");
        Console.WriteLine("");
        Console.WriteLine("              V1             Angle             V2");
        Console.WriteLine("    ----------------------  ------  ----------------------");

        v1[0] = 1.0;
        v1[1] = 1.0;
        v1[2] = 1.0;

        angle = 90.0;
        angle = Helpers.degrees_to_radians(angle);

        Geometry.vector_rotate_3d(v1, axis, angle, ref v2);

        Console.WriteLine("  " + v1[0].ToString("0.####").PadLeft(7)
                               + "  " + v1[1].ToString("0.####").PadLeft(7)
                               + "  " + v1[2].ToString("0.####").PadLeft(7)
                               + "  " + 90.0.ToString("0.####").PadLeft(7)
                               + "  " + v2[0].ToString("0.####").PadLeft(7)
                               + "  " + v2[1].ToString("0.####").PadLeft(7)
                               + "  " + v2[2].ToString("0.####").PadLeft(7) + "");

    }

    public static void vector_rotate_base_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_ROTATE_BASE_2D_TEST tests VECTOR_ROTATE_BASE_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double angle;
        double[] a_test = {30.0, -45.0, 270.0, 20.0};
        double[] p1;
        double[] p2 = new double[DIM_NUM];
        double[] pb = {10.0, 5.0};
        double[] p_test =
        {
            11.0, 5.0,
            10.0, 7.0,
            11.0, 6.0,
            10.0, 5.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("VECTOR_ROTATE_BASE_2D_TEST");
        Console.WriteLine("  VECTOR_ROTATE_BASE_2D rotates a vector (X1,Y1)");
        Console.WriteLine("  through an angle around a base point (XB,YB).");
        Console.WriteLine("");
        Console.WriteLine("        P1              PB       Angle          P2");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = p_test.Skip(+test * DIM_NUM).ToArray();

            angle = Helpers.degrees_to_radians(a_test[test]);

            Geometry.vector_rotate_base_2d(p1, pb, angle, ref p2);

            Console.WriteLine("  " + p1[0].ToString().PadLeft(7)
                                   + "  " + p1[1].ToString().PadLeft(7)
                                   + "  " + pb[0].ToString().PadLeft(7)
                                   + "  " + pb[1].ToString().PadLeft(7)
                                   + "  " + a_test[test].ToString().PadLeft(7)
                                   + "  " + p2[0].ToString().PadLeft(7)
                                   + "  " + p2[1].ToString().PadLeft(7) + "");
        }

    }

    public static void vector_separation_nd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_SEPARATION_ND_TEST tests VECTOR_SEPARATION_ND;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 5;

        int test1;
        int test2;
        double theta;
        double theta_deg;
        double[] v1;
        double[] v2;
        double[] v_test =
        {
            1.0, 0.0, 0.0,
            1.0, 2.0, 3.0,
            0.0, 0.0, 1.0,
            -3.0, 2.0, -1.0,
            -2.0, -4.0, -6.0
        };

        Console.WriteLine("");
        Console.WriteLine("VECTOR_SEPARATION_ND_TEST");
        Console.WriteLine("  VECTOR_SEPARATION_ND computes the separation angle");
        Console.WriteLine("  between two vectors.");
        Console.WriteLine("");
        Console.WriteLine("    -----Vector 1-----      -----Vector 2-----  ");
        Console.WriteLine("   Radians    Degrees");
        Console.WriteLine("");

        for (test1 = 0; test1 < TEST_NUM; test1++)
        {
            v1 = v_test.Skip(+test1 * DIM_NUM).ToArray();

            for (test2 = test1 + 1; test2 < TEST_NUM; test2++)
            {
                v2 = v_test.Skip(+test2 * DIM_NUM).ToArray();

                theta = Geometry.vector_separation_nd(DIM_NUM, v1, v2);

                theta_deg = Helpers.radians_to_degrees(theta);

                Console.WriteLine("  " + v1[0].ToString().PadLeft(7)
                                       + "  " + v1[1].ToString().PadLeft(7)
                                       + "  " + v1[2].ToString().PadLeft(7)
                                       + "  " + v2[0].ToString().PadLeft(7)
                                       + "  " + v2[1].ToString().PadLeft(7)
                                       + "  " + v2[2].ToString().PadLeft(7)
                                       + "  " + theta.ToString().PadLeft(7)
                                       + "  " + theta_deg.ToString().PadLeft(7) + "");
            }
        }
    }

}