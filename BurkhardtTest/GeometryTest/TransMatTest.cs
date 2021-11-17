using System;
using Burkardt.Geometry;
using Burkardt.Types;

namespace GeometryTest;

public static class TransMatTest
{
    public static void test204()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST204 tests TransMat.tmat_INIT, TransMat.tmat_ROT_AXIS, TransMat.tmat_ROT_VECTOR, TransMat.tmat_SCALE, TransMat.tmat_SHEAR, TransMat.tmat_TRANS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = new double[4 * 4];
        double angle;
        char axis1;
        string axis2;
        double[] b = new double[4 * 4];
        int i;
        double s;
        double[] v = new double[3];

        Console.WriteLine("");
        Console.WriteLine("TEST204");
        Console.WriteLine("  TMAT geometric transformation matrix routines:");
        Console.WriteLine("  TransMat.tmat_INIT initializes,");
        Console.WriteLine("  TransMat.tmat_ROT_AXIS for rotation about an axis,");
        Console.WriteLine("  TransMat.tmat_ROT_VECTOR for rotation about a vector,");
        Console.WriteLine("  TransMat.tmat_SCALE for scaling,");
        Console.WriteLine("  TransMat.tmat_SHEAR for shear,");
        Console.WriteLine("  TransMat.tmat_TRANS for translation.");
        //
        //  Initialization.
        //
        TransMat.tmat_init(a);

        Console.WriteLine("");
        Console.WriteLine("  Initial transformation matrix:");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + a[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + a[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + a[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + a[i + 3 * 4].ToString().PadLeft(10) + "");
        }

        //
        //  Rotation about an axis.
        //
        angle = 30.0;
        axis1 = 'x';
        TransMat.tmat_rot_axis(a, ref b, angle, axis1);

        Console.WriteLine("");
        Console.WriteLine("  Transformation matrix for");
        Console.WriteLine("  rotation about " + axis1 + " by " + angle + " degrees.");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 3 * 4].ToString().PadLeft(10) + "");
        }

        //
        //  Rotation about a vector.
        //
        angle = 30.0;
        v[0] = 1.0;
        v[1] = 2.0;
        v[2] = 3.0;
        TransMat.tmat_rot_vector(a, ref b, angle, v);

        Console.WriteLine("");
        Console.WriteLine("  Transformation matrix for");
        Console.WriteLine("  rotation about " + v[0] + "  " + v[1] + "  " + v[2]
                          + " by " + angle + " degrees.");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 3 * 4].ToString().PadLeft(10) + "");
        }

        //
        //  Scaling.
        //
        v[0] = 2.0;
        v[1] = 0.5;
        v[2] = 10.0;
        TransMat.tmat_scale(a, ref b, v);

        Console.WriteLine("");
        Console.WriteLine("  Transformation matrix for");
        Console.WriteLine("  scaling by " + v[0] + "  " + v[1] + "  " + v[2] + "");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 3 * 4].ToString().PadLeft(10) + "");
        }

        //
        //  Shear.
        //
        axis2 = "xy";
        s = 0.5;
        TransMat.tmat_shear(a, ref b, axis2, s);

        Console.WriteLine("");
        Console.WriteLine("  Transformation matrix for");
        Console.WriteLine("  " + axis2 + " shear coefficient of " + s + "");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 3 * 4].ToString().PadLeft(10) + "");
        }

        //
        //  Translation.
        //
        v[0] = 1.0;
        v[1] = 2.0;
        v[2] = 3.0;
        TransMat.tmat_trans(a, ref b, v);

        Console.WriteLine("");
        Console.WriteLine("  Transformation matrix for");
        Console.WriteLine("  translation by " + v[0] + "  " + v[1] + "  " + v[2] + "");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 1 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 2 * 4].ToString().PadLeft(10)
                                   + "  " + b[i + 3 * 4].ToString().PadLeft(10) + "");
        }
    }

    public static void test205()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST205 tests TransMat.tmat_MXP2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 4;

        double[] a = new double[4 * 4];
        double angle;
        char axis1;
        string axis2;
        double[] b = new double[4 * 4];
        double[] point = new double[DIM_NUM * N];
        double[] point2 = new double[DIM_NUM * N];
        double s;
        double[] v = new double[3];

        Console.WriteLine("");
        Console.WriteLine("TEST205");
        Console.WriteLine("  TransMat.tmat_MXP2 applies a geometric transformation");
        Console.WriteLine("  matrix to a set of points.");
        //
        //  Initialization.
        //
        point[0 + 0 * 3] = 1.0;
        point[1 + 0 * 3] = 0.0;
        point[2 + 0 * 3] = 0.0;

        point[0 + 1 * 3] = 0.0;
        point[1 + 1 * 3] = 1.0;
        point[2 + 1 * 3] = 0.0;

        point[0 + 2 * 3] = 0.0;
        point[1 + 2 * 3] = 0.0;
        point[2 + 2 * 3] = 1.0;

        point[0 + 3 * 3] = 1.0;
        point[1 + 3 * 3] = 1.0;
        point[2 + 3 * 3] = 1.0;

        typeMethods.r8mat_transpose_print(DIM_NUM, N, point, "  Points:");
        //
        //  Initialization of transformation matrix.
        //
        TransMat.tmat_init(a);
        //
        //  Rotation about an axis.
        //
        angle = 30.0;
        axis1 = 'x';
        TransMat.tmat_rot_axis(a, ref b, angle, axis1);

        TransMat.tmat_mxp2(b, point, ref point2, N);

        Console.WriteLine("");
        Console.WriteLine("  Rotation about " + axis1 + " by " + angle + " degrees.");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, point2, "  Transformed points:");
        //
        //  Rotation about a vector.
        //
        angle = 30.0;
        v[0] = 1.0;
        v[1] = 2.0;
        v[2] = 3.0;
        TransMat.tmat_rot_vector(a, ref b, angle, v);

        TransMat.tmat_mxp2(b, point, ref point2, N);

        Console.WriteLine("");
        Console.WriteLine("  Rotation about " + v[0]
                                              + "  " + v[1]
                                              + "  " + v[2]
                                              + " by " + angle + " degrees.");
        typeMethods.r8mat_transpose_print(DIM_NUM, N, point2, "  Transformed points:");
        //
        //  Scaling.
        //
        v[0] = 2.0;
        v[1] = 0.5;
        v[2] = 10.0;
        TransMat.tmat_scale(a, ref b, v);

        TransMat.tmat_mxp2(b, point, ref point2, N);

        Console.WriteLine("");
        Console.WriteLine("  Scaling by " + v[0]
                                          + "  " + v[1]
                                          + "  " + v[2] + "");
        typeMethods.r8mat_transpose_print(DIM_NUM, N, point2, "  Transformed points:");
        //
        //  Shear.
        //
        axis2 = "xy";
        s = 0.5;
        TransMat.tmat_shear(a, ref b, axis2, s);

        TransMat.tmat_mxp2(b, point, ref point2, N);

        Console.WriteLine("");
        Console.WriteLine("  " + axis2 + " shear coefficient of " + s + ":");
        typeMethods.r8mat_transpose_print(DIM_NUM, N, point2, "  Transformed points:");
        //
        //  Translation.
        //
        v[0] = 1.0;
        v[1] = 2.0;
        v[2] = 3.0;
        TransMat.tmat_trans(a, ref b, v);

        TransMat.tmat_mxp2(b, point, ref point2, N);

        Console.WriteLine("");
        Console.WriteLine("  Translation by " + v[0]
                                              + "  " + v[1]
                                              + "  " + v[2] + "");
        typeMethods.r8mat_transpose_print(DIM_NUM, N, point2, "  Transformed points:");

    }

}