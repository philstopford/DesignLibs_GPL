using System;
using System.Globalization;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace QuaternionsTest;

static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for QUATERNIONS_TEST.
        //
        //  Discussion:
        //
        //    QUATERNIONS_TEST tests the QUATERNIONS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("QUATERNIONS_TEST");
        Console.WriteLine("  Test the QUATERNIONS library.");

        q8_conjugate_test();
        q8_exponentiate_test();
        q8_inverse_test();
        q8_multiply_test();
        q8_multiply2_test();
        q8_norm_test();
        q8_normal_01_test();
        q8_transpose_print_test();

        r8_acos_test();

        r8mat_print_test();
        r8mat_print_some_test();

        r8vec_print_test();
        r8vec_uniform_01_new_test();

        rotation_axis_vector_test();
        rotation_axis2mat_test();
        rotation_axis2quat_test();

        rotation_mat_vector_test();
        rotation_mat2axis_test();
        rotation_mat2quat_test();

        rotation_quat_vector_test();
        rotation_quat2axis_test();
        rotation_quat2mat_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("QUATERNIONS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void q8_conjugate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_CONJUGATE_TEST tests Q8_CONJUGATE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_CONJUGATE_TEST");
        Console.WriteLine("  Q8_CONJUGATE conjugates a quaternion;");

        for (i = 1; i <= 5; i++)
        {
            double[] q1 = typeMethods.q8_normal_01(ref seed);
            double[] q2 = typeMethods.q8_conjugate(q1);

            Console.WriteLine("");
            typeMethods.q8_transpose_print(q1, "  q1 = q8_normal_01 ( seed ):");
            typeMethods.q8_transpose_print(q2, "  q2 = q8_conjugate ( q1 ):  ");

        }
    }

    private static void q8_exponentiate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_EXPONENTIATE_TEST tests Q8_EXPONENTIATE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_EXPONENTIATE_TEST");
        Console.WriteLine("  Q8_EXPONENTIATE exponentiates a quaternion");

        for (i = 1; i <= 5; i++)
        {
            double[] q1 = typeMethods.q8_normal_01(ref seed);
            double[] q2 = typeMethods.q8_exponentiate(q1);

            Console.WriteLine("");
            typeMethods.q8_transpose_print(q1, "  q1 = q8_normal_01 ( seed ):");
            typeMethods.q8_transpose_print(q2, "  q2 = q8_exponentiate ( q1 ):");

        }
    }

    private static void q8_inverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_INVERSE_TEST tests Q8_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_INVERSE_TEST");
        Console.WriteLine("  Q8_INVERSE inverts a quaternion");

        for (i = 1; i <= 5; i++)
        {
            double[] q1 = typeMethods.q8_normal_01(ref seed);
            double[] q2 = typeMethods.q8_inverse(q1);
            double[] q3 = typeMethods.q8_multiply(q1, q2);

            Console.WriteLine("");
            typeMethods.q8_transpose_print(q1, "  q1 = q8_normal_01 ( seed ):");
            typeMethods.q8_transpose_print(q2, "  q2 = q8_inverse ( q1 ):    ");
            typeMethods.q8_transpose_print(q3, "  q3 = q8_multiply ( q1, q2 ):    ");

        }
    }

    private static void q8_multiply_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_MULTIPLY_TEST tests Q8_MULTIPLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] q1;
        double[] q2;
        double[] q3;
        int seed;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_MULTIPLY_TEST");
        Console.WriteLine("  Q8_MULTIPLY multiplies two quaternions");

        for (i = 1; i <= 5; i++)
        {
            q1 = typeMethods.q8_normal_01(ref seed);
            q2 = typeMethods.q8_normal_01(ref seed);
            q3 = typeMethods.q8_multiply(q1, q2);

            Console.WriteLine("");
            typeMethods.q8_transpose_print(q1, "  q1 = q8_normal_01 ( seed ) :");
            typeMethods.q8_transpose_print(q2, "  q2 = q8_normal_01 ( seed ) :");
            typeMethods.q8_transpose_print(q3, "  q3 = q8_multiply ( q1, q2 ):");

        }

    }

    private static void q8_multiply2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_MULTIPLY2_TEST tests Q8_MULTIPLY2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] q1;
        double[] q2;
        double[] q3;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_MULTIPLY2_TEST");
        Console.WriteLine("  Q8_MULTIPLY2 multiplies two quaternions using a matrix");

        for (i = 1; i <= 5; i++)
        {
            q1 = typeMethods.q8_normal_01(ref seed);
            q2 = typeMethods.q8_normal_01(ref seed);
            q3 = typeMethods.q8_multiply2(q1, q2);

            Console.WriteLine("");
            typeMethods.q8_transpose_print(q1, "  q1 = q8_normal_01 ( seed )  :");
            typeMethods.q8_transpose_print(q2, "  q2 = q8_normal_01 ( seed )  :");
            typeMethods.q8_transpose_print(q3, "  q3 = q8_multiply2 ( q1, q2 ):");

        }
    }

    private static void q8_normal_01_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_NORMAL_01_TEST tests Q8_NORMAL_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_NORMAL_01_TEST");
        Console.WriteLine("  Q8_NORMAL_01 computes a normally distributed quaternion.");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            double[] q = typeMethods.q8_normal_01(ref seed);
            typeMethods.q8_transpose_print(q, "Sample quaternion:");
        }

    }

    private static void q8_norm_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_NORM_TEST tests Q8_NORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] q;
        int seed;
        double value = 0;

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_NORM_TEST");
        Console.WriteLine("  Q8_NORM computes the norm of a quaternion.");

        for (i = 1; i <= 5; i++)
        {
            Console.WriteLine("");
            q = typeMethods.q8_normal_01(ref seed);
            typeMethods.q8_transpose_print(q, "  q = q8_normal_01(seed):");
            value = typeMethods.q8_norm(q);
            Console.WriteLine("  q8_norm(q) = " + value + "");
        }
    }

    private static void q8_transpose_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Q8_TRANSPOSE_PRINT_TEST tests Q8_TRANSPOSE_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("Q8_TRANSPOSE_PRINT_TEST");
        Console.WriteLine("  Q8_TRANSPOSE_PRINT prints a quaternion 'transposed',");
        Console.WriteLine("  that is, writing it as a row vector.");

        double[] q = typeMethods.q8_normal_01(ref seed);

        typeMethods.q8_transpose_print(q, "  The quaternion:");

    }

    private static void r8_acos_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8_ACOS_TEST tests typeMethods.r8_ACOS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c;
        int test;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8_ACOS_TEST");
        Console.WriteLine("  typeMethods.r8_ACOS computes the arc-cosine of an angle.");
        Console.WriteLine("");
        Console.WriteLine("          C            typeMethods.r8_ACOS(C)        ACOS(C)");
        Console.WriteLine("");

        for (test = -1; test <= 13; test++)
        {
            c = (test - 6) / (double) 6;

            string cout = c.ToString().PadLeft(14) + "  "
                                                                               + typeMethods.r8_acos(c).ToString().PadLeft(14);

            switch (c)
            {
                case >= -1.0 and <= 1.0:
                    cout += "  " + Math.Acos(c).ToString().PadLeft(14);
                    break;
            }

            Console.WriteLine(cout);
        }
    }

    private static void r8mat_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8MAT_PRINT_TEST tests typeMethods.r8MAT_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 6;
        int N = 4;

        double[] a = new double[M * N];
        int i;
        int j;
        int m = M;
        int n = N;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8MAT_PRINT_TEST");
        Console.WriteLine("  typeMethods.r8MAT_PRINT prints an typeMethods.r8MAT.");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print(m, n, a, "  The matrix:");

    }

    private static void r8mat_print_some_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8MAT_PRINT_SOME_TEST tests typeMethods.r8MAT_PRINT_SOME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 6;
        int N = 4;

        double[] a = new double[M * N];
        int i;
        int j;
        int m = M;
        int n = N;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8MAT_PRINT_SOME_TEST");
        Console.WriteLine("  typeMethods.r8MAT_PRINT_SOME prints some of an typeMethods.r8MAT.");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print_some(m, n, a, 2, 1, 4, 2, "  Rows 2:4, Cols 1:2:");

    }

    private static void r8vec_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8VEC_PRINT_TEST tests typeMethods.r8VEC_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = {123.456, 0.000005, -1.0E+06, 3.14159265};
        int n = 4;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8VEC_PRINT_TEST");
        Console.WriteLine("  typeMethods.r8VEC_PRINT prints an typeMethods.r8VEC.");

        typeMethods.r8vec_print(n, a, "  The typeMethods.r8VEC:");
    }

    private static void r8vec_uniform_01_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    typeMethods.r8VEC_UNIFORM_01_NEW_TEST tests typeMethods.r8VEC_UNIFORM_01_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;

        int j;

        Console.WriteLine("");
        Console.WriteLine("typeMethods.r8VEC_UNIFORM_01_NEW_TEST");
        Console.WriteLine("  typeMethods.r8VEC_UNIFORM_01_NEW returns a random typeMethods.r8VEC");
        Console.WriteLine("  with entries in [ 0.0, 1.0 ]");

        int seed = 123456789;

        for (j = 1; j <= 3; j++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Input SEED = " + seed + "");
            Console.WriteLine("");

            double[] r = UniformRNG.r8vec_uniform_01_new(N, ref seed);

            typeMethods.r8vec_print(N, r, "  Random typeMethods.r8VEC:");

        }
    }

    private static void rotation_axis2mat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_AXIS2MAT_TEST tests ROTATION_AXIS2MAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] axis1 = {0.2361737, -0.8814124, -0.4090649};
        double[] axis2 = {0.0, 0.0, 2.0};
        double[] v1 = {1.0, 4.0, 10.0};
        double[] v2 = {1.0, 1.0, 1.0};

        Console.WriteLine("");
        Console.WriteLine("ROTATION_AXIS2MAT_TEST");
        Console.WriteLine("  ROTATION_AXIS2MAT converts a rotation axis to a matrix.");

        typeMethods.r8vec_print(3, v1, "  The vector V:");
        typeMethods.r8vec_print(3, axis1, "  The rotation axis:");
        double angle = 1.159804;
        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        double[] a = Rotation.rotation_axis2mat(axis1, angle);

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix A:");

        double[] w = typeMethods.r8mat_mv_new(3, 3, a, v1);

        typeMethods.r8vec_print(3, w, "  The rotated vector W = A * V:");

        //
        //  Test an axis vector that does not have unit length.
        //
        typeMethods.r8vec_print(3, v2, "  The vector V:");
        typeMethods.r8vec_print(3, axis2, "  The rotation axis:");

        angle = 90.0;
        angle = Helpers.degrees_to_radians(angle);

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        a = Rotation.rotation_axis2mat(axis2, angle);

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix A:");

        w = typeMethods.r8mat_mv_new(3, 3, a, v2);

        typeMethods.r8vec_print(3, w, "  The rotated vector W = A * V:");

    }

    private static void rotation_axis2quat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_AXIS2QUAT_TEST tests ROTATION_AXIS2QUAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle;
        double[] axis1 = {0.2361737, -0.8814124, -0.4090649};
        double[] axis2 = {0.0, 0.0, 2.0};
        double[] q;
        double[] v1 = {1.0, 4.0, 10.0};
        double[] v2 = {1.0, 1.0, 1.0};
        double[] w;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_AXIS2QUAT_TEST");
        Console.WriteLine("  ROTATION_AXIS2QUAT converts a rotation axis to a quaternion.");

        typeMethods.r8vec_print(3, v1, "  The vector V:");
        typeMethods.r8vec_print(3, axis1, "  The rotation axis:");

        angle = 1.159804;
        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        q = Rotation.rotation_axis2quat(axis1, angle);

        typeMethods.r8vec_print(4, q, "  The rotation quaternion Q:");

        w = Rotation.rotation_quat_vector(q, v1);

        typeMethods.r8vec_print(3, w, "  The rotated vector W:");

        //
        //  Another test of ROTATION_AXIS_VECTOR with an axis vector
        //  that does not have unit length.
        //
        typeMethods.r8vec_print(3, v2, "  The vector V:");
        typeMethods.r8vec_print(3, axis2, "  The rotation axis:");

        angle = 90.0;
        angle = Helpers.degrees_to_radians(angle);

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");
        q = Rotation.rotation_axis2quat(axis2, angle);

        typeMethods.r8vec_print(4, q, "  The rotation quaternion Q:");

        w = Rotation.rotation_quat_vector(q, v2);

        typeMethods.r8vec_print(3, w, "  The rotated vector W:");

    }

    private static void rotation_axis_vector_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_AXIS_VECTOR_TEST tests ROTATION_AXIS_VECTOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle;
        double[] axis1 = {0.2361737, -0.8814124, -0.4090649};
        double[] axis2 = {0.0, 0.0, 2.0};
        double[] v1 = {1.0, 4.0, 10.0};
        double[] v2 = {1.0, 1.0, 1.0};
        double[] w;

        angle = 1.159804;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_AXIS_VECTOR_TEST");
        Console.WriteLine("  ROTATION_AXIS_VECTOR applies an axis");
        Console.WriteLine("  rotation to a vector.");

        typeMethods.r8vec_print(3, v1, "  The vector:");

        typeMethods.r8vec_print(3, axis1, "  The rotation axis:");

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        w = Rotation.rotation_axis_vector(axis1, angle, v1);

        typeMethods.r8vec_print(3, w, "  The rotated vector:");

        //
        //  Another test of ROTATION_AXIS_VECTOR with an axis vector
        //  that does not have unit length.
        //
        typeMethods.r8vec_print(3, v2, "  The vector:");
        typeMethods.r8vec_print(3, axis2, "  The rotation axis:");

        angle = 90.0;
        angle = Helpers.degrees_to_radians(angle);

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        w = Rotation.rotation_axis_vector(axis2, angle, v2);

        typeMethods.r8vec_print(3, w, "  The rotated vector:");

    }

    private static void rotation_mat2axis_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_MAT2AXIS_TEST tests ROTATION_MAT2AXIS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            0.43301269, -0.5, 0.75,
            0.25, 0.86602539, 0.43301269,
            -0.86602539, 0.0, 0.5
        };
        double angle = 0;
        double[] axis;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_MAT2AXIS_TEST");
        Console.WriteLine("  ROTATION_MAT2AXIS computes a rotation axis");
        Console.WriteLine("  and angle from a rotation matrix.");
        Console.WriteLine("  ROTATION_AXIS2MAT computes a rotation matrix");
        Console.WriteLine("  from a rotation axis and angle.");

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix:");

        axis = new double[3];
        Rotation.rotation_mat2axis(a, axis, ref angle);

        typeMethods.r8vec_print(3, axis, "  The rotation axis:");

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        Rotation.rotation_axis2mat(axis, angle);

        typeMethods.r8mat_print(3, 3, a, "  The recovered rotation matrix:");

    }

    private static void rotation_mat2quat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_MAT2QUAT_TEST tests ROTATION_MAT2QUAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            0.43301269, -0.5, 0.75,
            0.25, 0.86602539, 0.43301269,
            -0.86602539, 0.0, 0.5
        };
        double[] a2;
        double[] q;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_MAT2QUAT_TEST");
        Console.WriteLine("  ROTATION_MAT2QUAT computes a quaternion");
        Console.WriteLine("  from a rotation matrix.");
        Console.WriteLine("  ROTATION_QUAT2MAT computes a rotation matrix");
        Console.WriteLine("  from a quaternion.");

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix:");

        q = Rotation.rotation_mat2quat(a);

        typeMethods.r8vec_print(4, q, "  The rotation quaternion Q:");

        a2 = Rotation.rotation_quat2mat(q);

        typeMethods.r8mat_print(3, 3, a2, "  The recovered rotation matrix:");

    }

    private static void rotation_mat_vector_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_MAT_VECTOR_TEST tests ROTATION_MAT_VECTOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a =
        {
            0.43301269, -0.5, 0.75,
            0.25, 0.86602539, 0.43301269,
            -0.86602539, 0.0, 0.5
        };
        double[] v = {1.0, 4.0, 10.0};
        double[] w;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_MAT_VECTOR_TEST");
        Console.WriteLine("  ROTATION_MAT_VECTOR applies a matrix");
        Console.WriteLine("  rotation to a vector.");

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix:");
        typeMethods.r8vec_print(3, v, "  The vector V:");

        w = Rotation.rotation_mat_vector(a, v);
        typeMethods.r8vec_print(3, w, "  The rotated vector W = A * V:");

    }

    private static void rotation_quat2axis_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_QUAT2AXIS_TEST tests ROTATION_QUAT2AXIS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle = 0;
        double[] axis;
        double[] q = {0.836516, 0.12941, -0.482963, -0.224144};
        double[] q2;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_QUAT2AXIS_TEST");
        Console.WriteLine("  ROTATION_QUAT2AXIS computes a rotation axis");
        Console.WriteLine("  and angle from a rotation quaternion.");
        Console.WriteLine("  ROTATION_AXIS2QUAT computes a rotation");
        Console.WriteLine("  quaternion from a rotation axis and angle.");

        typeMethods.r8vec_print(4, q, "  The rotation quaternion:");

        axis = new double[3];
        Rotation.rotation_quat2axis(q, axis, ref angle);

        typeMethods.r8vec_print(3, axis, "  The rotation axis:");

        Console.WriteLine("");
        Console.WriteLine("  The rotation angle is " + angle + "");

        q2 = Rotation.rotation_axis2quat(axis, angle);

        typeMethods.r8vec_print(4, q2, "  The recovered rotation quaternion:");

    }

    private static void rotation_quat2mat_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_QUAT2MAT_TEST tests ROTATION_QUAT2MAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] q = {0.836516, 0.12941, -0.482963, -0.224144};
        double[] q2;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_QUAT2MAT_TEST");
        Console.WriteLine("  ROTATION_QUAT2MAT computes a rotation axis");
        Console.WriteLine("  from a rotation quaternion.");
        Console.WriteLine("  ROTATION_MAT2QUAT computes a rotation");
        Console.WriteLine("  quaternion from a rotation matrix.");

        typeMethods.r8vec_print(4, q, "  The rotation quaternion:");

        a = Rotation.rotation_quat2mat(q);

        typeMethods.r8mat_print(3, 3, a, "  The rotation matrix:");

        q2 = Rotation.rotation_mat2quat(a);

        typeMethods.r8vec_print(4, q2, "  The recovered rotation quaternion:");
    }

    private static void rotation_quat_vector_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROTATION_QUAT_VECTOR_TEST tests ROTATION_QUAT_VECTOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //   06 August 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] q = {0.836516, 0.12941, -0.482963, -0.224144};
        double[] v = {1.0, 4.0, 10.0};
        double[] w;

        Console.WriteLine("");
        Console.WriteLine("ROTATION_QUAT_VECTOR_TEST");
        Console.WriteLine("  ROTATION_QUAT_VECTOR applies a quaternion");
        Console.WriteLine("  rotation to a vector.");

        typeMethods.r8vec_print(4, q, "  The rotation quaternion:");

        typeMethods.r8vec_print(3, v, "  The vector V:");

        w = Rotation.rotation_quat_vector(q, v);
        typeMethods.r8vec_print(3, w, "  The rotated vector:");

    }
}