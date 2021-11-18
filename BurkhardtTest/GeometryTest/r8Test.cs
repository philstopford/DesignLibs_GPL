using System;
using System.Globalization;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class r8Test
{
    public static void test0234()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0234 tests R8MAT_SOLVE_2D.
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
        double[] a;
        double[] b = null;
        double det = 0;
        int i;
        int n = 2;
        int seed;
        int test;
        int test_num = 5;
        double[] x;
        double[] x2;

        Console.WriteLine("");
        Console.WriteLine("TEST0234");
        Console.WriteLine("  R8MAT_SOLVE_2D solves 2D linear systems.");

        seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            a = UniformRNG.r8mat_uniform_01_new(n, n, ref seed);
            x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            typeMethods.r8mat_mv(n, n, a, x, ref b);

            x2 = typeMethods.r8mat_solve_2d(a, b, ref det);

            Console.WriteLine("");
            Console.WriteLine("  Solution / Computed:");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + x2[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

        }
    }

    public static void r8_acos_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ACOS_TEST tests R8_ACOS;
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
        int TEST_NUM = 9;

        double temp1;
        double temp2;
        int test;
        double x;
        double[] x_test =
        {
            5.0, 1.2, 1.0, 0.9, 0.5, 0.0, -0.9, -1.0, -1.01
        };

        Console.WriteLine("");
        Console.WriteLine("R8_ACOS_TEST");
        Console.WriteLine("  R8_ACOS computes an angle with a given cosine;");
        Console.WriteLine("");
        Console.WriteLine("  X  R8_ACOS(X) (Degrees)");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            x = x_test[test];

            temp1 = typeMethods.r8_acos(x);
            temp2 = Helpers.radians_to_degrees(temp1);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void r8_asin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ASIN_TEST tests R8_ASIN;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 9;

        double temp1;
        double temp2;
        int test;
        double x;
        double[] x_test =
        {
            5.0, 1.2, 1.0, 0.9, 0.5, 0.0, -0.9, -1.0, -1.01
        };

        Console.WriteLine("");
        Console.WriteLine("R8_ASIN_TEST");
        Console.WriteLine("  R8_ASIN computes an angle with a given sine;");
        Console.WriteLine("");
        Console.WriteLine("  X  R8_ASIN(X)  (Degrees)");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            x = x_test[test];

            temp1 = typeMethods.r8_asin(x);
            temp2 = Helpers.radians_to_degrees(temp1);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void r8_atan_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_ATAN_TEST tests R8_ATAN;
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
        int TEST_NUM = 8;

        double temp1;
        double temp2;
        double temp3;
        int test;
        double x;
        double[] x_test = {1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 0.0};
        double y;
        double[] y_test = {0.0, 1.0, 2.0, 0.0, -1.0, -1.0, -1.0, -1.0};

        Console.WriteLine("");
        Console.WriteLine("R8_ATAN_TEST");
        Console.WriteLine("  R8_ATAN computes an angle with a given tangent.");
        Console.WriteLine("");
        Console.WriteLine("  X, Y, ATAN(Y/X), ATAN2(Y,X), R8_ATAN(Y,X)");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            x = x_test[test];
            y = y_test[test];

            if (x != 0.0)
            {
                temp1 = Math.Atan(y / x);
            }
            else
            {
                temp1 = typeMethods.r8_huge();
            }

            temp2 = Math.Atan2(y, x);
            temp3 = typeMethods.r8_atan(y, x);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp3.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Repeat, but display answers in degrees.");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            x = x_test[test];
            y = y_test[test];

            if (x != 0.0)
            {
                temp1 = Helpers.radians_to_degrees(Math.Atan(y / x));
            }
            else
            {
                temp1 = typeMethods.r8_huge();
            }

            temp2 = Helpers.radians_to_degrees(Math.Atan2(y, x));
            temp3 = Helpers.radians_to_degrees(typeMethods.r8_atan(y, x));

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp2.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + temp3.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void test0243()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0243 tests R8VEC_ANY_NORMAL.
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
        int DIM_NUM = 10;
        int TEST_NUM = 5;

        int seed;
        int test;
        double[] v1;
        double v1_norm;
        double v1v2_dot;
        double[] v2;
        double v2_norm;

        Console.WriteLine("");
        Console.WriteLine("TEST0243");
        Console.WriteLine("  R8VEC_ANY_NORMAL computes a vector V2 that is normal");
        Console.WriteLine("  to a given vector V1.");
        Console.WriteLine("");
        Console.WriteLine("    Test    ||V1||      ||V2||        V1.V2");
        Console.WriteLine("");

        seed = 123456789;

        for (test = 0; test < TEST_NUM; test++)
        {
            v1 = UniformRNG.r8vec_uniform_01_new(DIM_NUM, ref seed);
            v1_norm = typeMethods.r8vec_norm(DIM_NUM, v1);
            v2 = typeMethods.r8vec_any_normal(DIM_NUM, v1);
            v2_norm = typeMethods.r8vec_norm(DIM_NUM, v2);
            v1v2_dot = typeMethods.r8vec_dot_product(DIM_NUM, v1, v2);
            Console.WriteLine("  " + test.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + v1_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v2_norm.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v1v2_dot.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void test0245()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0245 tests R8VEC_NORMAL_01.
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
        const int N_MAX = 1000;

        int i;
        int j;
        int n;
        int seed = 123456789;
        double[] x;
        double x_max;
        double x_mean;
        double x_min;
        double x_var;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST0245");
        Console.WriteLine("  R8VEC_NORMAL_01 computes a vector of normally");
        Console.WriteLine("  distributed random numbers.");
        Console.WriteLine("  Using initial random number seed = " + seed + "");
        //
        //  Test 1:
        //  Simply call 5 times for 1 value, and print.
        //
        Console.WriteLine("");
        Console.WriteLine("  Test #1: Call 5 times, 1 value each time.");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            n = 1;
            x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);
            for (j = 0; j < n; j++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }
        }

        //
        //  Test 2:
        //  Restore the random number seed, and repeat.
        //
        Console.WriteLine("");
        Console.WriteLine("  Test #2: Restore the random number seed.");
        Console.WriteLine("  call 5 times, 1 value each time.");
        Console.WriteLine("  The results should be identical.");
        Console.WriteLine("");

        n = -1;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            n = 1;
            x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);
            for (j = 0; j < n; j++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            }
        }

        //
        //  Test 3:
        //  Restore the random number seed, compute all 5 values at once.
        //
        Console.WriteLine("");
        Console.WriteLine("  Test #3: Restore the random number seed.");
        Console.WriteLine("  call 1 time for 5 values.");
        Console.WriteLine("  The results should be identical.");
        Console.WriteLine("");

        n = -1;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        seed = 123456789;

        n = 5;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        i = 0;
        for (j = 0; j < n; j++)
        {
            i += 1;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Test 4:
        //  Restore the random number seed, compute all 5 values at once.
        //
        Console.WriteLine("");
        Console.WriteLine("  Test #4: Restore the random number seed.");
        Console.WriteLine("  call for 2, 1, and 2 values.");
        Console.WriteLine("  The results should be identical.");
        Console.WriteLine("");

        n = -1;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        seed = 123456789;

        n = 2;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);
        i = 0;

        for (j = 0; j < n; j++)
        {
            i += 1;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        n = 1;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            i += 1;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        n = 2;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            i += 1;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        //
        //  Test 5:
        //  Determine the minimum, maximum, mean and variance.
        //
        n = N_MAX;
        x = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        x_min = typeMethods.r8vec_min(n, x);
        x_max = typeMethods.r8vec_max(n, x);
        x_mean = typeMethods.r8vec_mean(n, x);
        x_var = typeMethods.r8vec_variance(n, x);

        Console.WriteLine("");
        Console.WriteLine("  Test #5:");
        Console.WriteLine("  Number of samples was " + n + "");
        Console.WriteLine("  Minimum value was     " + x_min + "");
        Console.WriteLine("  Maximum value was     " + x_max + "");
        Console.WriteLine("  Average value was     " + x_mean + "");
        Console.WriteLine("  Variance was          " + x_var + "");
        Console.WriteLine("  Expected average      " + 0.0 + "");
        Console.WriteLine("  Expected variance     " + 1.0 + "");

    }

    public static void test1745()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1745 tests R8MAT_SOLVE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;
        int RHS_NUM = 2;

        double[] a =
        {
            1.0, 4.0, 7.0,
            2.0, 5.0, 8.0,
            3.0, 6.0, 0.0,
            14.0, 32.0, 23.0,
            7.0, 16.0, 7.0
        };
        int i;
        int info;
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST1745");
        Console.WriteLine("  R8MAT_SOLVE solves linear systems.");
        //
        //  Print out the matrix to be inverted.
        //
        typeMethods.r8mat_print(N, N + RHS_NUM, a, "  The linear system:");
        //
        //  Solve the systems.
        //
        info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  The input matrix was singular.");
            Console.WriteLine("  The solutions could not be computed.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("  The computed solutions:");
        Console.WriteLine("");
        for (i = 0; i < N; i++)
        {
            string cout = "";
            for (j = N; j < N + RHS_NUM; j++)
            {
                cout += "  " + a[i + j * N].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test1746()

        //****************************************************************************80********
        //
        //  Purpose:
        //
        //    TEST1746 tests R8MAT_INVERSE_3D.
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
        double[] a =
        {
            3.0, 2.0, 0.0,
            2.0, 2.0, 1.0,
            1.0, 1.0, 1.0
        };
        double[] b;
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST1746");
        Console.WriteLine("  R8MAT_INVERSE_3D inverts a 3 by 3 matrix.");

        Console.WriteLine("");
        Console.WriteLine("  Matrix A:");
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            Console.WriteLine("  " + a[i + 0 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + a[i + 1 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + a[i + 2 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        b = typeMethods.r8mat_inverse_3d(a);

        Console.WriteLine("");
        Console.WriteLine("  Inverse matrix B:");
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            Console.WriteLine("  " + b[i + 0 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b[i + 1 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b[i + 2 * 3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }
}