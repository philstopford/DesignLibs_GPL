using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace LegendreProductPolynomialTest;

public static class r8Test
{
    public static void r8mat_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT_TEST tests R8MAT_PRINT.
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
        const int M = 6;
        const int N = 4;

        double[] a = new double[M * N];
        int j;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PRINT_TEST");
        Console.WriteLine("  R8MAT_PRINT prints an R8MAT.");

        for (j = 0; j < N; j++)
        {
            int i;
            for (i = 0; i < M; i++)
            {
                a[i + j * M] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print(M, N, a, "  The R8MAT:");
    }

    public static void r8mat_print_some_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
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
        const int M = 6;
        const int N = 4;

        double[] a = new double[M * N];
        int j;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PRINT_SOME_TEST");
        Console.WriteLine("  R8MAT_PRINT_SOME prints some of an R8MAT.");

        for (j = 0; j < N; j++)
        {
            int i;
            for (i = 0; i < M; i++)
            {
                a[i + j * M] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print_some(M, N, a, 2, 1, 4, 2, "  The R8MAT, rows 2:4, cols 1:2:");
    }

    public static void r8mat_uniform_ab_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UNIFORM_AB_NEW_TEST tests R8MAT_UNIFORM_AB_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 5;
        const int N = 4;

        const double b = 2.0E+00;
        const double c = 10.0E+00;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_UNIFORM_AB_NEW_TEST");
        Console.WriteLine("  R8MAT_UNIFORM_AB_NEW returns a random R8MAT in [A,B].");
        Console.WriteLine("");

        double[] a = UniformRNG.r8mat_uniform_ab_new(M, N, b, c, ref seed);

        typeMethods.r8mat_print(M, N, a, "  The random R8MAT:");
    }

    public static void r8vec_permute_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PERMUTE_TEST tests R8VEC_PERMUTE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 5;
        int[] p = {1, 3, 4, 0, 2};
        double[] x = {1.0, 2.0, 3.0, 4.0, 5.0};

        Console.WriteLine("");
        Console.WriteLine("R8VEC_PERMUTE_TEST");
        Console.WriteLine("  R8VEC_PERMUTE permutes an R8VEC in place.");

        typeMethods.r8vec_print(n, x, "  Original Array X[]:");

        typeMethods.i4vec_print(n, p, "  Permutation Vector P[]:");

        typeMethods.r8vec_permute(n, p, ref x);

        typeMethods.r8vec_print(n, x, "  Permuted array X[P[]]:");
    }

    public static void r8vec_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_PRINT_TEST tests R8VEC_PRINT.
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
        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST1335");
        Console.WriteLine("  R8VEC_PRINT prints an R8VEC.");

        typeMethods.r8vec_print(n, a, "  The R8VEC:");
    }

    public static void r8vec_uniform_ab_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_UNIFORM_AB_NEW_TEST tests R8VEC_UNIFORM_AB_NEW.
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

        const double a = 10.0;
        const double b = 20.0;
        int j;

        Console.WriteLine("");
        Console.WriteLine("R8VEC_UNIFORM_AB_NEW_TEST");
        Console.WriteLine("  R8VEC_UNIFORM returns a random R8VEC");
        Console.WriteLine("  with entries in a given range [ A, B ]");
        Console.WriteLine("");
        Console.WriteLine("  For this problem:");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("");

        int seed = 123456789;

        for (j = 1; j <= 3; j++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Input SEED = " + seed + "");
            Console.WriteLine("");

            double[] r = UniformRNG.r8vec_uniform_ab_new(N, a, b, ref seed);

            typeMethods.r8vec_print(N, r, "  Random R8VEC:");
        }
    }
}