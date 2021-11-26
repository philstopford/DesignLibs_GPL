﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace LegendreProductPolynomialTest;

public static class i4Test
{
    public static void i4_choose_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_CHOOSE_TEST tests I4_CHOOSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4_CHOOSE_TEST");
        Console.WriteLine("  I4_CHOOSE evaluates C(N,K).");
        Console.WriteLine("");
        Console.WriteLine("       N       K     CNK");

        for (n = 0; n <= 4; n++)
        {
            Console.WriteLine("");
            int k;
            for (k = 0; k <= n; k++)
            {
                int cnk = typeMethods.i4_choose(n, k);

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + k.ToString().PadLeft(6) + "  "
                                  + cnk.ToString().PadLeft(6) + "");
            }
        }

    }

    public static void i4_uniform_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int a = -100;
        const int b = 200;
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4_UNIFORM_AB_TEST");
        Console.WriteLine("  I4_UNIFORM_AB computes pseudorandom values");
        Console.WriteLine("  in an interval [A,B].");

        Console.WriteLine("");
        Console.WriteLine("  The lower endpoint A = " + a + "");
        Console.WriteLine("  The upper endpoint B = " + b + "");
        Console.WriteLine("  The initial seed is " + seed + "");
        Console.WriteLine("");

        for (i = 1; i <= 20; i++)
        {
            int j = UniformRNG.i4_uniform_ab(a, b, ref seed);

            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + j.ToString().PadLeft(8) + "");
        }

    }

    public static void i4vec_permute_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PERMUTE_TEST tests I4VEC_PERMUTE.
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
        const int n = 12;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PERMUTE_TEST");
        Console.WriteLine("  I4VEC_PERMUTE reorders an integer vector");
        Console.WriteLine("  according to a given permutation.");
        Console.WriteLine("  Using initial random number seed = " + seed + "");

        const int b = 0;
        int[] a = UniformRNG.i4vec_uniform_ab_new(n, b, n, ref seed);

        typeMethods.i4vec_print(n, a, "  A, before rearrangement:");

        int[] p = typeMethods.perm_uniform_new(n, ref seed);

        typeMethods.i4vec_print(n, p, "  Permutation vector P:");

        typeMethods.i4vec_permute(n, p, ref a);

        typeMethods.i4vec_print(n, a, "  A, after rearrangement:");

    }

    public static void i4vec_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PRINT_TEST tests I4VEC_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4;
        int[] v = {91, 92, 93, 94};

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PRINT_TEST");
        Console.WriteLine("  I4VEC_PRINT prints an I4VEC");

        typeMethods.i4vec_print(n, v, "  Here is the I4VEC:");

    }

    public static void i4vec_sort_heap_index_a_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SORT_HEAP_INDEX_A_TEST tests I4VEC_SORT_HEAP_INDEX_A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 20;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_HEAP_INDEX_A_TEST");
        Console.WriteLine("  I4VEC_SORT_HEAP_INDEX_A creates an ascending");
        Console.WriteLine("  sort index for an I4VEC.");

        const int b = 0;
        const int c = 3 * n;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, b, c, ref seed);

        typeMethods.i4vec_print(n, a, "  Unsorted array:");

        int[] indx = typeMethods.i4vec_sort_heap_index_a(n, a);

        typeMethods.i4vec_print(n, indx, "  Sort vector INDX:");

        Console.WriteLine("");
        Console.WriteLine("       I   INDX(I)  A(INDX(I))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(8) + "  "
                              + indx[i].ToString().PadLeft(8) + "  "
                              + a[indx[i]].ToString().PadLeft(8) + "");
        }

    }

    public static void i4vec_sum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SUM_TEST tests I4VEC_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("I4VEC_SUM_TEST");
        Console.WriteLine("  I4VEC_SUM sums the entries of an I4VEC.");

        const int n = 5;
        const int lo = 0;
        const int hi = 10;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, a, "  The vector:");

        int s = typeMethods.i4vec_sum(n, a);
        Console.WriteLine("");
        Console.WriteLine("  The vector entries sum to " + s + "");


    }

    public static void i4vec_uniform_ab_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_UNIFORM_AB_NEW_TEST tests I4VEC_UNIFORM_AB_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int a = -100;
        const int b = 200;
        const int n = 20;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_UNIFORM_AB_NEW_TEST");
        Console.WriteLine("  I4VEC_UNIFORM_AB_NEW computes pseudorandom values");
        Console.WriteLine("  in an interval [A,B].");

        Console.WriteLine("");
        Console.WriteLine("  The lower endpoint A = " + a + "");
        Console.WriteLine("  The upper endpoint B = " + b + "");
        Console.WriteLine("  The initial seed is " + seed + "");
        Console.WriteLine("");

        int[] v = UniformRNG.i4vec_uniform_ab_new(n, a, b, ref seed);

        typeMethods.i4vec_print(n, v, "  The random vector:");

    }
}