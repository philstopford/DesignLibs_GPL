using System;
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
        int cnk;
        int k;
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4_CHOOSE_TEST");
        Console.WriteLine("  I4_CHOOSE evaluates C(N,K).");
        Console.WriteLine("");
        Console.WriteLine("       N       K     CNK");

        for (n = 0; n <= 4; n++)
        {
            Console.WriteLine("");
            for (k = 0; k <= n; k++)
            {
                cnk = typeMethods.i4_choose(n, k);

                Console.WriteLine("  "
                                  + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + k.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                  + cnk.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
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
        int a = -100;
        int b = 200;
        int i;
        int j;
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
            j = UniformRNG.i4_uniform_ab(a, b, ref seed);

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
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
        int[] a;
        int b;
        int c;
        int n = 12;
        int[] p;
        int seed;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PERMUTE_TEST");
        Console.WriteLine("  I4VEC_PERMUTE reorders an integer vector");
        Console.WriteLine("  according to a given permutation.");
        Console.WriteLine("  Using initial random number seed = " + seed + "");

        b = 0;
        c = n;
        a = UniformRNG.i4vec_uniform_ab_new(n, b, c, ref seed);

        typeMethods.i4vec_print(n, a, "  A, before rearrangement:");

        p = typeMethods.perm_uniform_new(n, ref seed);

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
        int n = 4;
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
        int[] a;
        int b;
        int c;
        int i;
        int[] indx;
        int n = 20;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_HEAP_INDEX_A_TEST");
        Console.WriteLine("  I4VEC_SORT_HEAP_INDEX_A creates an ascending");
        Console.WriteLine("  sort index for an I4VEC.");

        b = 0;
        c = 3 * n;
        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(n, b, c, ref seed);

        typeMethods.i4vec_print(n, a, "  Unsorted array:");

        indx = typeMethods.i4vec_sort_heap_index_a(n, a);

        typeMethods.i4vec_print(n, indx, "  Sort vector INDX:");

        Console.WriteLine("");
        Console.WriteLine("       I   INDX(I)  A(INDX(I))");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + indx[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + a[indx[i]].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
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
        int[] a;
        int hi;
        int lo;
        int n;
        int s;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SUM_TEST");
        Console.WriteLine("  I4VEC_SUM sums the entries of an I4VEC.");

        n = 5;
        lo = 0;
        hi = 10;
        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, a, "  The vector:");

        s = typeMethods.i4vec_sum(n, a);
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
        int a = -100;
        int b = 200;
        int n = 20;
        int seed = 123456789;
        int[] v;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_UNIFORM_AB_NEW_TEST");
        Console.WriteLine("  I4VEC_UNIFORM_AB_NEW computes pseudorandom values");
        Console.WriteLine("  in an interval [A,B].");

        Console.WriteLine("");
        Console.WriteLine("  The lower endpoint A = " + a + "");
        Console.WriteLine("  The upper endpoint B = " + b + "");
        Console.WriteLine("  The initial seed is " + seed + "");
        Console.WriteLine("");

        v = UniformRNG.i4vec_uniform_ab_new(n, a, b, ref seed);

        typeMethods.i4vec_print(n, v, "  The random vector:");

    }
}