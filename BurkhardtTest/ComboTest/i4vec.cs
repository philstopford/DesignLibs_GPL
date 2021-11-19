using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ComboTest;

internal partial class Program
{
    private static void i4vec_backtrack_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_BACKTRACK_TEST tests I4VEC_BACKTRACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int found_num;
        int i;
        int indx;
        int k;
        int n = 8;
        int maxstack = 100;
        int[] ncan = new int[8];
        int nstack;
        int[] stacks = new int[100];
        int t;
        int total;
        int[] w =  {
                15, 22, 14, 26, 32, 9, 16, 8
            }
            ;
        int[] x = new int[8];

        Console.WriteLine("");
        Console.WriteLine("I4VEC_BACKTRACK_TEST");
        Console.WriteLine("  I4VEC_BACKTRACK uses backtracking, seeking an I4VEC X of");
        Console.WriteLine("  N values which satisfies some condition.");

        Console.WriteLine("");
        Console.WriteLine("  In this demonstration, we have 8 integers W(I).");
        Console.WriteLine("  We seek all subsets that sum to 53.");
        Console.WriteLine("  X(I) is 0 or 1 if the entry is skipped or used.");
        Console.WriteLine("");

        t = 53;

        for (i = 0; i < n; i++)
        {
            x[i] = 0;
        }

        indx = 0;
        k = 0;
        nstack = 0;
        for (i = 0; i < n; i++)
        {
            ncan[i] = 0;
        }

        found_num = 0;

        for (;;)
        {
            typeMethods.i4vec_backtrack(n, maxstack, stacks, ref x, ref indx, ref k, ref nstack, ref ncan);

            if (indx == 1)
            {
                found_num += 1;
                string cout = "  " + found_num.ToString().PadLeft(2) + "   ";

                total = typeMethods.i4vec_dot_product(n, w, x);
                cout += "  " + total.ToString().PadLeft(3) + ":  ";

                for (i = 0; i < n; i++)
                {
                    switch (x[i])
                    {
                        case 1:
                            cout += "  " + w[i].ToString().PadLeft(2);
                            break;
                    }
                }

                Console.WriteLine(cout);
            }
            //
            //  Given that we've chose X(1:K-1), what are our choices for X(K)?
            //
            //     if T < TOTAL, 
            //       no choices
            //     if T = TOTAL, 
            //       X(K) = 0
            //     if T > TOTAL and K < N, 
            //       X(k) = 0
            //       If ( W(K)+TOTAL <= T ) X(K) = 1
            //     If T > TOTAL and K = N,
            //       If ( W(K) + TOTAL) = T ) X(K) = 1
            //
            else if (indx == 2)
            {
                total = typeMethods.i4vec_dot_product(k - 1, w, x);

                if (t < total)
                {
                    ncan[k - 1] = 0;
                }
                else if (t == total)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0;
                    nstack += 1;
                }
                else if (total < t && k < n)
                {
                    ncan[k - 1] += 1;
                    stacks[nstack] = 0;
                    nstack += 1;

                    if (total + w[k - 1] <= t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1;
                        nstack += 1;
                    }
                }
                else if (total < t && k == n)
                {
                    if (total + w[k - 1] == t)
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 1;
                        nstack += 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("  Done!");
                break;
            }
        }
    }

    private static void i4vec_dot_product_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_DOT_PRODUCT_TEST tests I4VEC_DOT_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a;
        int[] b;
        int d;
        int hi;
        int lo;
        int n = 5;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_DOT_PRODUCT_TEST");
        Console.WriteLine("  I4VEC_DOT_PRODUCT computes the dot product of two I4VECs.");

        lo = 0;
        hi = 10;
        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, a, "  The vector A:");

        b = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, b, "  The vector B:");

        d = typeMethods.i4vec_dot_product(n, a, b);
        Console.WriteLine("");
        Console.WriteLine("  The dot product is " + d + "");
    }

    private static void i4vec_part1_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PART1_NEW_TEST tests I4VEC_PART1_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int npart;
        int[] x;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART1_NEW_TEST:");
        Console.WriteLine("  I4VEC_PART1_NEW partitions an integer N into NPART parts.");

        n = 17;
        npart = 5;

        Console.WriteLine("");
        Console.WriteLine("  Partition N = " + n + " into NPART = " + npart + " parts:");
        Console.WriteLine("");

        x = typeMethods.i4vec_part1_new(n, npart);

        typeMethods.i4vec_print(npart, x, "  The partition:");
    }

    private static void i4vec_part2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PART2_TEST tests I4VEC_PART2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int npart;
        int[] x;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART2_TEST:");
        Console.WriteLine("  I4VEC_PART2 partitions an integer N into NPART parts.");

        n = 17;
        npart = 5;
        x = new int[npart];

        Console.WriteLine("");
        Console.WriteLine("  Partition N = " + n + " into NPART = " + npart + " parts:");
        Console.WriteLine("");

        typeMethods.i4vec_part2(n, npart, ref x);

        typeMethods.i4vec_print(npart, x, "  The partition:");
    }

    private static void i4vec_part2_new_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PART2_NEW_TEST tests I4VEC_PART2_NEW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int npart;
        int[] x;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART2_NEW_TEST:");
        Console.WriteLine("  I4VEC_PART2_NEW partitions an integer N into NPART parts.");

        n = 17;
        npart = 5;

        Console.WriteLine("");
        Console.WriteLine("  Partition N = " + n + " into NPART = " + npart + " parts:");
        Console.WriteLine("");

        x = typeMethods.i4vec_part2_new(n, npart);

        typeMethods.i4vec_print(npart, x, "  The partition:");
    }

    private static void i4vec_search_binary_a_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SEARCH_BINARY_A_TEST tests I4VEC_SEARCH_BINARY_A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =  {
                0, 1, 1, 2, 3, 4, 5, 6, 7, 8
            }
            ;
        int b;
        int index;
        int n = 10;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SEARCH_BINARY_A_TEST");
        Console.WriteLine("  I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC.");

        typeMethods.i4vec_print(n, a, "  Ascending sorted I4VEC:");

        b = 5;

        Console.WriteLine("");
        Console.WriteLine("  Now search for an instance of the value " + b + "");

        index = typeMethods.i4vec_search_binary_a(n, a, b);

        Console.WriteLine("");
        switch (index)
        {
            case 0:
                Console.WriteLine("  The value does not occur.");
                break;
            default:
                Console.WriteLine("  The value occurs at index = " + index + "");
                break;
        }
    }

    private static void i4vec_search_binary_d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SEARCH_BINARY_D_TEST tests I4VEC_SEARCH_BINARY_D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int[] a =  {
                8, 7, 6, 5, 4, 3, 2, 1, 1, 0
            }
            ;
        int b;
        int index;
        int n = N;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SEARCH_BINARY_D_TEST");
        Console.WriteLine("  I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC.");

        typeMethods.i4vec_print(n, a, "  Descending sorted I4VEC:");

        b = 5;

        Console.WriteLine("");
        Console.WriteLine("  Now search for an instance of the value " + b + "");

        index = typeMethods.i4vec_search_binary_d(n, a, b);

        Console.WriteLine("");
        switch (index)
        {
            case 0:
                Console.WriteLine("  The value does not occur.");
                break;
            default:
                Console.WriteLine("  The value occurs at index = " + index + "");
                break;
        }
    }

    private static void i4vec_sort_insert_a_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SORT_INSERT_A_TEST tests I4VEC_SORT_INSERT_A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =  {
                6, 7, 1, 0, 4, 3, 2, 1, 5, 8
            }
            ;
        int n = 10;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_INSERT_A_TEST");
        Console.WriteLine("  I4VEC_SORT_INSERT_A ascending sorts an I4VEC.");

        typeMethods.i4vec_print(n, a, "  Before ascending sort:");

        typeMethods.i4vec_sort_insert_a(n, ref a);

        typeMethods.i4vec_print(n, a, "  After ascending sort:");
    }

    private static void i4vec_sort_insert_d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SORT_INSERT_D_TEST tests I4VEC_SORT_INSERT_D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =  {
                6, 7, 1, 0, 4, 3, 2, 1, 5, 8
            }
            ;
        int n = 10;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_INSERT_D_TEST");
        Console.WriteLine("  I4VEC_SORT_INSERT_D descending sorts an I4VEC.");

        typeMethods.i4vec_print(n, a, "  Before descending sort:");

        typeMethods.i4vec_sort_insert_d(n, ref a);

        typeMethods.i4vec_print(n, a, "  After descending sort:");
    }

    private static void i4vec_uniform_ab_new_test()

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