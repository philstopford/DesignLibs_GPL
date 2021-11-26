using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ComboTest;

internal static partial class Program
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
        int i;
        const int n = 8;
        const int maxstack = 100;
        int[] ncan = new int[8];
        int[] stacks = new int[100];
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

        const int t = 53;

        for (i = 0; i < n; i++)
        {
            x[i] = 0;
        }

        int indx = 0;
        int k = 0;
        int nstack = 0;
        for (i = 0; i < n; i++)
        {
            ncan[i] = 0;
        }

        int found_num = 0;

        for (;;)
        {
            typeMethods.i4vec_backtrack(n, maxstack, stacks, ref x, ref indx, ref k, ref nstack, ref ncan);

            int total;
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

                switch (total)
                {
                    case > t:
                        ncan[k - 1] = 0;
                        break;
                    case t:
                        ncan[k - 1] += 1;
                        stacks[nstack] = 0;
                        nstack += 1;
                        break;
                    case < t when k < n:
                    {
                        ncan[k - 1] += 1;
                        stacks[nstack] = 0;
                        nstack += 1;

                        if (total + w[k - 1] > t)
                        {
                            continue;
                        }

                        ncan[k - 1] += 1;
                        stacks[nstack] = 1;
                        nstack += 1;
                        break;
                    }
                    case < t when k == n:
                    {
                        if (total + w[k - 1] != t)
                        {
                            continue;
                        }

                        ncan[k - 1] += 1;
                        stacks[nstack] = 1;
                        nstack += 1;
                        break;
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
        const int n = 5;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_DOT_PRODUCT_TEST");
        Console.WriteLine("  I4VEC_DOT_PRODUCT computes the dot product of two I4VECs.");

        const int lo = 0;
        const int hi = 10;
        int seed = 123456789;

        int[] a = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, a, "  The vector A:");

        int[] b = UniformRNG.i4vec_uniform_ab_new(n, lo, hi, ref seed);
        typeMethods.i4vec_print(n, b, "  The vector B:");

        int d = typeMethods.i4vec_dot_product(n, a, b);
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
        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART1_NEW_TEST:");
        Console.WriteLine("  I4VEC_PART1_NEW partitions an integer N into NPART parts.");

        const int n = 17;
        const int npart = 5;

        Console.WriteLine("");
        Console.WriteLine("  Partition N = " + n + " into NPART = " + npart + " parts:");
        Console.WriteLine("");

        int[] x = typeMethods.i4vec_part1_new(n, npart);

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
        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART2_TEST:");
        Console.WriteLine("  I4VEC_PART2 partitions an integer N into NPART parts.");

        const int n = 17;
        const int npart = 5;
        int[] x = new int[npart];

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
        Console.WriteLine("");
        Console.WriteLine("I4VEC_PART2_NEW_TEST:");
        Console.WriteLine("  I4VEC_PART2_NEW partitions an integer N into NPART parts.");

        const int n = 17;
        const int npart = 5;

        Console.WriteLine("");
        Console.WriteLine("  Partition N = " + n + " into NPART = " + npart + " parts:");
        Console.WriteLine("");

        int[] x = typeMethods.i4vec_part2_new(n, npart);

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
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SEARCH_BINARY_A_TEST");
        Console.WriteLine("  I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC.");

        typeMethods.i4vec_print(n, a, "  Ascending sorted I4VEC:");

        const int b = 5;

        Console.WriteLine("");
        Console.WriteLine("  Now search for an instance of the value " + b + "");

        int index = typeMethods.i4vec_search_binary_a(n, a, b);

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
        const int N = 10;

        int[] a =  {
                8, 7, 6, 5, 4, 3, 2, 1, 1, 0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SEARCH_BINARY_D_TEST");
        Console.WriteLine("  I4VEC_SEARCH_BINARY_D searches a descending sorted I4VEC.");

        typeMethods.i4vec_print(N, a, "  Descending sorted I4VEC:");

        const int b = 5;

        Console.WriteLine("");
        Console.WriteLine("  Now search for an instance of the value " + b + "");

        int index = typeMethods.i4vec_search_binary_d(N, a, b);

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
        const int n = 10;

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
        const int n = 10;

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