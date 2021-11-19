using System;
using Burkardt.Types;

namespace PolPakTest;

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
        //    02 June 2007
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
        Console.WriteLine("   N     K    CNK");
        Console.WriteLine("");

        for (n = 0; n <= 4; n++)
        {
            for (k = 0; k <= n; k++)
            {
                cnk = typeMethods.i4_choose(n, k);

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + k.ToString().PadLeft(6) + "  "
                                  + cnk.ToString().PadLeft(6) + "");
            }
        }

    }

    public static void i4_factor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTOR_TEST tests I4_FACTOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i = 0;
        int j = 0;
        int maxfactor = 10;
        int n = 0;
        int[] n_test = { 60, 664048, 8466763 };
        int nfactor = 0;
        int nleft = 0;
        int[] factor = new int[10];
        int[] power = new int[10];

        Console.WriteLine("");
        Console.WriteLine("I4_FACTOR_TEST:");
        Console.WriteLine("  I4_FACTOR tries to factor an I4");

        for (i = 0; i < 3; i++)
        {
            n = n_test[i];
            typeMethods.i4_factor(n, maxfactor, ref nfactor, ref factor, ref power, ref nleft);
            Console.WriteLine("");
            Console.WriteLine("  Factors of N = " + n + "");
            for (j = 0; j < nfactor; j++)
            {
                Console.WriteLine("    " + factor[j] + "^" + power[j] + "");
            }

            if (nleft != 1)
            {
                Console.WriteLine("  Unresolved factor NLEFT = " + nleft + "");
            }
        }

    }

    public static void i4_factorial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL_TEST tests I4_FACTORIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fn = 0;
        int fn2 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("I4_FACTORIAL_TEST:");
        Console.WriteLine("  I4_FACTORIAL evaluates the factorial function.");
        Console.WriteLine("");
        Console.WriteLine("     X       Exact F       I4_FACTORIAL(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            typeMethods.i4_factorial_values(ref n_data, ref n, ref fn);

            if (n_data == 0)
            {
                break;
            }

            fn2 = typeMethods.i4_factorial(n);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + fn.ToString().PadLeft(12) + "  "
                              + fn2.ToString().PadLeft(12) + "");

        }

    }

    public static void i4_factorial2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FACTORIAL2_TEST tests I4_FACTORIAL2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fn = 0;
        int fn2 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("I4_FACTORIAL2_TEST:");
        Console.WriteLine("  I4_FACTORIAL2 evaluates the double factorial function.");
        Console.WriteLine("");
        Console.WriteLine("   N   Exact  I4_FACTORIAL2(N)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            typeMethods.i4_factorial2_values(ref n_data, ref n, ref fn);

            if (n_data == 0)
            {
                break;
            }

            fn2 = typeMethods.i4_factorial2(n);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + fn.ToString().PadLeft(8) + "  "
                              + fn2.ToString().PadLeft(8) + "");
        }

    }

    public static void i4_is_fibonacci_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_FIBONACCI_TEST tests I4_IS_FIBONACCI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int i4;
        int[] i4_test = { -13, 0, 1, 8, 10, 50, 55, 100, 144, 200 };
        bool l;
        int test_num = 10;

        Console.WriteLine("");
        Console.WriteLine("I4_IS_FIBONACCI_TEST");
        Console.WriteLine("  I4_IS_FIBONACCI determines if an I4 is a Fibonacci number.");
        Console.WriteLine("");
        Console.WriteLine("   I4       T/F");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            i4 = i4_test[i];
            l = typeMethods.i4_is_fibonacci(i4);

            Console.WriteLine(" " + i4.ToString().PadLeft(4)
                                  + "  " + l.ToString().PadLeft(1) + "");
        }

    }

    public static void i4_is_triangular_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_IS_TRIANGULAR_TEST tests I4_IS_TRIANGULAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        bool l;

        Console.WriteLine("");
        Console.WriteLine("I4_IS_TRIANGULAR_TEST");
        Console.WriteLine("  I4_IS_TRIANGULAR returns 0 or 1 depending on");
        Console.WriteLine("  whether I is triangular.");
        Console.WriteLine("");
        Console.WriteLine("   I  =>   0/1");
        Console.WriteLine("");

        for (i = 0; i <= 20; i++)
        {
            l = typeMethods.i4_is_triangular(i);

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + l.ToString().PadLeft(1) + "");
        }

    }

    public static void i4_partition_distinct_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_DISTINCT_COUNT_TEST tests I4_PARTITION_DISTINCT_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int c2 = 0;
        int n = 0;
        int n_data;
        const int N_MAX = 20;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_DISTINCT_COUNT_TEST:");
        Console.WriteLine("  For the number of partitions of an integer");
        Console.WriteLine("  into distinct parts,");
        Console.WriteLine("  I4_PARTITION_DISTINCT_COUNT computes any value.");
        Console.WriteLine("");
        Console.WriteLine("     N       Exact F    Q(N)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Partition.partition_distinct_count_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            if (N_MAX < n)
            {
                continue;
            }

            c2 = typeMethods.i4_partition_distinct_count(n);

            Console.WriteLine("  "
                              + n.ToString().PadLeft(10) + "  "
                              + c.ToString().PadLeft(10) + "  "
                              + c2.ToString().PadLeft(10) + "");

        }

    }

    public static void i4_to_triangle_lower_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_TRIANGLE_LOWER_TEST tests I4_TO_TRIANGLE_LOWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i = 0;
        int j = 0;
        int k = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_TRIANGLE_LOWER_TEST");
        Console.WriteLine("  I4_TO_TRIANGLE_LOWER converts a linear index to a");
        Console.WriteLine("  lower triangular one.");
        Console.WriteLine("");
        Console.WriteLine("     K  => I     J");
        Console.WriteLine("");

        for (k = 0; k <= 20; k++)
        {
            typeMethods.i4_to_triangle_lower(k, ref i, ref j);

            Console.WriteLine("  " + k.ToString().PadLeft(4)
                                   + "  " + i.ToString().PadLeft(4)
                                   + "  " + j.ToString().PadLeft(4) + "");
        }

    }

}