using System;
using Burkardt;
using Burkardt.Values;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class i4Test
{
    public static void i4_bclr_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BCLR_TEST tests I4_BCLR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i4;
        int[] i4_test =  {
                101, -31
            }
            ;
        int j;
        int pos;
        int test;
        int test_num = 2;

        Console.WriteLine("");
        Console.WriteLine("I4_BCLR_TEST");
        Console.WriteLine("  I4_BCLR sets a given bit to 0.");

        for (test = 0; test < test_num; test++)
        {
            i4 = i4_test[test];

            Console.WriteLine("");
            Console.WriteLine("  Working on I4 = " + i4 + "");
            Console.WriteLine("");
            Console.WriteLine("       Pos     I4_BCLR(I4,POS)");
            Console.WriteLine("");

            for (pos = 0; pos < 32; pos++)
            {
                j = typeMethods.i4_bclr(i4, pos);

                Console.WriteLine("  " + pos.ToString().PadLeft(8)
                                       + "  " + j.ToString().PadLeft(12) + "");
            }
        }
    }

    public static void i4_bset_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BSET_TEST tests I4_BSET.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i4;
        int[] i4_test =  {
                101, -31
            }
            ;
        int j;
        int pos;
        int test;
        int test_num = 2;

        Console.WriteLine("");
        Console.WriteLine("I4_BSET_TEST");
        Console.WriteLine("  I4_BSET sets a given bit to 1.");

        for (test = 0; test < test_num; test++)
        {
            i4 = i4_test[test];

            Console.WriteLine("");
            Console.WriteLine("  Working on I4 = " + i4 + "");
            Console.WriteLine("");
            Console.WriteLine("       Pos     I4_BSET(I4,POS)");
            Console.WriteLine("");

            for (pos = 0; pos < 32; pos++)
            {
                j = typeMethods.i4_bset(i4, pos);

                Console.WriteLine("  " + pos.ToString().PadLeft(8)
                                       + "  " + j.ToString().PadLeft(12) + "");
            }
        }
    }

    public static void i4_btest_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_BTEST_TEST tests I4_BTEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i4;
        int[] i4_test =  {
                101, -31
            }
            ;
        bool j;
        int pos;
        int test;

        Console.WriteLine("");
        Console.WriteLine("I4_BTEST_TEST");
        Console.WriteLine("  I4_BTEST reports whether a given bit is 0 or 1.");

        for (test = 0; test < 2; test++)
        {
            i4 = i4_test[test];

            Console.WriteLine("");
            Console.WriteLine("  Analyze the integer I4 = " + i4 + "");
            Console.WriteLine("");
            Console.WriteLine("       Pos     I4_BTEST(I4,POS)");
            Console.WriteLine("");

            for (pos = 0; pos <= 31; pos++)
            {
                j = typeMethods.i4_btest(i4, pos);

                Console.WriteLine("  " + pos.ToString().PadLeft(8)
                                       + "  " + j + "");
            }
        }
    }

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
        //    05 September 2005
        //
    {
        int FACTOR_MAX = 10;

        int[] factor = new int[FACTOR_MAX];
        int factor_num = 0;
        int i;
        int n;
        int nleft = 0;
        int[] power = new int[FACTOR_MAX];

        Console.WriteLine("");
        Console.WriteLine("I4_FACTOR_TEST");
        Console.WriteLine("  I4_FACTOR factors an integer,");

        n = 2 * 2 * 17 * 37;

        Console.WriteLine("");
        Console.WriteLine("  The integer is " + n + "");

        typeMethods.i4_factor(n, FACTOR_MAX, ref factor_num, ref factor, ref power, ref nleft);

        Console.WriteLine("");
        Console.WriteLine("  Prime representation:");
        Console.WriteLine("");
        Console.WriteLine("  I  FACTOR(I)  POWER(I)");
        Console.WriteLine("");

        if (Math.Abs(nleft) != 1)
        {
            Console.WriteLine("  " + 0.ToString().PadLeft(6)
                                   + "  " + nleft.ToString().PadLeft(6)
                                   + "  " + "(Unfactored portion)");
        }

        for (i = 0; i < factor_num; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString().PadLeft(6)
                                   + "  " + factor[i].ToString().PadLeft(6)
                                   + "  " + power[i].ToString().PadLeft(6) + "");
        }
    }

    public static void i4_fall_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_FALL_TEST tests I4_FALL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int f1 = 0;
        int f2 = 0;
        int m = 0;
        int n = 0;
        int n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_FALL_TEST");
        Console.WriteLine("  I4_FALL evaluates the falling factorial function.");
        Console.WriteLine("");
        Console.WriteLine("         M         N     Exact  I4_Fall(M,N)");

        n_data = 0;

        while (true)
        {
            typeMethods.i4_fall_values(ref n_data, ref m, ref n, ref f1);

            if (n_data == 0)
            {
                break;
            }

            f2 = typeMethods.i4_fall(m, n);

            Console.WriteLine("  " + m.ToString().PadLeft(8)
                                   + "  " + n.ToString().PadLeft(8)
                                   + "  " + f1.ToString().PadLeft(8)
                                   + "  " + f2.ToString().PadLeft(8) + "");
        }
    }

    public static void i4_gcd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_GCD_TEST tests I4_GCD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4_GCD_TEST");
        Console.WriteLine("  I4_GCD computes the greatest common divisor");
        Console.WriteLine("  of two integers.");

        Console.WriteLine("");
        Console.WriteLine("     I     J    I4_GCD(I,J)");
        Console.WriteLine("");

        seed = 123456789;

        for (k = 1; k <= 15; k++)
        {
            i = UniformRNG.i4_uniform_ab(-5, 15, ref seed);
            j = UniformRNG.i4_uniform_ab(1, 15, ref seed);

            Console.WriteLine(i.ToString().PadLeft(4) + "  "
                                                      + j.ToString().PadLeft(4) + "  "
                                                      + typeMethods.i4_gcd(i, j).ToString().PadLeft(4) + "");
        }
    }

    public static void i4_huge_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_HUGE_TEST tests I4_HUGE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("I4_HUGE_TEST");
        Console.WriteLine("  I4_HUGE returns a huge integer.");
        Console.WriteLine("");
        Console.WriteLine("  I4_HUGE() = " + typeMethods.i4_huge() + "");
    }

    public static void i4_log_10_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_LOG_10_TEST tests I4_LOG_10.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 21;

        int i;
        int[] x =  {
                0, 1, 2, 3, 9, 10, 11, 99, 100, 101, 999, 1000, 1001,
                -1, -2, -3, -9, -10, -11, -99, -101
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("I4_LOG_10_TEST");
        Console.WriteLine("  I4_LOG_10: whole part of log base 10,");
        Console.WriteLine("");
        Console.WriteLine("     X I4_LOG_10");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(x[i].ToString().PadLeft(6) + "  "
                                                         + Math.Log10(x[i]).ToString().PadLeft(6) + "");

        }
    }

    public static void i4_modp_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MODP_TEST tests I4_MODP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 4;

        int[] ndivid =  {
                50, -50, 50, -50
            }
            ;
        int nmult;
        int nrem;
        int[] number =  {
                107, 107, -107, -107
            }
            ;
        int test;

        Console.WriteLine("");
        Console.WriteLine("I4_MODP_TEST");
        Console.WriteLine("  I4_MODP factors a number");
        Console.WriteLine("  into a multiple and a remainder.");
        Console.WriteLine("");
        Console.WriteLine("    Number   Divisor  Multiple Remainder");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            nrem = typeMethods.i4_modp(number[test], ndivid[test]);
            nmult = number[test] / ndivid[test];

            Console.WriteLine("  " + number[test].ToString().PadLeft(10)
                                   + "  " + ndivid[test].ToString().PadLeft(10)
                                   + "  " + nmult.ToString().PadLeft(10)
                                   + "  " + nrem.ToString().PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Repeat using C++ percent operator:");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            nrem = number[test] % ndivid[test];
            nmult = number[test] / ndivid[test];

            Console.WriteLine("  " + number[test].ToString().PadLeft(10)
                                   + "  " + ndivid[test].ToString().PadLeft(10)
                                   + "  " + nmult.ToString().PadLeft(10)
                                   + "  " + nrem.ToString().PadLeft(10) + "");
        }
    }

    public static void i4_moebius_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_MOEBIUS_TEST tests I4_MOEBIUS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c1 = 0;
        int c2;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("I4_MOEBIUS_TEST");
        Console.WriteLine("  I4_MOEBIUS evaluates the Moebius function.");
        Console.WriteLine("");
        Console.WriteLine("         N     Exact  I4_Moebius(N)");

        n_data = 0;

        while (true)
        {
            Moebius.moebius_values(ref n_data, ref n, ref c1);

            if (n_data == 0)
            {
                break;
            }

            c2 = typeMethods.i4_moebius(n);

            Console.WriteLine("  " + n.ToString().PadLeft(8)
                                   + "  " + c1.ToString().PadLeft(8)
                                   + "  " + c2.ToString().PadLeft(8) + "");
        }
    }

    public static void i4_partition_conj_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_CONJ_TEST tests I4_PARTITION_CONJ.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 14;
        int NPART1 = 4;

        int[] a1 =  {
                2, 5, 1, 4
            }
            ;
        int[] a2 = new int[N];
        int[] mult1 =  {
                1, 1, 3, 1
            }
            ;
        int[] mult2 = new int[N];
        int npart2 = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_CONJ_TEST");
        Console.WriteLine("  I4_PARTITION_CONJ conjugates an integer partition.");
        Console.WriteLine("");
        Console.WriteLine("  Original partition:");
        Console.WriteLine("");

        typeMethods.i4_partition_print(N, NPART1, a1, mult1);

        typeMethods.i4_partition_conj(N, a1, mult1, NPART1, ref a2, ref mult2, ref npart2);

        Console.WriteLine("");
        Console.WriteLine("  Conjugate partition:");
        Console.WriteLine("");

        typeMethods.i4_partition_print(N, npart2, a2, mult2);
    }

    public static void i4_partition_count_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_COUNT_TEST tests I4_PARTITION_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 20;

        int n = 0;
        int n_data;
        int p = 0;
        int[] p2 = new int[N_MAX + 1];

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_COUNT_TEST");
        Console.WriteLine("  I4_PARTITION_COUNT counts partitions of an integer.");

        n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("   N     Exact     Count");
        Console.WriteLine("");

        for (;;)
        {
            typeMethods.i4_partition_count_values(ref n_data, ref n, ref p);

            if (n_data == 0)
            {
                break;
            }

            string cout = n.ToString().PadLeft(4) + "  "
                                                  + p.ToString().PadLeft(10) + "  ";

            if (n <= N_MAX)
            {
                typeMethods.i4_partition_count(n, ref p2);
                cout += p2[n].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);

        }

    }

    public static void i4_partition_count2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_COUNT2_TEST tests I4_PARTITION_COUNT2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 20;

        int n = 0;
        int n_data;
        int p = 0;
        int[] p2;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_COUNT2_TEST");
        Console.WriteLine("  I4_PARTITION_COUNT2 counts partitions of an integer.");

        n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("   N     Exact     Count");
        Console.WriteLine("");

        for (;;)
        {
            typeMethods.i4_partition_count_values(ref n_data, ref n, ref p);

            if (n_data == 0)
            {
                break;
            }

            string cout = "  "
                          + n.ToString().PadLeft(4) + "  "
                          + p.ToString().PadLeft(10) + "  ";

            if (n <= N_MAX)
            {
                p2 = typeMethods.i4_partition_count2(n);
                cout += "  "
                        + "      "
                        + p2[n].ToString().PadLeft(10);
            }

            Console.WriteLine(cout);
        }

    }

    public static void i4_partition_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_NEXT_TEST tests I4_PARTITION_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 7;

        int[] a = new int[N];
        bool done = false;
        int[] mult = new int [N];
        int npart = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_NEXT_TEST");
        Console.WriteLine("  I4_PARTITION_NEXT generates partitions of an integer.");
        Console.WriteLine("  Here N = " + N + "");
        Console.WriteLine("");

        done = true;

        for (;;)
        {
            typeMethods.i4_partition_next(ref done, a, mult, N, ref npart);

            if (done)
            {
                break;
            }

            typeMethods.i4_partition_print(N, npart, a, mult);

        }
    }

    public static void i4_partition_next2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_NEXT2_TEST tests I4_PARTITION_NEXT2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 7;

        int[] a = new int[N];
        bool more = false;
        int[] mult = new int[N];
        int npart = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_NEXT2_TEST");
        Console.WriteLine("  I4_PARTITION_NEXT2 produces partitions of an integer.");
        Console.WriteLine("");

        more = false;

        for (;;)
        {
            typeMethods.i4_partition_next2(N, ref a, ref mult, ref npart, ref more);

            typeMethods.i4_partition_print(N, npart, a, mult);

            if (!more)
            {
                break;
            }

        }
    }

    public static void i4_partition_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_PRINT_TEST tests I4_PARTITION_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = {
                2, 5, 1, 4
            }
            ;
        int[] mult =  {
                1, 1, 3, 1
            }
            ;
        int n;
        int npart;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_PRINT_TEST");
        Console.WriteLine("  I4_PARTITION_PRINT prints an integer partition.");
        Console.WriteLine("");

        n = 14;
        npart = 4;
        typeMethods.i4_partition_print(n, npart, a, mult);
    }

    public static void i4_partition_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITION_RANDOM_TEST tests I4_PARTITION_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 8;

        int[] a = new int[N];
        int i;
        int j;
        int[] mult = new int[N];
        int npart = 0;
        int seed;
        int[] table;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITION_RANDOM_TEST");
        Console.WriteLine("  I4_PARTITION_RANDOM generates a random partition.");
        Console.WriteLine("");

        seed = 123456789;
        //
        //  Call once just to get the partition table.
        //
        table = typeMethods.i4_partition_count2(N);

        Console.WriteLine("");
        Console.WriteLine("  The number of partitions of N.");
        Console.WriteLine("");
        Console.WriteLine("     N    Number of partitions");
        Console.WriteLine("");

        for (j = 0; j < N; j++)
        {
            Console.WriteLine((j + 1).ToString().PadLeft(6) + "  "
                                                            + table[j].ToString().PadLeft(6) + "");
        }

        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            typeMethods.i4_partition_random(N, table, ref seed, ref a, ref mult, ref npart);

            typeMethods.i4_partition_print(N, npart, a, mult);

        }
    }

    public static void i4_partitions_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_PARTITIONS_NEXT_TEST tests I4_PARTITIONS_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int[] m = new int[3];
        int msum;
        int s = 3;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("I4_PARTITIONS_NEXT_TEST");
        Console.WriteLine("  I4_PARTITIONS_NEXT produces the next");
        Console.WriteLine("  nondecreasing partitions of an integer, and");
        Console.WriteLine("  if necessary, increments the integer to keep on going.");

        i = 0;
        m[0] = 0;
        m[1] = 0;
        m[2] = 0;

        Console.WriteLine("");
        Console.WriteLine("   I Sum    Partition");
        Console.WriteLine("");
        msum = typeMethods.i4vec_sum(s, m);
        cout = "  " + i.ToString().PadLeft(2)
                    + "  " + msum.ToString().PadLeft(2) + "    ";
        for (j = 0; j < s; j++)
        {
            cout += m[j].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        for (i = 1; i <= 15; i++)
        {
            typeMethods.i4_partitions_next(s, ref m);
            msum = typeMethods.i4vec_sum(s, m);
            cout = "  " + i.ToString().PadLeft(1)
                        + "  " + msum.ToString().PadLeft(1) + "    ";
            for (j = 0; j < s; j++)
            {
                cout += m[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  You can start from any legal partition.");
        Console.WriteLine("  Here, we restart at ( 2, 1, 0 ).");

        i = 0;
        m[0] = 2;
        m[1] = 1;
        m[2] = 0;

        Console.WriteLine("");
        Console.WriteLine("   I Sum    Partition");
        Console.WriteLine("");
        msum = typeMethods.i4vec_sum(s, m);
        cout = "  " + i.ToString().PadLeft(2)
                    + "  " + msum.ToString().PadLeft(2) + "    ";
        for (j = 0; j < s; j++)
        {
            cout += m[j].ToString().PadLeft(2);
        }

        Console.WriteLine(cout);

        for (i = 1; i <= 15; i++)
        {
            typeMethods.i4_partitions_next(s, ref m);
            msum = typeMethods.i4vec_sum(s, m);
            cout = "  " + i.ToString().PadLeft(2)
                        + "  " + msum.ToString().PadLeft(2) + "    ";
            for (j = 0; j < s; j++)
            {
                cout += m[j].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }
    }

    public static void i4_rise_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_RISE_TEST tests I4_RISE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int f1 = 0;
        int f2 = 0;
        int m = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("I4_RISE_TEST");
        Console.WriteLine("  I4_RISE evaluates the rising factorial function.");
        Console.WriteLine("");
        Console.WriteLine("         M         N     Exact  I4_RISE(M,N)");

        n_data = 0;

        while (true)
        {
            typeMethods.i4_rise_values(ref n_data, ref m, ref n, ref f1);

            if (n_data == 0)
            {
                break;
            }

            f2 = typeMethods.i4_rise(m, n);

            Console.WriteLine("  " + m.ToString().PadLeft(8)
                                   + "  " + n.ToString().PadLeft(8)
                                   + "  " + f1.ToString().PadLeft(8)
                                   + "  " + f2.ToString().PadLeft(8) + "");
        }
    }

    public static void i4_sqrt_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SQRT_TEST tests I4_SQRT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int q = 0;
        int r = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_SQRT_TEST");
        Console.WriteLine("  I4_SQRT computes the square root of an integer.");
        Console.WriteLine("");
        Console.WriteLine("       N  Sqrt(N) Remainder");
        Console.WriteLine("");

        for (n = -5; n <= 20; n++)
        {
            typeMethods.i4_sqrt(n, ref q, ref r);

            Console.WriteLine(n.ToString().PadLeft(9) + "  "
                                                      + q.ToString().PadLeft(9) + "  "
                                                      + r.ToString().PadLeft(9) + "");
        }
    }

    public static void i4_sqrt_cf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_SQRT_CF_TEST tests I4_SQRT_CF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAX_TERM = 100;

        int[] b = new int[MAX_TERM + 1];
        int i = 0;
        int n = 0;
        int n_term = 0;

        Console.WriteLine("");
        Console.WriteLine("I4_SQRT_CF_TEST");
        Console.WriteLine("  I4_SQRT_CF computes the continued fraction form");
        Console.WriteLine("  of the square root of an integer.");
        Console.WriteLine("");
        Console.WriteLine("   N  Period  Whole  Repeating Part");
        Console.WriteLine("");

        for (n = 1; n <= 20; n++)
        {
            typeMethods.i4_sqrt_cf(n, MAX_TERM, ref n_term, ref b);
            string cout = n.ToString().PadLeft(5) + "  "
                                                  + n_term.ToString().PadLeft(5) + "  ";
            for (i = 0; i <= n_term; i++)
            {
                cout += b[i].ToString().PadLeft(5) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

    public static void i4_to_dvec_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_DVEC_TEST tests I4_TO_DVEC;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] dvec = new int[6];
        int i;
        int i1;
        int i2;
        int n;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_DVEC_TEST");
        Console.WriteLine("  I4_TO_DVEC converts an I4 to a DVEC;");
        Console.WriteLine("");
        Console.WriteLine("         I4 => DVEC => I4");
        Console.WriteLine("");

        seed = 123456789;
        i1 = UniformRNG.i4_uniform_ab(-10000, 10000, ref seed);

        n = 6;
        typeMethods.i4_to_dvec(i1, n, ref dvec);

        i2 = typeMethods.dvec_to_i4(n, ref dvec);

        string cout = "  " + i1.ToString().PadLeft(6) + "  ";
        for (i = n - 1; 0 <= i; i--)
        {
            cout += dvec[i].ToString().PadLeft(2);
        }

        Console.WriteLine(cout + "  " + i2.ToString().PadLeft(6) + "");
    }

    public static void i4_to_i4poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_I4POLY_TEST tests I4_TO_I4POLY.;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DEGREE_MAX = 5;
        int TEST_NUM = 9;

        int[] a = new int[DEGREE_MAX + 1];
        int base_;
        int[] base_test =  {
                2, 2, 2, 3, 4, 5, 6, 23, 24
            }
            ;
        int degree = 0;
        int i;
        int intval;
        int intval2;
        int[] intval_test =  {
                1, 6, 23, 23, 23, 23, 23, 23, 23
            }
            ;
        int test;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_I4POLY_TEST");
        Console.WriteLine("  I4_TO_I4POLY converts an integer to a polynomial");
        Console.WriteLine("  in a given base;");
        Console.WriteLine("");
        Console.WriteLine("       I    BASE  DEGREE  Coefficients");
        Console.WriteLine("");
        for (test = 0; test < TEST_NUM; test++)
        {
            intval = intval_test[test];
            base_ = base_test[test];
            typeMethods.i4_to_i4poly(intval, base_, DEGREE_MAX, ref degree, ref a);
            cout = "  "
                   + intval.ToString().PadLeft(6) + "  "
                   + base_.ToString().PadLeft(6) + "  "
                   + degree.ToString().PadLeft(6) + "  ";
            for (i = 0; i <= degree; i++)
            {
                cout += a[i].ToString().PadLeft(6) + "  ";
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Now let I4_TO_I4POLY convert I to a polynomial,");
        Console.WriteLine("  use I4POLY_TO_I4 to evaluate it, and compare.");
        Console.WriteLine("");
        Console.WriteLine("       I    I2");
        Console.WriteLine("");
        for (test = 0; test < TEST_NUM; test++)
        {
            intval = intval_test[test];
            base_ = base_test[test];
            typeMethods.i4_to_i4poly(intval, base_, DEGREE_MAX, ref degree, ref a);
            intval2 = typeMethods.i4poly_to_i4(degree, a, base_);
            Console.WriteLine("  "
                              + intval.ToString().PadLeft(6) + "  "
                              + intval2.ToString().PadLeft(6) + "");
        }
    }

    public static void i4_to_van_der_corput_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_VAN_DER_CORPUT_TEST tests I4_TO_VAN_DER_CORPUT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double h;
        int i;
        int j;
        int p;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_VAN_DER_CORPUT_TEST");
        Console.WriteLine("  I4_TO_VAN_DER_CORPUT computes the elements ");
        Console.WriteLine("  of a van der Corput sequence.");
        Console.WriteLine("  The sequence depends on the prime number used");
        Console.WriteLine("  as a base.");
        Console.WriteLine("");
        string cout = "Base: ";
        for (j = 1; j <= 5; j++)
        {
            p = Burkardt.Function.Prime.prime(j);
            cout += p.ToString().PadLeft(10) + "  ";
        }

        Console.WriteLine(cout);
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            cout = i.ToString().PadLeft(4) + "  ";
            for (j = 1; j <= 5; j++)
            {
                p = Burkardt.Function.Prime.prime(j);
                h = Burkardt.Sequence.VanDerCorput.i4_to_van_der_corput(i, p);
                cout += h.ToString().PadLeft(10) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

    public static void i4mat_01_rowcolsum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_02_ROWCOLSUM_TEST tests I4MAT_01_ROWCOLSUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 5;
        int N = 5;

        int[] a = new int[M * N];
        int[] c =  {
                2, 2, 2, 2, 1
            }
            ;
        bool error = false;
        int[] r =  {
                3, 2, 2, 1, 1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("I4MAT_01_ROWCOLSUM_TEST");
        Console.WriteLine("  I4MAT_01_ROWCOLSUM constructs a 01 matrix with");
        Console.WriteLine("  given row and column sums.");

        typeMethods.i4vec1_print(M, r, "  The rowsum vector:");
        typeMethods.i4vec1_print(N, c, "  The columnsum vector: ");

        typeMethods.i4mat_01_rowcolsum(M, N, r, c, ref a, ref error);

        switch (error)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  I4MAT_01_ROWCOLSUM returned error flag.");
                break;
            default:
                typeMethods.i4mat_print(M, N, a, "  The rowcolsum matrix:");
                break;
        }
    }

    public static void i4mat_u1_inverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_U1_INVERSE_TEST tests I4MAT_U1_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        int[] a =  {
                1, 0, 0, 0, 0, 0,
                1, 1, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0,
                0, 0, 1, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                75, 0, 0, 0, 1, 1
            }
            ;
        int[] b = new int [N * N];

        Console.WriteLine("");
        Console.WriteLine("I4MAT_U1_INVERSE_TEST");
        Console.WriteLine("  I4MAT_U1_INVERSE inverts a unit upper triangular matrix.");

        typeMethods.i4mat_print(N, N, a, "  The input matrix:");

        typeMethods.i4mat_u1_inverse(N, a, ref b);

        typeMethods.i4mat_print(N, N, b, "  The inverse matrix:");
    }

    public static void i4mat_perm0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_PERM0_TEST tests I4MAT_PERM0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 9;

        int[] a = new int[N * N];
        int i;
        int j;
        int[] p =  {
                1,2,8,5,6,7,4,3,0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("I4MAT_PERM0_TEST");
        Console.WriteLine("  I4MAT_PERM0 reorders an integer matrix in place.");
        Console.WriteLine("  The rows and columns use the same permutation.");

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * N] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.i4mat_print(N, N, a, "  The input matrix:");

        Permutation.perm0_print(N, p, "  The row and column permutation:");

        typeMethods.i4mat_perm0(N, ref a, p);

        typeMethods.i4mat_print(N, N, a, "  The permuted matrix:");
    }

    public static void i4mat_2perm0_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4MAT_2PERM0_TEST tests I4MAT_2PERM0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 9;
        int N = 7;

        int[] a = new int[M * N];
        int i;
        int j;
        int[] p =  {
                1,2,8,5,6,7,4,3,0
            }
            ;
        int[] q =  {
                2,3,4,5,6,0,1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("I4MAT_2PERM0_TEST");
        Console.WriteLine("  I4MAT_2PERM0 reorders an integer matrix in place.");
        Console.WriteLine("  Rows and columns use different permutations.");

        for (i = 0; i < M; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[i + j * M] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.i4mat_print(M, N, a, "  The input matrix:");

        Permutation.perm0_print(M, p, "  The row permutation:");

        Permutation.perm0_print(N, q, "  The column permutation:");

        typeMethods.i4mat_2perm0(M, N, a, p, q);

        typeMethods.i4mat_print(M, N, a, "  The permuted matrix:");
    }

    public static void i4poly_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_TEST test I4POLY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 6;

        int[] a = new int[N];
        int iopt = 0;
        int test = 0;
        int val = 0;
        int x0 = 0;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_TEST");
        Console.WriteLine("  I4POLY converts between power sum, factorial");
        Console.WriteLine("  and Taylor forms, and can evaluate a polynomial");
        Console.WriteLine("");

        for (test = 1; test <= 6; test++)
        {
            switch (test)
            {
                case 1:
                    iopt = -3;
                    break;
                case 2:
                    iopt = -2;
                    break;
                case 3:
                    iopt = -1;
                    x0 = 2;
                    break;
                case 4:
                    iopt = 0;
                    x0 = 2;
                    break;
                case 5:
                    iopt = 6;
                    x0 = 2;
                    break;
                case 6:
                    iopt = 6;
                    x0 = -2;
                    break;
            }

            a[0] = 0;
            a[1] = 0;
            a[2] = 0;
            a[3] = 0;
            a[4] = 0;
            a[5] = 1;

            switch (test)
            {
                case 1:
                    typeMethods.i4vec1_print(N, a, "  All calls have input A as follows:");
                    break;
            }

            typeMethods.i4poly(N, ref a, x0, iopt, ref val);

            Console.WriteLine("");
            Console.WriteLine("  Option IOPT = " + iopt + "");

            switch (iopt)
            {
                case >= -1:
                    Console.WriteLine("  X0 = " + x0 + "");
                    break;
            }

            switch (iopt)
            {
                case -3:
                case -2:
                case > 0:
                    typeMethods.i4vec1_print(N, a, "  Output array:");
                    break;
                case -1:
                case 0:
                    Console.WriteLine("  Value = " + val + "");
                    break;
            }
        }
    }

    public static void i4poly_add_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_ADD_TEST tests I4POLY_ADD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a =  {
                0, 1, 2, 3, 4, 5
            }
            ;
        int[] b =  {
                1, -2, 7, 8, 0, -5
            }
            ;
        int[] c;
        int na = 5;
        int nb = 5;
        int nc;
        int nc2;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_ADD_TEST");
        Console.WriteLine("  I4POLY_ADD adds two polynomials.");

        typeMethods.i4poly_print(na, a, "  Polynomial A:");

        typeMethods.i4poly_print(nb, b, "  Polynomial B:");

        c = typeMethods.i4poly_add(na, a, nb, b);

        nc = Math.Max(na, nb);

        nc2 = typeMethods.i4poly_degree(nc, c);

        typeMethods.i4poly_print(nc2, c, "  Polynomial C = A+B:");
    }

    public static void i4poly_cyclo_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_CYCLO_TEST tests I4POLY_CYCLO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 10;

        int[] phi = new int [N_MAX + 1];
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_CYCLO_TEST");
        Console.WriteLine("  I4POLY_CYCLO computes cyclotomic polynomials.");

        for (n = 0; n <= N_MAX; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("");

            typeMethods.i4poly_cyclo(n, phi);

            typeMethods.i4poly_print(n, phi, "  The cyclotomic polynomial:");
        }

    }

    public static void i4poly_degree_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_DEGREE_TEST tests I4POLY_DEGREE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = {
                0, 1, 0, 3, 4, 0, 6, 7, 0, 0, 0
            }
            ;
        int degree;
        int n = 4;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_DEGREE_TEST");
        Console.WriteLine("  I4POLY_DEGREE determines the degree of an I4POLY.");

        typeMethods.i4poly_print(n, a, "  The polynomial:");

        degree = typeMethods.i4poly_degree(n, a);

        Console.WriteLine("");
        Console.WriteLine("The polyomial degree = " + degree + "");
    }

    public static void i4poly_dif_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_DIF_TEST tests I4POLY_DIF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[11];
        int[] b;
        int d = 0;
        int na = 0;
        int test_num = 2;
        int test;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_DIF_TEST");
        Console.WriteLine("  I4POLY_DIF computes derivatives of an I4POLY.");
        Console.WriteLine("");
        //
        //  1: Differentiate X^3 + 2*X^2 - 5*X - 6 once.
        //  2: Differentiate X^4 + 3*X^3 + 2*X^2 - 2  3 times.
        //
        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    na = 3;
                    d = 1;
                    a[0] = -6;
                    a[1] = -5;
                    a[2] = 2;
                    a[3] = 1;
                    break;
                case 2:
                    na = 4;
                    d = 3;
                    a[0] = -2;
                    a[1] = 5;
                    a[2] = 2;
                    a[3] = 3;
                    a[4] = 1;
                    break;
            }

            typeMethods.i4poly_print(na, a, "  The polynomial A:");

            Console.WriteLine("");
            Console.WriteLine("  Differentiate A " + d + " times.");

            b = typeMethods.i4poly_dif(na, a, d);
            typeMethods.i4poly_print(na - d, b, "  The derivative, B:");
        }
    }

    public static void i4poly_div_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_DIV_TEST tests I4POLY_DIV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = new int[11];
        int[] b = new int[11];
        int na = 0;
        int nb = 0;
        int nq = 0;
        int nr = 0;
        int ntest = 2;
        int[] q = new int[11];
        int[] r = new int[11];
        int test = 0;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_DIV_TEST");
        Console.WriteLine("  I4POLY_DIV computes the quotient and");
        Console.WriteLine("  remainder for polynomial division.");
        Console.WriteLine("");
        //
        //  1: Divide X^3 + 2*X^2 - 5*X - 6  by X-2.  
        //     Quotient is 3+4*X+X^2, remainder is 0.
        //
        //  2: Divide X^4 + 3*X^3 + 2*X^2 - 2  by  X^2 + X - 3.
        //     Quotient is X^2 + 2*X + 3, remainder 8*X + 7.
        //
        for (test = 1; test <= ntest; test++)
        {
            switch (test)
            {
                case 1:
                    na = 3;
                    a[0] = -6;
                    a[1] = -5;
                    a[2] = 2;
                    a[3] = 1;

                    nb = 1;
                    b[0] = -2;
                    b[1] = 1;
                    break;
                case 2:
                    na = 4;
                    a[0] = -2;
                    a[1] = 5;
                    a[2] = 2;
                    a[3] = 3;
                    a[4] = 1;
                    nb = 2;
                    b[0] = -3;
                    b[1] = 1;
                    b[2] = 1;
                    break;
            }

            typeMethods.i4poly_print(na, a, "  The polynomial to be divided, A:");
            typeMethods.i4poly_print(nb, b, "  The divisor polynomial, B:");

            typeMethods.i4poly_div(na, a, nb, b, ref nq, ref q, ref nr, ref r);

            typeMethods.i4poly_print(nq, q, "  The quotient polynomial, Q:");
            typeMethods.i4poly_print(nr, r, "  The remainder polynomial, R:");
        }
    }

    public static void i4poly_mul_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_MUL_TEST tests I4POLY_MUL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int MAXN = 5;

        int[] a = new int[MAXN + 1];
        int[] b = new int[MAXN + 1];
        int[] c = new int[MAXN + 1];
        int na = 0;
        int nb = 0;
        int ntest = 2;
        int test = 0;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_MUL_TEST");
        Console.WriteLine("  I4POLY_MUL multiplies two polynomials.");
        Console.WriteLine("");
        //
        //  1: Multiply (1+X) times (1-X).  Answer is 1-X^2.
        //  2: Multiply (1+2*X+3*X^2) by (1-2*X). Answer is 1 + 0*X - X^2 - 6*X^3
        //
        for (test = 1; test <= ntest; test++)
        {
            switch (test)
            {
                case 1:
                    na = 1;
                    a[0] = 1;
                    a[1] = 1;
                    nb = 1;
                    b[0] = 1;
                    b[1] = -1;
                    break;
                case 2:
                    na = 2;
                    a[0] = 1;
                    a[1] = 2;
                    a[2] = 3;
                    nb = 1;
                    b[0] = 1;
                    b[1] = -2;
                    break;
            }

            typeMethods.i4poly_mul(na, a, nb, b, ref c);

            typeMethods.i4poly_print(na, a, "  The factor A:");

            typeMethods.i4poly_print(nb, b, "  The factor B:");

            typeMethods.i4poly_print(na + nb, c, "  The product C = A*B:");
        }
    }

    public static void i4poly_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_PRINT_TEST tests I4POLY_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a = {
                -2, 5, 2, 3, 1
            }
            ;
        int n = 4;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_PRINT_TEST");
        Console.WriteLine("  I4POLY_PRINT prints an I4POLY.");

        typeMethods.i4poly_print(n, a, "  The polynomial:");
    }

    public static void i4poly_to_i4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4POLY_TO_I4_TEST tests I4POLY_TO_I4;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DEGREE_MAX = 5;
        int TEST_NUM = 9;

        int[] a = new int[DEGREE_MAX + 1];
        int base_ = 0;
        int[] base_test =  {
                2, 2, 2, 3, 4, 5, 6, 23, 24
            }
            ;
        int degree = 0;
        int i = 0;
        int intval = 0;
        int intval2 = 0;
        int[] intval_test =  {
                1, 6, 23, 23, 23, 23, 23, 23, 23
            }
            ;
        int test = 0;

        Console.WriteLine("");
        Console.WriteLine("I4POLY_TO_I4_TEST");
        Console.WriteLine("  I4POLY_TO_I4 evaluates an integer polynomial");
        Console.WriteLine("  at a given point;");
        Console.WriteLine("");
        Console.WriteLine("       I    BASE  DEGREE  Coefficients");
        Console.WriteLine("");
        for (test = 0; test < TEST_NUM; test++)
        {
            intval = intval_test[test];
            base_ = base_test[test];
            typeMethods.i4_to_i4poly(intval, base_, DEGREE_MAX, ref degree, ref a);
            string cout = "  "
                          + intval.ToString().PadLeft(6) + "  "
                          + base_.ToString().PadLeft(6) + "  "
                          + degree.ToString().PadLeft(6) + "  ";
            for (i = 0; i <= degree; i++)
            {
                cout += a[i].ToString().PadLeft(6) + "  ";
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Now let I4_TO_I4POLY convert I to a polynomial,");
        Console.WriteLine("  use I4POLY_TO_I4 to evaluate it, and compare.");
        Console.WriteLine("");
        Console.WriteLine("       I    I2");
        Console.WriteLine("");
        for (test = 0; test < TEST_NUM; test++)
        {
            intval = intval_test[test];
            base_ = base_test[test];
            typeMethods.i4_to_i4poly(intval, base_, DEGREE_MAX, ref degree, ref a);
            intval2 = typeMethods.i4poly_to_i4(degree, a, base_);
            Console.WriteLine("  "
                              + intval.ToString().PadLeft(6) + "  "
                              + intval2.ToString().PadLeft(6) + "");
        }
    }

    public static void i4vec_backtrack_test()

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

    public static void i4vec_descends_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_DESCENDS_TEST tests I4VEC_DESCENDS;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        int[] a;
        int i;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_DESCENDS_TEST");
        Console.WriteLine("  I4VEC_DESCENDS is true if an integer vector decreases.");
        Console.WriteLine("");

        seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            a = UniformRNG.i4vec_uniform_ab_new(N, 1, N, ref seed);

            typeMethods.i4vec1_print(N, a, "  The integer array to search:");

            if (typeMethods.i4vec_descends(N, ref a))
            {
                Console.WriteLine("  The preceding vector is descending.");
            }
            else
            {
                Console.WriteLine("  The preceding vector is not descending.");
            }
        }
    }

    public static void i4vec_frac_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_FRAC_TEST tests I4VEC_FRAC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int[] a;
        int afrac;
        int k;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_FRAC_TEST");
        Console.WriteLine("  I4VEC_FRAC: K-th smallest integer vector entry.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 1, 2 * N, ref seed);

        typeMethods.i4vec1_print(N, a, "  The integer array to search:");

        Console.WriteLine("");
        Console.WriteLine("     K   K-th smallest");
        Console.WriteLine("");

        for (k = 1; k <= N; k++)
        {
            afrac = typeMethods.i4vec_frac(N, ref a, k);

            Console.WriteLine(k.ToString().PadLeft(6) + "  "
                                                      + afrac.ToString().PadLeft(6) + "");

        }
    }

    public static void i4vec_index_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_INDEX_TEST tests I4VEC_INDEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        int[] a;
        int aval;
        int first;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_INDEX_TEST");
        Console.WriteLine("  I4VEC_INDEX returns the index of the first occurrence");
        Console.WriteLine("  of a given value in an integer vector.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 1, N / 2, ref seed);

        aval = a[N / 2];

        typeMethods.i4vec1_print(N, a, "  The integer array to search:");

        first = typeMethods.i4vec_index(N, a, aval);

        Console.WriteLine("");
        Console.WriteLine("  The value searched for is " + aval + "");
        Console.WriteLine("  The index of first occurrence is " + first + "");
    }

    public static void i4vec_maxloc_last_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_MAXLOC_LAST_TEST tests I4VEC_MAXLOC_LAST;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        int[] a;
        int last;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_MAXLOC_LAST_TEST");
        Console.WriteLine("  I4VEC_MAXLOC_LAST: index of the last maximal");
        Console.WriteLine("  entry in an integer vector.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 1, N / 4, ref seed);

        typeMethods.i4vec1_print(N, a, "  The integer array to search:");

        last = typeMethods.i4vec_maxloc_last(N, a);

        Console.WriteLine("");
        Console.WriteLine("  Index of last maximal entry is " + last + "");
    }

    public static void i4vec_pairwise_prime_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_PAIRWISE_PRIME_TEST tests I4VEC_PAIRWISE_PRIME;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        int[] a;
        int i;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_PAIRWISE_PRIME_TEST");
        Console.WriteLine("  I4VEC_PAIRWISE_PRIME is true if an integer vector");
        Console.WriteLine("  is pairwise prime.");
        Console.WriteLine("");

        seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            a = UniformRNG.i4vec_uniform_ab_new(N, 1, N, ref seed);

            typeMethods.i4vec1_print(N, a, "  The array to check:");

            if (typeMethods.i4vec_pairwise_prime(N, a))
            {
                Console.WriteLine("  The preceding vector is pairwise prime.");
            }
            else
            {
                Console.WriteLine("  The preceding vector is not pairwise prime.");
            }
        }
    }

    public static void i4vec_reverse_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_REVERSE_TEST tests I4VEC_REVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int[] a;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_REVERSE_TEST");
        Console.WriteLine("  I4VEC_REVERSE reverses an integer vector.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 1, N, ref seed);

        typeMethods.i4vec1_print(N, a, "  The integer array:");

        typeMethods.i4vec_reverse(N, ref a);

        typeMethods.i4vec1_print(N, a, "  The reversed integer array:");
    }

    public static void i4vec_sort_bubble_a_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SORT_BUBBLE_A_TEST tests I4VEC_SORT_BUBBLE_A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        int[] a;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_BUBBLE_A_TEST");
        Console.WriteLine("  I4VEC_SORT_BUBBLE_A ascending sorts an integer vector");
        Console.WriteLine("  using bubble sort.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 0, 3 * N, ref seed);

        typeMethods.i4vec1_print(N, a, "  Unsorted array:");

        typeMethods.i4vec_sort_bubble_a(N, ref a);

        typeMethods.i4vec1_print(N, a, "  Sorted array:");
    }

    public static void i4vec_sort_heap_index_d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_SORT_HEAP_INDEX_D_TEST tests I4VEC_SORT_HEAP_INDEX_D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 20;

        int[] a;
        int i;
        int[] indx;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_SORT_HEAP_INDEX_D_TEST");
        Console.WriteLine("  I4VEC_SORT_HEAP_INDEX_D descending index-sorts");
        Console.WriteLine("  an integer vector using heap sort.");
        Console.WriteLine("");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new(N, 0, 3 * N, ref seed);

        typeMethods.i4vec1_print(N, a, "  Unsorted array:");

        indx = typeMethods.i4vec_sort_heap_index_d(N, a);

        Console.WriteLine("");
        Console.WriteLine("     I  INDX[I]  A[INDX[I]-1]");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + indx[i].ToString().PadLeft(6) + "  "
                                                      + a[indx[i]].ToString().PadLeft(6) + "");
        }
    }

    public static void i4vec_transpose_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_TRANSPOSE_PRINT_TEST tests I4VEC_TRANSPOSE_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("I4VEC_TRANSPOSE_PRINT_TEST");
        Console.WriteLine("  I4VEC_TRANSPOSE_PRINT prints an integer vector");
        Console.WriteLine("  with 5 entries to a row, and an optional title.");

        n = 12;
        a = typeMethods.i4vec_indicator1_new(n);

        typeMethods.i4vec_transpose_print(n, a, "  My array:  ");

    }

    public static void i4vec_uniform_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC_UNIFORM_AB_TEST tests I4VEC_UNIFORM_AB_NEW.
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
        Console.WriteLine("I4VEC_UNIFORM_AB_TEST");
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