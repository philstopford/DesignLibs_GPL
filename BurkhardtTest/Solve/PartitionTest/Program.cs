using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace PartitionTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PARTITION_PROBLEM_TEST.
        //
        //  Discussion:
        //
        //    PARTITION_PROBLEM_TEST tests the PARTITION_PROBLEM library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 0;
        int test;
        int test_num = 5;
        int[] w = null;
        int[] w1 = { 19, 17, 13, 9, 6 };
        int[] w2 = { 484, 114, 205, 288, 506, 503, 201, 127, 410 };
        int[] w3 = { 771, 121, 281, 854, 885, 734, 486, 1003, 83, 62 };
        int[] w4 = { 2, 10, 3, 8, 5, 7, 9, 5, 3, 2 };
        int[] w5 = { 3, 4, 3, 1, 3, 2, 3, 2, 1 };

        Console.WriteLine("");
        Console.WriteLine("PARTITION_PROBLEM_TEST:");
        Console.WriteLine("  Test the PARTITION_PROBLEM library.");
        //
        //  Find individual solutions.
        //
        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    n = 5;
                    w = typeMethods.i4vec_copy_new(n, w1);
                    break;
                case 2:
                    n = 9;
                    w = typeMethods.i4vec_copy_new(n, w2);
                    break;
                case 3:
                    n = 10;
                    w = typeMethods.i4vec_copy_new(n, w3);
                    break;
                case 4:
                    n = 10;
                    w = typeMethods.i4vec_copy_new(n, w4);
                    break;
                case 5:
                    n = 9;
                    w = typeMethods.i4vec_copy_new(n, w5);
                    break;
            }

            partition_brute_test(n, w);
        }

        //
        //  Count solutions.
        //
        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    n = 5;
                    w = typeMethods.i4vec_copy_new(n, w1);
                    break;
                case 2:
                    n = 9;
                    w = typeMethods.i4vec_copy_new(n, w2);
                    break;
                case 3:
                    n = 10;
                    w = typeMethods.i4vec_copy_new(n, w3);
                    break;
                case 4:
                    n = 10;
                    w = typeMethods.i4vec_copy_new(n, w4);
                    break;
                case 5:
                    n = 9;
                    w = typeMethods.i4vec_copy_new(n, w5);
                    break;
            }

            partition_count_test(n, w);
        }

        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PARTITION_PROBLEM_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void partition_brute_test(int n, int[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_BRUTE_TEST tests PARTITION_BRUTE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of weights.
        //
        //    Input, int W[N], a set of weights.
        //
    {
        int[] c;
        int discrepancy = -1;
        int i;
        int w0_sum;
        int w1_sum;

        Console.WriteLine("");
        Console.WriteLine("PARTITION_BRUTE_TEST:");
        Console.WriteLine("  PARTITION_BRUTE_TEST partitions a set of N integers W so that the subsets");
        Console.WriteLine("  have equal sums.");

        c = new int[n];

        Partition.partition_brute(n, w, ref c, ref discrepancy);

        Console.WriteLine("");
        Console.WriteLine("     I        W0        W1");
        Console.WriteLine("");
        w0_sum = 0;
        w1_sum = 0;
        for (i = 0; i < n; i++)
        {
            switch (c[i])
            {
                case 0:
                    w0_sum += w[i];
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                    break;
                default:
                    w1_sum += w[i];
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + "        "
                                           + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                    break;
            }
        }

        Console.WriteLine("        --------  --------");
        Console.WriteLine("  " + "    "
                               + "  " + w0_sum.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                               + "  " + w1_sum.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        Console.WriteLine("");
        Console.WriteLine("  Discrepancy = " + discrepancy.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

    }

    private static void partition_count_test(int n, int[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_COUNT_TEST tests PARTITION_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of weights.
        //
        //    Input, int W[N], a set of weights.
        //
    {
        int count;
        int i;

        Console.WriteLine("");
        Console.WriteLine("PARTITION_COUNT_TEST:");
        Console.WriteLine("  PARTITION_COUNT counts the number of exact solutions");
        Console.WriteLine("  of the partition problem.");

        count = Partition.partition_count(n, w);

        Console.WriteLine("");
        Console.WriteLine("     I        W");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + w[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of solutions = " + count + "");

    }
}