using System;
using Burkardt.Types;

namespace SubsetTestNS;

public static class TupleTest
{
    public static void tuple_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT_TEST tests TUPLE_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 2;

        const int m1 = 2;
        const int m2 = 4;
        int[] x = new int[N];

        Console.WriteLine("");
        Console.WriteLine("TUPLE_NEXT_TEST");
        Console.WriteLine("  TUPLE_NEXT returns the next \"tuple\", that is,");
        Console.WriteLine("  a vector of N integers, each between M1 and M2.");
        Console.WriteLine("");
        Console.WriteLine("  M1 = " + m1 + "");
        Console.WriteLine("  M2 = " + m2 + "");
        Console.WriteLine("  N = " + N + "");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            BTuple.tuple_next(m1, m2, N, ref rank, ref x);

            if (rank == 0)
            {
                break;
            }

            string cout = rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < N; i++)
            {
                cout += x[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

        }
    }

    public static void tuple_next_fast_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT_FAST_TEST tests TUPLE_NEXT_FAST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] base_ = new int[2];
        const int m = 3;
        const int n = 2;
        int[] x = new int[2];

        Console.WriteLine("");
        Console.WriteLine("TUPLE_NEXT_FAST_TEST");
        Console.WriteLine("  TUPLE_NEXT_FAST returns the next \"tuple\", that is,");
        Console.WriteLine("  a vector of N integers, each between 1 and M.");
        Console.WriteLine("");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("");
        //
        //  Initialize.
        //
        int rank = -1;
        BTuple.tuple_next_fast(m, n, rank, ref base_, ref x);

        int rank_hi = (int)Math.Pow(m, n);

        for (rank = 0; rank < rank_hi; rank++)
        {
            BTuple.tuple_next_fast(m, n, rank, ref base_, ref x);

            string cout = rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < n; i++)
            {
                cout += x[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

    public static void tuple_next_ge_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT_GE_TEST tests TUPLE_NEXT_GE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 3;

        const int m = 3;
        int[] x = new int[N];

        Console.WriteLine("");
        Console.WriteLine("TUPLE_NEXT_GE_TEST");
        Console.WriteLine("  TUPLE_NEXT_GE returns the next nondecreasting \"tuple\",");
        Console.WriteLine("  that is, a vector of N integers, each between 1 and M,");
        Console.WriteLine("  with the additional property that the digits never decrease");
        Console.WriteLine("  reading from left to right.");
        Console.WriteLine("");
        Console.WriteLine("  M = " + m + "");
        Console.WriteLine("  N = " + N + "");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            BTuple.tuple_next_ge(m, N, ref rank, ref x);

            if (rank == 0)
            {
                break;
            }

            string cout = rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < N; i++)
            {
                cout += x[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

        }

    }

    public static void tuple_next2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT2_TEST tests TUPLE_NEXT2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 3;

        int[] x = new int[N];
        int[] xmin = { 2, 3, 8 };
        int[] xmax = { 4, 3, 5 };

        Console.WriteLine("");
        Console.WriteLine("TUPLE_NEXT2_TEST");
        Console.WriteLine("  TUPLE_NEXT2 returns the next \"tuple\",");
        Console.WriteLine("  that is, a vector of N integers.");
        Console.WriteLine("  Each position in the vector has a separate min and max.");
        Console.WriteLine("  reading from left to right.");
        Console.WriteLine("");
        Console.WriteLine("  N = " + N + "");
        Console.WriteLine("");
        typeMethods.i4vec1_print(N, xmin, "  The minimum values:");
        typeMethods.i4vec1_print(N, xmax, "  The maximum values:");

        Console.WriteLine("");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            BTuple.tuple_next2(N, xmin, xmax, ref x, ref rank);

            if (rank == 0)
            {
                break;
            }

            string cout = rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < N; i++)
            {
                cout += x[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);

        }
    }
}