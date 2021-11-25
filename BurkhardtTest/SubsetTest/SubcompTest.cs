using System;
using Burkardt.Composition;

namespace SubsetTestNS;

public static class SubcompTest
{
    public static void subcomp_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCOMP_NEXT_TEST tests SUBCOMP_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 6;
        const int k = 3;
        int[] a = new int[k];
        bool more = false;
        int h = 0;
        int t = 0;

        SubCompData data = new();

        Console.WriteLine("");
        Console.WriteLine("SUBCOMP_NEXT_TEST");
        Console.WriteLine("  SUBCOMP_NEXT generates subcompositions.");
        Console.WriteLine("");
        Console.WriteLine("  Seek all subcompositions of N = " + n + "");
        Console.WriteLine("  using K = " + k + " parts.");
        Console.WriteLine("");
        Console.WriteLine("     #   Sum");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            SubComp.subcomp_next(ref data, n, k, ref a, ref more, ref h, ref t);

            int total = 0;
            int i;
            for (i = 0; i < k; i++)
            {
                total += a[i];
            }

            rank += 1;
            string cout = "  " + rank.ToString().PadLeft(4)
                               + "  " + total.ToString().PadLeft(4)
                               + "  ";

            for (i = 0; i < k; i++)
            {
                cout += a[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            if (!more)
            {
                break;
            }
        }
    }

    public static void subcompnz_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCOMPNZ_NEXT_TEST tests SUBCOMPNZ_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 December 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 6;
        const int k = 3;
        int[] a = new int[k];
        bool more = false;
        int h = 0;
        int t = 0;
        int n2 = 0;
        bool more2 = false;
        CompNZData data = new();

        Console.WriteLine("");
        Console.WriteLine("SUBCOMPNZ_NEXT_TEST");
        Console.WriteLine("  SUBCOMPNZ_NEXT generates subcompositions using nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("  Seek all subcompositions of N = " + n + "");
        Console.WriteLine("  using K = " + k + " nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("     #   Sum");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            SubComp.subcompnz_next(ref data, n, k, ref a, ref more, ref h, ref t, ref n2, ref more2);

            int total = 0;
            int i;
            for (i = 0; i < k; i++)
            {
                total += a[i];
            }

            rank += 1;
            string cout = "  " + rank.ToString().PadLeft(4)
                               + "  " + total.ToString().PadLeft(4)
                               + "  ";

            for (i = 0; i < k; i++)
            {
                cout += a[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            if (!more)
            {
                break;
            }
        }
    }

    public static void subcompnz2_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCOMPNZ2_NEXT_TEST tests SUBCOMPNZ2_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n_lo = 5;
        const int n_hi = 7;
        const int k = 3;
        int[] a = new int[k];
        bool more = false;
        int h = 0;
        int t = 0;
        int n2 = 0;
        bool more2 = false;
        CompNZData data = new();

        Console.WriteLine("");
        Console.WriteLine("SUBCOMPNZ2_NEXT_TEST");
        Console.WriteLine("  SUBCOMPNZ2_NEXT generates subcompositions using nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("  Seek all subcompositions of N");
        Console.WriteLine("  using K = " + k + " nonzero parts.");
        Console.WriteLine("");
        Console.WriteLine("  Here N is in the range " + n_lo + " <= N <= " + n_hi + "");
        Console.WriteLine("");
        Console.WriteLine("     #     N");
        Console.WriteLine("");

        int rank = 0;

        for (;;)
        {
            SubComp.subcompnz2_next(ref data, n_lo, n_hi, k, ref a, ref more, ref h, ref t, ref n2, ref more2);

            int n = 0;
            int i;
            for (i = 0; i < k; i++)
            {
                n += a[i];
            }

            rank += 1;
            string cout = "  " + rank.ToString().PadLeft(4)
                               + "  " + n.ToString().PadLeft(4)
                               + "  ";

            for (i = 0; i < k; i++)
            {
                cout += a[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);

            if (!more)
            {
                break;
            }
        }
    }
}