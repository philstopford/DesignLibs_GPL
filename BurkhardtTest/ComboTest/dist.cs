using System;
using Burkardt.RankingNS;

namespace ComboTest;

internal partial class Program
{
    private static void dist_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIST_ENUM_TEST tests DIST_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;

        Console.WriteLine("");
        Console.WriteLine("DIST_ENUM_TEST");
        Console.WriteLine("  DIST_ENUM enumerates distributions of N indistinguishable");
        Console.WriteLine("  objects among M distinguishable slots:");
        Console.WriteLine("");
        Console.WriteLine("      N:      0       1       2       3       4       5");
        Console.WriteLine("   M");

        for (m = 0; m <= 10; m++)
        {
            string cout = "  " + m.ToString().PadLeft(2)
                               + ":  ";
            for (n = 0; n <= 5; n++)
            {
                cout += "  " + Ranking.dist_enum(m, n).ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void dist_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIST_NEXT_TEST tests DIST_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int idist;
        int k;
        int leftmost;
        int m;
        bool more;
        int[] q;

        k = 3;
        m = 5;
        q = new int[k];
        leftmost = 0;
        more = false;

        Console.WriteLine("");
        Console.WriteLine("DIST_NEXT_TEST");
        Console.WriteLine("  DIST_NEXT produces the next distribution of M indistinguishable");
        Console.WriteLine("  objects among K distinguishable slots.");
        Console.WriteLine("");
        Console.WriteLine("  Number of:");
        Console.WriteLine("    indistinguishable objects = " + m + "");
        Console.WriteLine("    distinguishable slots =     " + k + "");
        Console.WriteLine("    distributions is            " + Ranking.dist_enum(k, m) + "");
        Console.WriteLine("");

        idist = 0;

        for (;;)
        {
            Ranking.dist_next(k, m, ref q, ref leftmost, ref more);

            if (!more)
            {
                break;
            }

            idist += 1;
            string cout = "  " + idist.ToString().PadLeft(4);
            for (i = 0; i < k; i++)
            {
                cout += "  " + q[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }
    }
}