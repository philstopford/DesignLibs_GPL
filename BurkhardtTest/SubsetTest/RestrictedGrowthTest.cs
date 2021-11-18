using System;
using Burkardt.SolveNS;

namespace SubsetTestNS;

public static class RestrictedGrowthTest
{
    public static void regro_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REGRO_NEXT_TEST tests REGRO_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 4;

        bool done;
        int i;
        int rank;
        int[] v = new int[N];
        int[] vmax = new int[N];

        Console.WriteLine("");
        Console.WriteLine("REGRO_NEXT_TEST");
        Console.WriteLine("  REGRO_NEXT generates all restricted growth");
        Console.WriteLine("  functions.");
        Console.WriteLine("");

        rank = 0;

        done = true;

        for (;;)
        {
            RestrictedGrowth.regro_next(ref done, N, ref v, ref vmax);

            if (done)
            {
                break;
            }

            rank += 1;
            string cout = "  "
                          + rank.ToString(CultureInfo.InvariantCulture).PadLeft(3) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += v[i].ToString(CultureInfo.InvariantCulture).PadLeft(1) + "  ";
            }

            Console.WriteLine(cout);

        }
    }

}