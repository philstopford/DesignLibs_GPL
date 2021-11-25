﻿using System;
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
        const int N = 4;

        int[] v = new int[N];
        int[] vmax = new int[N];

        Console.WriteLine("");
        Console.WriteLine("REGRO_NEXT_TEST");
        Console.WriteLine("  REGRO_NEXT generates all restricted growth");
        Console.WriteLine("  functions.");
        Console.WriteLine("");

        int rank = 0;

        bool done = true;

        for (;;)
        {
            RestrictedGrowth.regro_next(ref done, N, ref v, ref vmax);

            if (done)
            {
                break;
            }

            rank += 1;
            string cout = "  "
                          + rank.ToString().PadLeft(3) + "  ";
            int i;
            for (i = 0; i < N; i++)
            {
                cout += v[i].ToString().PadLeft(1) + "  ";
            }

            Console.WriteLine(cout);

        }
    }

}