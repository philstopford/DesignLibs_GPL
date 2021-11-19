using System;
using Burkardt.RankingNS;

namespace ComboTest;

internal partial class Program
{
    private static void mountain_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOUNTAIN_TEST tests MOUNTAIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 5;
        int x;
        int y;

        Console.WriteLine("");
        Console.WriteLine("MOUNTAIN_TEST");
        Console.WriteLine("  MOUNTAIN computes mountain numbers.");
        Console.WriteLine("");
        Console.WriteLine("  Y  MXY");
        Console.WriteLine("");

        for ( y = 0; y <= n; y++ )
        {
            string cout = "  " + y.ToString().PadLeft(2) + "   ";

            for ( x = 0; x <= 2 * n; x++ )
            {
                cout += "  " + Ranking.mountain ( n, x, y ).ToString().PadLeft(4);
            }
            Console.WriteLine(cout);
        }
    }
}