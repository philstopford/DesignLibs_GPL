using System;
using Burkardt.Sequence;

namespace SubsetTestNS;

public static class PerrinTest
{
    public static void perrin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERRIN_TEST tests PERRIN;
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
        const int N = 10;

        int i;
        int[] p = new int[N];

        Console.WriteLine("");
        Console.WriteLine("PERRIN_TEST");
        Console.WriteLine("  PERRIN computes the Perrin numbers.");
        Console.WriteLine("");

        Perrin.perrin(N, ref p);

        Console.WriteLine("");
        Console.WriteLine("   N    P(N)");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(4) + "  "
                              + p[i].ToString().PadLeft(6) + "");
        }
    }
}