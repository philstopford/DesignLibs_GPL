using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class motzkinTest
{
    public static void motzkin_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOTZKIN_TEST tests MOTZKIN.
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
        int N = 10;

        int[] a = new int[N + 1];
        int i;

        Console.WriteLine("");
        Console.WriteLine("MOTZKIN_TEST");
        Console.WriteLine("  MOTZKIN computes the Motzkin numbers A(0:N).");
        Console.WriteLine("  A(N) counts the paths from (0,0) to (N,0).");
        Console.WriteLine("");
        Console.WriteLine("  I,  A(I)");
        Console.WriteLine("");

        Motzkin.motzkin(N, ref a);

        for (i = 0; i <= N; i++)
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }
}