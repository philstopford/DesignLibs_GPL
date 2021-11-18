using System;
using Burkardt.SubsetNS;

namespace PolPakTest;

public static class combTest
{
    public static void comb_row_next_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMB_ROW_NEXT_TEST tests COMB_ROW_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 10;

        int[] c = new int[N_MAX+1];
        int i;
        int n;

        Console.WriteLine("");
        Console.WriteLine("COMB_ROW_NEXT_TEST");
        Console.WriteLine("  COMB_ROW_NEXT computes the next row of Pascal's triangle.");
        Console.WriteLine("");

        for ( n = 0; n <= N_MAX; n++ )
        {
            Comb.comb_row_next ( n, ref c );
            string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  ";
            for ( i = 0; i <= n; i++ )
            {
                cout += c[i].ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }
            Console.WriteLine(cout);
        }

    }

}