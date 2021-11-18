using System;
using Burkardt;

namespace SubsetTestNS;

public static class EulerTest
{
    public static void euler_row_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EULER_ROW_TEST tests EULER_ROW.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 9;

        int i;
        int[] ieuler = new int[N_MAX+1];
        int n;

        Console.WriteLine("");
        Console.WriteLine("EULER_ROW_TEST");
        Console.WriteLine("  EULER_ROW gets rows of the Euler triangle.");
        Console.WriteLine("");

        for ( n = 0; n <= N_MAX; n++ )
        {
            Permutation.euler_row ( n, ref ieuler );

            string cout = "";
            for ( i = 0; i <= n; i++ )
            {
                cout += ieuler[i].ToString(CultureInfo.InvariantCulture).PadLeft(7) + "  ";
            }
            Console.WriteLine(cout);
        }
 
    }
}