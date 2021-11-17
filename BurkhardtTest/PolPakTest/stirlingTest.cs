using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class stirlingTest
{
    public static void stirling1_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING1_TEST tests STIRLING1.
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
        int i;
        int j;
        int m = 8;
        int n = 8;
        int[] s1;

        Console.WriteLine("");
        Console.WriteLine("STIRLING1_TEST");
        Console.WriteLine("  STIRLING1: Stirling numbers of first kind.");
        Console.WriteLine("  Get rows 1 through " + m + "");
        Console.WriteLine("");

        s1 = Stirling.stirling1 ( m, n );

        for ( i = 0; i < m; i++ )
        {
            string cout = (i+1).ToString().PadLeft(6) + "  ";
            for ( j = 0; j < n; j++ )
            {
                cout += s1[i+j*m].ToString().PadLeft(6) + "  ";
            }
            Console.WriteLine(cout);
        }

    }

    public static void stirling2_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING2_TEST tests STIRLING2.
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
        int  M = 8;
        int  N = 8;

        int i;
        int j;
        int[] s2;

        Console.WriteLine("");
        Console.WriteLine("STIRLING2_TEST");
        Console.WriteLine("  STIRLING2: Stirling numbers of second kind.");
        Console.WriteLine("  Get rows 1 through " + M + "");
        Console.WriteLine("");

        s2 = Stirling.stirling2 ( M, N );

        for ( i = 0; i < M; i++ )
        {
            string cout = (i+1).ToString().PadLeft(6) + "  ";
            for ( j = 0; j < N; j++ )
            {
                cout += s2[i+j*M].ToString().PadLeft(6) + "  ";
            }
            Console.WriteLine(cout);
        }

    }

}