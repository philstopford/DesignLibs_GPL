using System;
using Burkardt.Sequence;

namespace PolPakTest;

public static class collatzTest
{
    public static void collatz_count_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COLLATZ_COUNT_TEST tests COLLATZ_COUNT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int count = 0;
        int count2 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("COLLATZ_COUNT_TEST:");
        Console.WriteLine("  COLLATZ_COUNT(N) counts the length of the");
        Console.WriteLine("  Collatz sequence beginning with N.");
        Console.WriteLine("");
        Console.WriteLine("       N       COUNT(N)     COUNT(N)");
        Console.WriteLine("              (computed)    (table)");
        Console.WriteLine("");

        n_data = 0;

        for ( ; ; )
        {
            Burkardt.Values.Collatz.collatz_count_values ( ref n_data, ref n, ref count );

            if ( n_data == 0 )
            {
                break;
            }

            count2 = Collatz.collatz_count ( n );

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + count.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + count2.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

    }

    public static void collatz_count_max_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COLLATZ_COUNT_MAX_TEST tests COLLATZ_COUNT_MAX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i_max = 0;
        int j_max = 0;
        int n;

        Console.WriteLine("");
        Console.WriteLine("COLLATZ_COUNT_MAX_TEST:");
        Console.WriteLine("  COLLATZ_COUNT_MAX(N) returns the length of the");
        Console.WriteLine("  longest Collatz sequence from 1 to N.");
        Console.WriteLine("");
        Console.WriteLine("         N     I_MAX     J_MAX");
        Console.WriteLine("");

        n = 10;

        while ( n <= 100000 )
        {
            Collatz.collatz_count_max ( n, ref i_max, ref j_max );

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + i_max.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + j_max.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

            n *= 10;
        }

    }

}