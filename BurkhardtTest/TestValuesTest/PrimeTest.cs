using System;
using Burkardt.Values;

namespace TestValuesTest;

public class PrimeTest
{
    public static void prime_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRIME_VALUES_TEST tests PRIME_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 0;
        int n_data;
        int p = 0;
        Console.WriteLine("");
        Console.WriteLine("PRIME_VALUES_TEST:");
        Console.WriteLine("  PRIME_VALUES returns values of");
        Console.WriteLine("  the prime function.");
        Console.WriteLine("");
        Console.WriteLine("           N          P[N]");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Prime.prime_values ( ref n_data, ref n, ref p );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + n.ToString().PadLeft(12) + "  "
                              + p.ToString().PadLeft(12) + "");
        }
    }
}