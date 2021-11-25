using System;
using Burkardt.Function;

namespace SubsetTestNS;

public static class PrimeTest
{
    public static void prime_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRIME_TEST tests PRIME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("PRIME_TEST");
        Console.WriteLine("  PRIME returns primes from a table.");

        const int n = -1;
        int prime_max = Prime.prime ( n );
        Console.WriteLine("");
        Console.WriteLine("  Number of primes stored is " + prime_max + "");
        Console.WriteLine("");
        Console.WriteLine("     I    Prime(I)");
        Console.WriteLine("");
        for ( i = 1; i <= 10; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(4) + "  "
                              + Prime.prime ( i ).ToString().PadLeft(6) + "");
        }
        Console.WriteLine("");
        for ( i = prime_max - 10; i <= prime_max; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString().PadLeft(4) + "  "
                              + Prime.prime ( i ).ToString().PadLeft(6) + "");
        }
    }
}