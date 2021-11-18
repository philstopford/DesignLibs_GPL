using System;
using Burkardt.Function;

namespace PolPakTest;

public static class primeTest
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
        int n;
        int prime_max;

        Console.WriteLine("");
        Console.WriteLine("PRIME_TEST");
        Console.WriteLine("  PRIME returns primes from a table.");

        n = -1;
        prime_max = Prime.prime ( n );
        Console.WriteLine("");
        Console.WriteLine("  Number of primes stored is " + prime_max + "");
        Console.WriteLine("");
        Console.WriteLine("     I    Prime(I)");
        Console.WriteLine("");
        for ( i = 1; i <= 10; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + Prime.prime ( i ).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
        Console.WriteLine("");
        for ( i = prime_max - 10; i <= prime_max; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + Prime.prime ( i ).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

    }

    public static void phi_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PHI_TEST tests PHI.
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
        int c = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("PHI_TEST");
        Console.WriteLine("  PHI computes the PHI function.");
        Console.WriteLine("");
        Console.WriteLine("  N   Exact   PHI(N)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Phi.phi_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + c.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + Prime.phi(n).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }

    }
}