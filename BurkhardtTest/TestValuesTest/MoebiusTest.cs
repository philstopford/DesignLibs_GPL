using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class MoebiusTest
{

    public static void moebius_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    MOEBIUS_VALUES_TEST tests MOEBIUS_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fn = 0;
        int n = 0;
        Console.WriteLine("");
        Console.WriteLine("MOEBIUS_VALUES_TEST:");
        Console.WriteLine("  MOEBIUS_VALUES returns values of");
        Console.WriteLine("  the Moebius function.");
        Console.WriteLine("");
        Console.WriteLine("     N         MU(N)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Moebius.moebius_values(ref n_data, ref n, ref fn);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(8) + "  "
                              + fn.ToString().PadLeft(12) + "");
        }
    }
}