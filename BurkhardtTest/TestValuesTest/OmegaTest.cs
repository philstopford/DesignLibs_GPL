﻿using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class OmegaTest
{
    public static void omega_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    OMEGA_VALUES_TEST tests OMEGA_VALUES.
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
        int fn = 0;
        int n = 0;
        Console.WriteLine("");
        Console.WriteLine("OMEGA_VALUES_TEST:");
        Console.WriteLine("  OMEGA_VALUES returns values of");
        Console.WriteLine("  the Omega function.");
        Console.WriteLine("");
        Console.WriteLine("     N           OMEGA(N)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Omega.omega_values(ref n_data, ref n, ref fn);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(12) + n + "  "
                              + fn.ToString().PadLeft(10) + fn + "");
        }
    }
}