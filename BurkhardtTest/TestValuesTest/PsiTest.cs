﻿using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class PsiTest
{
    public static void psi_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    PSI_VALUES_TEST tests PSI_VALUES.
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
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("PSI_VALUES_TEST");
        Console.WriteLine("  PSI_VALUES stores values of");
        Console.WriteLine("  the PSI function.");
        Console.WriteLine("");
        Console.WriteLine("      X            PSI(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Psi.psi_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}