﻿using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class CatalanTest
{
    public static void catalan_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CATALAN_VALUES_TEST tests CATALAN_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int n = 0;
        Console.WriteLine("");
        Console.WriteLine("CATALAN_VALUES_TEST:");
        Console.WriteLine("  CATALAN_VALUES returns values of ");
        Console.WriteLine("  the Catalan numbers.");
        Console.WriteLine("");
        Console.WriteLine("     N        C(N)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            Catalan.catalan_values ( ref n_data, ref n, ref c );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + c.ToString().PadLeft(10) + "");
        }
    }

}