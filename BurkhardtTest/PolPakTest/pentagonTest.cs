using System;
using Burkardt.Function;

namespace PolPakTest;

public static class pentagonTest
{
    public static void pentagon_num_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PENTAGON_NUM_TEST tests PENTAGON_NUM.
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
        int n;

        Console.WriteLine("");
        Console.WriteLine("PENTAGON_NUM_TEST");
        Console.WriteLine("  PENTAGON_NUM computes the pentagonal numbers.");
        Console.WriteLine("");

        for ( n = 1; n <= 10; n++ )
        {
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + Pentagon.pentagon_num ( n ).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

    }
}