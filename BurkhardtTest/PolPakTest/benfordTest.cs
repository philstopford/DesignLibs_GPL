using System;
using Burkardt.Probability;

namespace PolPakTest;

public static class benfordTest
{
    public static void benford_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BENFORD_TEST tests BENFORD.
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

        Console.WriteLine("");
        Console.WriteLine("BENFORD_TEST");
        Console.WriteLine("  BENFORD(I) is the Benford probability of the");
        Console.WriteLine("  initial digit sequence I.");
        Console.WriteLine("");
        Console.WriteLine("     I  BENFORD(I)");
        Console.WriteLine("");

        for ( i = 1; i <= 9; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)              + "  "
                              + Benford.benford_pdf ( i ).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

}