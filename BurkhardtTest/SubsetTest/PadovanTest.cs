using System;
using Burkardt.Sequence;

namespace SubsetTestNS;

public class PadovanTest
{
    public static void padovan_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PADOVAN_TEST tests PADOVAN;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 10;

        int i;
        int[] p = new int[N];

        Console.WriteLine("");
        Console.WriteLine("PADOVAN_TEST");
        Console.WriteLine("  PADOVAN computes the Padovan numbers.");
        Console.WriteLine("");

        Padovan.padovan ( N, ref p );

        Console.WriteLine("");
        Console.WriteLine("   N    P(N)");
        Console.WriteLine("");

        for ( i = 0; i < N; i++ )
        {
            Console.WriteLine("  "
                              + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)    + "  "
                              + p[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
            
    }
}