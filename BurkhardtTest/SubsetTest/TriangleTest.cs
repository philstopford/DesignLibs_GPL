using System;
using Burkardt.Types;

namespace SubsetTestNS;

public static class TriangleTest
{
    public static void subtriangle_next_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBTRIANGLE_NEXT_TEST tests SUBTRIANGLE_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4;
        int rank = 0;

        bool more = false;
        int i1 = 0;
        int j1 = 0;
        int i2 = 0;
        int j2 = 0;
        int i3 = 0;
        int j3 = 0;

        Console.WriteLine("");
        Console.WriteLine("SUBTRIANGLE_NEXT_TEST");
        Console.WriteLine("  SUBTRIANGLE_NEXT generates the indices of subtriangles");
        Console.WriteLine("  in a triangle whose edges were divided into N subedges.");
        Console.WriteLine("");
        Console.WriteLine("  For this test, N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("  Rank    I1  J1    I2  J2    I3  J3");
        Console.WriteLine("");

        for ( ; ; )
        {
            typeMethods.subtriangle_next ( n, ref more, ref i1, ref j1, ref i2, ref j2, ref i3, ref j3 );

            rank += 1;

            Console.WriteLine("  " + rank.ToString().PadLeft(4) + "  "
                              + "  " + i1.ToString().PadLeft(2)
                              + "  " + j1.ToString().PadLeft(2) + "  "
                              + "  " + i2.ToString().PadLeft(2)
                              + "  " + j2.ToString().PadLeft(2) + "  "
                              + "  " + i3.ToString().PadLeft(2)
                              + "  " + j3.ToString().PadLeft(2) + ""); 

            if ( !more )
            {
                break;
            }

        }
    }

}