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
        int i1;
        int i2;
        int i3;
        int j1;
        int j2;
        int j3;
        bool more;
        int n;
        int rank;

        n = 4;
        rank = 0;

        more = false;
        i1 = 0;
        j1 = 0;
        i2 = 0;
        j2 = 0;
        i3 = 0;
        j3 = 0;

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