using System;
using Burkardt;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class TriangTest
{
    public static void triang_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANG_TEST tests TRIANG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 10;

        int[] a = {
            1,0,1,0,1,0,1,0,0,1, 
            0,1,0,0,1,0,0,0,0,0, 
            0,0,1,0,1,0,1,0,0,1, 
            0,1,1,1,1,1,1,1,0,1, 
            0,0,0,0,1,0,0,0,0,0, 
            0,1,0,0,1,1,1,0,0,0, 
            0,0,0,0,1,0,1,0,0,0, 
            0,1,0,0,1,1,1,1,0,1, 
            0,0,0,0,0,0,0,0,0,0, 
            0,0,0,0,1,0,1,0,0,1 };
        int[] p = new int[N];

        Console.WriteLine("");
        Console.WriteLine("TRIANG_TEST");
        Console.WriteLine("  TRIANG relabels elements for a partial ordering,");

        typeMethods.i4mat_print ( N, N, a, "  The input matrix:" );
 
        Triang.triang ( N, a, ref p );
 
        Permutation.perm0_print ( N, p, "  The new ordering:" );

        typeMethods.i4mat_2perm0 ( N, N, a, p, p );
 
        typeMethods.i4mat_print ( N, N, a, "  The reordered matrix:" );
            
    }

}