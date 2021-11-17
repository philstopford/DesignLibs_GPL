using System;
using Burkardt;
using Burkardt.SortNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SubsetTestNS;

public static class SortHeapExternalTest
{
    public static void sort_heap_external_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SORT_HEAP_EXTERNAL_TEST tests SORT_HEAP_EXTERNAL.
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
        int[] a;
        int i;
        int indx;
        int isgn;
        int j;
        int n = 20;
        int seed;
        int temp;
        SortHeapExternalData data = new();

        Console.WriteLine("");
        Console.WriteLine("SORT_HEAP_EXTERNAL_TEST");
        Console.WriteLine("  SORT_HEAP_EXTERNAL sorts objects externally.");

        seed = 123456789;

        a = UniformRNG.i4vec_uniform_ab_new ( n, 1, n, ref seed );

        typeMethods.i4vec1_print ( n, a, "  Before sorting:" );

        indx = 0;
        i = 0;
        j = 0;
        isgn = 0;

        for ( ; ; )
        {
            Sort.sort_heap_external (ref data, n, ref indx, ref i, ref j, isgn );
 
            if ( indx < 0 )
            {

                if ( a[i-1] <= a[j-1] )
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else if ( 0 < indx )
            {
                temp = a[i-1];
                a[i-1] = a[j-1];
                a[j-1] = temp;
            }
            else
            {
                break;
            }

        }

        typeMethods.i4vec1_print ( n, a, "  After sorting:" ); 
    }

}