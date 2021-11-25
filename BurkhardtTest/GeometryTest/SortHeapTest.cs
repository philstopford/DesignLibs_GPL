using System;
using Burkardt.SortNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class SortHeapTest
{
    public static void test180 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST180 tests SORT_HEAP_EXTERNAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 20;

        int[] a = new int[N];
        SortHeapExternalData data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST180");
        Console.WriteLine("  SORT_HEAP_EXTERNAL sorts objects externally.");

        int indx = 0;
        int i;
        int j = 0;
        int isgn = 0;
        int seed = 123456789;

        for ( i = 0; i < N; i++ )
        {
            a[i] = UniformRNG.i4_uniform_ab ( 1, N, ref seed );
        }

        typeMethods.i4vec_print ( N, a, "  Unsorted array:" );
        //
        //void the sort routine over and over.
        //
        for ( ;; )
        {
            Sort.sort_heap_external ( ref data, N, ref indx, ref i, ref j, isgn );
            //
            //  If the return value of INDX is negative, we're asked to compare
            //  array elements I and J;
            //
            if ( indx < 0 )
            {
                if ( a[i-1] <= a[j-1] )
                {
                    isgn = -1;
                }
                else
                {
                    isgn = 1;
                }
            }
            //
            //  ...and if the return value of INDX is positive, we're asked to switch
            //  array elements I and J;
            //
            else if ( 0 < indx )
            {
                (a[i - 1], a[j - 1]) = (a[j - 1], a[i - 1]);
                //
                //  ...and if the return value of INDX is 0, we're done.
                //
            }
            else
            {
                break;
            }
        }

        typeMethods.i4vec_print ( N, a, "  Sorted array:" );
    }

}