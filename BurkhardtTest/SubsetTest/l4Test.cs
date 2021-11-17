using System;
using Burkardt.Types;

namespace SubsetTestNS;

public static class l4Test
{
    public static void l4vec_next_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L4VEC_NEXT_TEST tests L4VEC_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool done;
        int i;
        bool[] l4vec = new bool[3];
        int n = 3;

        Console.WriteLine("");
        Console.WriteLine("L4VEC_NEXT_TEST");
        Console.WriteLine("  L4VEC_NEXT generates logical vectors in order.");
        Console.WriteLine("");

        for ( i = 0; i < n; i++ )
        {
            l4vec[i] = false;
        }
 
        done = false;

        for ( ; ; )
        {
            string cout = "  ";
            for ( i = 0; i < n; i++ )
            {
                cout += l4vec[i];
            }
            Console.WriteLine(cout);

            if ( done )
            {
                break;
            }

            typeMethods.l4vec_next ( n, ref l4vec );

            done = true;
            for ( i = 0; i < n; i++ )
            {
                done = l4vec[i] switch
                {
                    false => false,
                    _ => done
                };
            }
        }
    }

}