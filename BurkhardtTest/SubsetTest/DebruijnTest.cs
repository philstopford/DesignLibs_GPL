﻿using System;
using Burkardt.Sequence;

namespace SubsetTestNS;

public static class DebruijnTest
{
    public static void debruijn_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //     DEBRUIJN_TEST tests DEBRUIJN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NUM_TEST = 3;

        int[] mtest = { 2, 3, 2 };
        int[] ntest = { 3, 3, 4 };
        int[] string_ = new int[28];
        int test;

        Console.WriteLine("");
        Console.WriteLine("DEBRUIJN_TEST");
        Console.WriteLine("  DEBRUIJN computes a de Bruijn string.");

        for ( test = 0; test < NUM_TEST; test++ )
        {
            int m = mtest[test];
            int n = ntest[test];

            Console.WriteLine("");
            Console.WriteLine("  The alphabet size is M = " + m + "");
            Console.WriteLine("  The string length is N = " + n + "");

            Debruijn.debruijn ( m, n, ref string_ );

            int ihi = (int)Math.Pow ( m, n );

            Console.WriteLine("");
            string cout = "  ";
            int i;
            for ( i = 0; i < ihi; i++ )
            {
                cout += string_[i];
            }
            Console.WriteLine(cout);

        }
    }

}