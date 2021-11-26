using System;
using Burkardt.Types;

namespace LegendreProductPolynomialTest;

public static class permTest
{
    public static void perm_uniform_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_UNIFORM_TEST tests PERM_UNIFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 10;
        int test;

        Console.WriteLine("");
        Console.WriteLine("PERM_UNIFORM_TEST");
        Console.WriteLine("  PERM_UNIFORM randomly selects a permutation.");
        Console.WriteLine("");

        int seed = 123456789;

        for ( test = 1; test <= 5; test++ )
        {
            int[] p = typeMethods.perm_uniform_new ( n, ref seed );
            string cout = "  ";
            int i;
            for ( i = 0; i < n; i++ )
            {
                cout += p[i].ToString().PadLeft(4);
            }
            Console.WriteLine(cout);
        }
    }
}