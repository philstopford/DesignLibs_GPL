using System;
using Burkardt.SimplexNS;

namespace PolPakTest;

public static class simplexTest
{
    public static void simplex_num_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_NUM_TEST tests SIMPLEX_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int m;
        int n;
        int value;

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_NUM_TEST");
        Console.WriteLine("  SIMPLEX_NUM computes the N-th simplex number");
        Console.WriteLine("  in M dimensions.");
        Console.WriteLine("");
        Console.WriteLine("      M: 0     1     2     3     4     5");
        Console.WriteLine("   N");

        for (n = 0; n <= 10; n++)
        {
            string cout = "  " + n.ToString().PadLeft(2);
            for (m = 0; m <= 5; m++)
            {
                value = Simplex.simplex_num(m, n);
                cout += "  " + value.ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

}