using System;
using Burkardt;
using Burkardt.Sequence;

namespace SubsetTestNS
{
    public static class SchroederTest
    {
        public static void schroeder_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCHROEDER_TEST tests SCHROEDER;
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
            int N = 10;

            int i;
            int[] s = new int[N];

            Console.WriteLine("");
            Console.WriteLine("SCHROEDER_TEST");
            Console.WriteLine("  SCHROEDER computes the Schroeder numbers.");

            Schroeder.schroeder ( N, ref s );

            Console.WriteLine("");
            Console.WriteLine("   N    S(N)");
            Console.WriteLine("");

            for ( i = 0; i < N; i++ )
            {
                Console.WriteLine("  "
                                        + (i+1).ToString().PadLeft(4)  + "  "
                                        + s[i].ToString().PadLeft(6) + "");
            }
        }

    }
}