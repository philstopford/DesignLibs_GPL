using System;
using Burkardt.TestValues;

namespace SubsetTestNS
{
    public static class BellTest
    {
        public static void bell_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BELL_TEST tests BELL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c = 0;
            int[] c2;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("BELL_TEST");
            Console.WriteLine("  BELL computes Bell numbers.");
            Console.WriteLine("");
            Console.WriteLine("  N  exact C(I)  computed C(I)");
            Console.WriteLine("");

            n_data = 0;

            for ( ; ; )
            {
                Bell.bell_values ( ref n_data, ref n, ref c );

                if ( n_data == 0 )
                {
                    break;
                }

                c2 = new int[n+1];

                Bell.bell ( n, ref c2 );

                Console.WriteLine("  "
                                     + n.ToString().PadLeft(4) + "  "
                                     + c.ToString().PadLeft(8) + "  "
                                     + c2[n].ToString().PadLeft(8) + "");

            }
        }
    }
}