using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class hailTest
    {
        public static void hail_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAIL_TEST tests HAIL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine("HAIL_TEST");
            Console.WriteLine("  HAIL(I) computes the length of the hail sequence");
            Console.WriteLine("  for I, also known as the 3*N+1 sequence.");
            Console.WriteLine("");
            Console.WriteLine("  I,  HAIL(I)");
            Console.WriteLine("");

            for ( i = 1; i <= 20; i++ )
            {
                Console.WriteLine("  "
                                              + i.ToString().PadLeft(4)          + "  "
                                              + Hail.hail ( i ).ToString().PadLeft(6) + "");
            }

        }
    }
}