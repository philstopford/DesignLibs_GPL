using System;
using Burkardt;
using Burkardt.Function;

namespace SubsetTestNS
{
    public static class PentEnumTest
    {
        public static void pent_enum_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PENT_ENUM_TEST tests PENT_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int i;

            Console.WriteLine("");
            Console.WriteLine("PENT_ENUM_TEST");
            Console.WriteLine("  PENT_ENUM counts points in pentagons.");
            Console.WriteLine("");
            Console.WriteLine("   N  Pent(N)");
            Console.WriteLine("");

            for ( i = 0; i <= N; i++ )
            {
                Console.WriteLine(i.ToString().PadLeft(4)               + "  "
                                  + Pentagon.pent_enum ( i ).ToString().PadLeft(6) + "");
            }
        }
    }
}