using System;
using Burkardt;
using Burkardt.Function;

namespace SubsetTestNS
{
    public static class FrobeniusTest
    {
        public static void frobenius_number_order2_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   FROBENIUS_NUMBER_ORDER2_TEST tests FROBENIUS_NUMBER_ORDER2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c1 = 0;
            int c2 = 0;
            int f1 = 0;
            int f2 = 0;
            int n_data = 0;

            Console.WriteLine("");
            Console.WriteLine("FROBENIUS_NUMBER_ORDER2_TEST");
            Console.WriteLine("  FROBENIUS_NUMBER_ORDER2 computes Frobenius numbers of order 2.");
            Console.WriteLine("");
            Console.WriteLine("        C1        C1   exact F  comput F");
            Console.WriteLine("");

            n_data = 0;

            for ( ; ; )
            {
                Burkardt.TestValues.Frobenius.frobenius_number_order2_values ( ref n_data, ref c1, ref c2, ref f1 );

                if ( n_data == 0 )
                {
                    break;
                }

                f2 = Frobenius.frobenius_number_order2 ( c1, c2 );

                Console.WriteLine("  " + c1.ToString().PadLeft(8)
                                  + "  " + c2.ToString().PadLeft(8)
                                  + "  " + f1.ToString().PadLeft(8)
                                  + "  " + f2.ToString().PadLeft(8) + "");
            }
        }

    }
}