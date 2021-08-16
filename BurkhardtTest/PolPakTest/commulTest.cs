using System;
using Burkardt.PolynomialNS;

namespace PolPakTest
{
    public static class commulTest
    {
        public static void commul_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMMUL_TEST tests COMMUL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] factor = new int[4];
            int i;
            int ncomb;
            int nfactor;

            Console.WriteLine("");
            Console.WriteLine("COMMUL_TEST");
            Console.WriteLine("  COMMUL computes a multinomial coefficient.");
            Console.WriteLine("");

            n = 8;
            nfactor = 2;
            factor[0] = 6;
            factor[1] = 2;
            ncomb = Multinomial.commul ( n, nfactor, factor );
            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of factors = " + factor + "");
            for ( i = 0; i < nfactor; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                  + "  " + factor[i].ToString().PadLeft(8) + "");
            }
            Console.WriteLine("  Value of coefficient = " + ncomb + "");

            n = 8;
            nfactor = 3;
            factor[0] = 2;
            factor[1] = 2;
            factor[2] = 4;
            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of factors = " + factor + "");
            for ( i = 0; i < nfactor; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                  + "  " + factor[i].ToString().PadLeft(8) + "");
            }
            Console.WriteLine("  Value of coefficient = " + ncomb + "");

            n = 13;
            nfactor = 4;
            factor[0] = 5;
            factor[1] = 3;
            factor[2] = 3;
            factor[3] = 2;
            ncomb = Multinomial.commul ( n, nfactor, factor );
            Console.WriteLine("");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  Number of factors = " + factor + "");
            for ( i = 0; i < nfactor; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                  + "  " + factor[i].ToString().PadLeft(8) + "");
            }
            Console.WriteLine("  Value of coefficient = " + ncomb + "");

        }

    }
}