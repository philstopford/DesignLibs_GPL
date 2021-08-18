using System;
using Burkardt.PolynomialNS;

namespace PolPakTest
{
    public static class bellTest
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
            //    02 June 2007
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
                Burkardt.Values.Bell.bell_values ( ref n_data, ref n, ref c );

                if ( n_data == 0 )
                {
                    break;
                }

                c2 = new int[n+1];

                Burkardt.Values.Bell.bell ( n, ref c2 );

                Console.WriteLine("  "
                                         + n.ToString().PadLeft(4)     + "  "
                                         + c.ToString().PadLeft(8)     + "  "
                                         + c2[n].ToString().PadLeft(8) + "");
                
            }

        }

        public static void bell_poly_coef_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BELL_POLY_COEF_TEST tests BELL_POLY_COEF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int []c;
            int i;
            int n;
            int n_max = 10;

            Console.WriteLine("");
            Console.WriteLine("BELL_POLY_COEF_TEST");
            Console.WriteLine("  BELL_POLY_COEF returns the coefficients of a");
            Console.WriteLine("  Bell polynomial.");
            Console.WriteLine("");
            Console.WriteLine("  Table of polynomial coefficients:");
            Console.WriteLine("");

            for ( n = 0; n <= n_max; n++ )
            {
                c = Bell.bell_poly_coef ( n );
                string cout = "  " + n.ToString().PadLeft(2) + ":  ";
                for ( i = 0; i <= n; i++ )
                {
                    cout += "  " + c[i].ToString().PadLeft(6);
                }
                Console.WriteLine(cout);
            }

        }

    }
}