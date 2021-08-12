using System;
using Burkardt;

namespace SubsetTest
{
    public static class MorseThueTest
    {
        public static void morse_thue_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MORSE_THUE_TEST tests MORSE_THUE.
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
            int N = 100;

            int i;
            int ihi;
            int ilo;
            int[] s = new int[N+1];

            Console.WriteLine("");
            Console.WriteLine("MORSE_THUE_TEST");
            Console.WriteLine("  MORSE_THUE computes the Morse-Thue numbers.");
            Console.WriteLine("");

            for ( i = 0; i <= N; i++ )
            {
                s[i] = MorseThue.morse_thue ( i );
            }

            for ( ilo = 0; ilo <= N; ilo = ilo + 10 )
            {
                string cout = "  ";
                ihi = Math.Min ( ilo + 9, N );
                for ( i = ilo; i <= ihi; i++ )
                {
                    cout += s[i].ToString().PadLeft(1);
                }
                Console.WriteLine(cout);
            }
            
        }

    }
}