namespace SubsetTest
{
    public static class DebruijnTest
    {
        public static void debruijn_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //     DEBRUIJN_TEST tests DEBRUIJN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
# define NUM_TEST 3

            int i;
            int ihi;
            int m;
            int mtest[NUM_TEST] = { 2, 3, 2 };
            int n;
            int ntest[NUM_TEST] = { 3, 3, 4 };
            int string[28];
            int test;

            Console.WriteLine("");
            Console.WriteLine("DEBRUIJN_TEST");
            Console.WriteLine("  DEBRUIJN computes a de Bruijn string.");

            for ( test = 0; test < NUM_TEST; test++ )
            {
                m = mtest[test];
                n = ntest[test];

                Console.WriteLine("");
                Console.WriteLine("  The alphabet size is M = " + m + "");
                Console.WriteLine("  The string length is N = " + n + "");

                debruijn ( m, n, string );

                ihi = i4_power ( m, n );

                Console.WriteLine("");
                Console.WriteLine("  ";
                for ( i = 0; i < ihi; i++ )
                {
                    Console.WriteLine(setw(1) + string[i];
                }
                Console.WriteLine("");

            }

            return;
            UM_TEST
        }

    }
}