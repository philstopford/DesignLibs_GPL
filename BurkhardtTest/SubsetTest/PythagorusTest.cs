using System;
using Burkardt;

namespace SubsetTest
{
    public static class PythagorusTest
    {
        public static void pythag_triple_next_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYTHAG_TRIPLE_NEXT_TEST tests PYTHAG_TRIPLE_NEXT;
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int a = 0;
            int b = 0;
            int c = 0;
            int d = 0;
            int e = 0;
            int i = 0;
            int j = 0;
            int k = 0;

            Console.WriteLine("");
            Console.WriteLine("PYTHAG_TRIPLE_NEXT_TEST");
            Console.WriteLine("  PYTHAG_TRIPLE_NEXT computes the next");
            Console.WriteLine("  Pythagorean triple.");
            Console.WriteLine("");
            Console.WriteLine("     I       J       A       B       C A^2+B^2     C^2");
            Console.WriteLine("");

            i = 0;
            j = 0;

            for ( k = 0; k <= 20; k++ )
            {
                Pythagorus.pythag_triple_next ( ref i, ref j, ref a, ref b, ref c );

                d = a * a + b * b;
                e = c * c;

                Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                  + j.ToString().PadLeft(6) + "  "
                                  + a.ToString().PadLeft(6) + "  "
                                  + b.ToString().PadLeft(6) + "  "
                                  + c.ToString().PadLeft(6) + "  "
                                  + d.ToString().PadLeft(6) + "  "
                                  + e.ToString().PadLeft(6) + "");
            }
        }
    }
}