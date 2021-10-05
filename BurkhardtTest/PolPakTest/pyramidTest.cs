using System;
using Burkardt.PyramidNS;

namespace PolPakTest
{
    public static class pyramidTest
    {
        public static void pyramid_num_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_NUM_TEST tests PYRAMID_NUM.
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
            int n;

            Console.WriteLine("");
            Console.WriteLine("PYRAMID_NUM_TEST");
            Console.WriteLine("  PYRAMID_NUM computes the pyramidal numbers.");
            Console.WriteLine("");

            for (n = 1; n <= 10; n++)
            {
                Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + Pyramid.pyramid_num(n).ToString().PadLeft(6) + "");
            }

        }

        public static void pyramid_square_num_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PYRAMID_SQUARE_NUM_TEST tests PYRAMID_SQUARE_NUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;

            Console.WriteLine("");
            Console.WriteLine("PYRAMID_SQUARE_NUM_TEST");
            Console.WriteLine("  PYRAMID_SQUARE_NUM computes the pyramidal square numbers.");
            Console.WriteLine("");

            for (n = 1; n <= 10; n++)
            {
                Console.WriteLine("  "
                              + n.ToString().PadLeft(6) + "  "
                              + Pyramid.pyramid_square_num(n).ToString().PadLeft(6) + "");
            }

        }

    }
}