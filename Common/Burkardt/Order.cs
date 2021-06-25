using System;

namespace Burkardt
{
    public static class Order
    {
        public static int order_from_level_135 ( int l )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
            //
            //  Discussion:
            //
            //    Clenshaw Curtis rules, and some others, often use the following
            //    scheme:
            //
            //    L: 0  1  2  3   4   5
            //    N: 1  3  5  9  17  33 ... 2^L+1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int L, the level, which should be 0 or greater.
            //
            //    Output, int ORDER_FROM_LEVEL_135, the order.
            //
        {
            int n = 0;

            if ( l < 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("ORDER_FROM_LEVEL_135 - Fatal error!");
                Console.WriteLine("  Illegal input value of L!");
                return ( 1 );
            }
            else if ( l == 0 )
            {
                n = 1;
            }
            else
            {
                n = (int)Math.Pow ( 2, l ) + 1;
            }
            return n;
        }
    }
}