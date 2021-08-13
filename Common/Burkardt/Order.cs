using System;

namespace Burkardt.OrderNS
{
    public static class Order
    {
        public static bool order_check ( int order )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ORDER_CHECK checks the value of ORDER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ORDER, the requested order.
            //
            //    Output, bool ORDER_CHECK, is TRUE if the requested order is acceptable.
            //
        {
            bool check;

            check = ( order ==   1 ) ||
                    ( order ==   3 ) ||
                    ( order ==   7 ) ||
                    ( order ==  15 ) ||
                    ( order ==  31 ) ||
                    ( order ==  63 ) ||
                    ( order == 127 ) ||
                    ( order == 255 ) ||
                    ( order == 511 );

            return check;
        }
        
        public static int order_code ( string code )

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    ORDER_CODE returns the order for each element.
            //
            //  List:
            //
            //    CODE  Order  Definition
            //    ----  -----  ----------
            //    Q4     4     4 node linear Lagrange/serendipity quadrilateral;
            //    Q8     8     8 node quadratic serendipity quadrilateral;
            //    Q9     9     9 node quadratic Lagrange quadrilateral;
            //    Q12   12     12 node cubic serendipity quadrilateral;
            //    Q16   16     16 node cubic Lagrange quadrilateral;
            //    QL     6     6 node linear/quadratic quadrilateral;
            //    T3     3     3 node linear triangle;
            //    T4     4     4 node bubble triangle;
            //    T6     6     6 node quadratic triangle;
            //    T10   10     10 node cubic triangle.
            // 
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 March 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string CODE, the code for the element.
            //
            //    Output, int ORDER_CODE, the order of the element.
            //
        {
            int value;

            if ( code == "Q4" )
            {
                value = 4;
            }
            else if ( code == "Q8" )
            {
                value = 8;
            }
            else if ( code == "Q9" )
            {
                value = 9;
            }
            else if ( code == "Q12" )
            {
                value = 12;
            }
            else if ( code == "Q16" )
            {
                value = 16;
            }
            else if ( code == "QL" )
            {
                value = 6;
            }
            else if ( code == "T3" )
            {
                value = 3;
            }
            else if ( code == "T4" )
            {
                value = 4;
            }
            else if ( code == "T6" )
            {
                value = 6;
            }
            else if ( code == "T10" )
            {
                value = 10;
            }
            else
            {
                value = -1;
            }

            return value;
        }
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