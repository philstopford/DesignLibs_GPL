using System;

namespace Burkardt.OrderNS;

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
        bool check = order is 1 or 3 or 7 or 15 or 31 or 63 or 127 or 255 or 511;

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
        int value = code switch
        {
            "Q4" => 4,
            "Q8" => 8,
            "Q9" => 9,
            "Q12" => 12,
            "Q16" => 16,
            "QL" => 6,
            "T3" => 3,
            "T4" => 4,
            "T6" => 6,
            "T10" => 10,
            _ => -1
        };

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

        switch (l)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("ORDER_FROM_LEVEL_135 - Fatal error!");
                Console.WriteLine("  Illegal input value of L!");
                return 1;
            case 0:
                n = 1;
                break;
            default:
                n = (int)Math.Pow ( 2, l ) + 1;
                break;
        }
        return n;
    }
}