namespace Burkardt.FEM
{
    public class Order
    {
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
    }
}