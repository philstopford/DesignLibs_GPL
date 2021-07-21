namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int r8vec_order_type(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries of the array.
            //
            //    Input, double X[N], the array to be checked.
            //
            //    Output, int R8VEC_ORDER_TYPE, order indicator:
            //    -1, no discernable order;
            //    0, all entries are equal;
            //    1, ascending order;
            //    2, strictly ascending order;
            //    3, descending order;
            //    4, strictly descending order.
            //
        {
            int i;
            int order;
            //
            //  Search for the first value not equal to X(0).
            //
            i = 0;

            for (;;)
            {
                i = i + 1;
                if (n - 1 < i)
                {
                    order = 0;
                    return order;
                }

                if (x[0] < x[i])
                {
                    if (i == 1)
                    {
                        order = 2;
                        break;
                    }
                    else
                    {
                        order = 1;
                        break;
                    }
                }
                else if (x[i] < x[0])
                {
                    if (i == 1)
                    {
                        order = 4;
                        break;
                    }
                    else
                    {
                        order = 3;
                        break;
                    }
                }
            }

            //
            //  Now we have a "direction".  Examine subsequent entries.
            //
            for (;;)
            {
                i = i + 1;
                if (n - 1 < i)
                {
                    break;
                }

                if (order == 1)
                {
                    if (x[i] < x[i - 1])
                    {
                        order = -1;
                        break;
                    }
                }
                else if (order == 2)
                {
                    if (x[i] < x[i - 1])
                    {
                        order = -1;
                        break;
                    }
                    else if (x[i] == x[i - 1])
                    {
                        order = 1;
                    }
                }
                else if (order == 3)
                {
                    if (x[i - 1] < x[i])
                    {
                        order = -1;
                        break;
                    }
                }
                else if (order == 4)
                {
                    if (x[i - 1] < x[i])
                    {
                        order = -1;
                        break;
                    }
                    else if (x[i] == x[i - 1])
                    {
                        order = 3;
                    }
                }
            }

            return order;
        }

    }
}