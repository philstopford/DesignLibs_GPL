using System;

namespace Burkardt.Types
{
    public class BTupleData
    {
        public int[] base_;
    }

    public static class BTuple
    {
        public static void tuple_next_fast(ref BTupleData data, int m, int n, int rank, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
            //
            //  Discussion:
            //
            //    The elements are N vectors.  Each entry is constrained to lie
            //    between 1 and M.  The elements are produced one at a time.
            //    The first element is
            //      (1,1,...,1)
            //    and the last element is
            //      (M,M,...,M)
            //    Intermediate elements are produced in lexicographic order.
            //
            //  Example:
            //
            //    N = 2,
            //    M = 3
            //
            //    INPUT        OUTPUT
            //    -------      -------
            //    Rank          X
            //    ----          ----
            //   -1            -1 -1
            //
            //    0             1  1
            //    1             1  2
            //    2             1  3
            //    3             2  1
            //    4             2  2
            //    5             2  3
            //    6             3  1
            //    7             3  2
            //    8             3  3
            //    9             1  1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the maximum entry in each component.
            //    M must be greater than 0.
            //
            //    Input, int N, the number of components.
            //    N must be greater than 0.
            //
            //    Input, integer RANK, indicates the rank of the tuples.
            //    Typically, 0 <= RANK < N**M; values larger than this are legal
            //    and meaningful, and are equivalent to the corresponding value
            //    MOD N**M.  If RANK < 0, this indicates that this is the first
            //    call for the given values of (M,N).  Initialization is done,
            //    and X is set to a dummy value.
            //
            //    Output, int X[N], the next tuple, or a dummy value if initialization
            //    is being done.
            //
        {
            int i;
            //
            if (rank < 0)
            {
                if (m <= 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                    Console.WriteLine("  The value M <= 0 is not legal.");
                    Console.WriteLine("  M = " + m + "");
                    return;
                }

                if (n <= 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                    Console.WriteLine("  The value N <= 0 is not legal.");
                    Console.WriteLine("  N = " + n + "");
                    return;
                }

                data.base_[n - 1] = 1;
                for (i = n - 2; 0 <= i; i--)
                {
                    data.base_[i] = data.base_[i + 1] * m;
                }

                for (i = 0; i < n; i++)
                {
                    x[i] = -1;
                }
            }
            else
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = ((rank / data.base_[i]) % m) + 1;
                }
            }
        }

        public static void tuple_next(int m1, int m2, int n, ref int rank, ref int[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TUPLE_NEXT computes the next element of a tuple space.
            //
            //  Discussion:
            //
            //    The elements are N vectors.  Each entry is constrained to lie
            //    between M1 and M2.  The elements are produced one at a time.
            //    The first element is
            //      (M1,M1,...,M1),
            //    the second element is
            //      (M1,M1,...,M1+1),
            //    and the last element is
            //      (M2,M2,...,M2)
            //    Intermediate elements are produced in lexicographic order.
            //
            //  Example:
            //
            //    N = 2, M1 = 1, M2 = 3
            //
            //    INPUT        OUTPUT
            //    -------      -------
            //    Rank  X      Rank   X
            //    ----  ---    -----  ---
            //    0     * *    1      1 1
            //    1     1 1    2      1 2
            //    2     1 2    3      1 3
            //    3     1 3    4      2 1
            //    4     2 1    5      2 2
            //    5     2 2    6      2 3
            //    6     2 3    7      3 1
            //    7     3 1    8      3 2
            //    8     3 2    9      3 3
            //    9     3 3    0      0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M1, M2, the minimum and maximum entries.
            //
            //    Input, int N, the number of components.
            //
            //    Input/output, int *RANK, counts the elements.
            //    On first call, set RANK to 0.  Thereafter, the output value of RANK
            //    will indicate the order of the element returned.  When there are no
            //    more elements, RANK will be returned as 0.
            //
            //    Input/output, int X[N], on input the previous tuple.
            //    On output, the next tuple.
            //
        {
            int i;
            int j;

            if (m2 < m1)
            {
                rank = 0;
                return;
            }

            if (rank <= 0)
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = m1;
                }

                rank = 1;
            }
            else
            {
                rank = rank + 1;
                i = n - 1;

                for (;;)
                {

                    if (x[i] < m2)
                    {
                        x[i] = x[i] + 1;
                        break;
                    }

                    x[i] = m1;

                    if (i == 0)
                    {
                        rank = 0;
                        for (j = 0; j < n; j++)
                        {
                            x[j] = m1;
                        }

                        break;
                    }

                    i = i - 1;
                }
            }
        }
    }
}