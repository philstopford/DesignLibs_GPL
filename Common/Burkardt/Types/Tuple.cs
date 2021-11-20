using System;

namespace Burkardt.Types;

public class BTupleData
{
    public int[] base_;
}

public static class BTuple
{
    public static void tuple_next_fast(int m, int n, int rank, ref int[] base_, ref int[] x)

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
        //    04 June 2015
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
        //    Input, int RANK, indicates the rank of the tuples.
        //    Typically, 0 <= RANK < N^M; values larger than this are legal
        //    and meaningful, and are equivalent to the corresponding value
        //    MOD N^M.  If RANK < 0, this indicates that this is the first
        //    call for the given values of (M,N).  Initialization is done,
        //    and X is set to a dummy value.
        //
        //    Input/output, int BASE[N], a bookkeeping array needed by 
        //    this procedure.  The user should allocate space for this array, but
        //    should not alter it between successive calls.
        //
        //    Output, int X[N], the next tuple, or a dummy value if initialization
        //    is being done.
        //
    {
        int i;

        switch (rank)
        {
            case < 0 when m <= 0:
                Console.WriteLine("");
                Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                Console.WriteLine("  The value M <= 0 is not legal.");
                Console.WriteLine("  M = " + m + "");
                return;
            case < 0 when n <= 0:
                Console.WriteLine("");
                Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                Console.WriteLine("  The value N <= 0 is not legal.");
                Console.WriteLine("  N = " + n + "");
                return;
            case < 0:
            {
                base_[n - 1] = 1;
                for (i = n - 2; 0 <= i; i--)
                {
                    base_[i] = base_[i + 1] * m;
                }

                for (i = 0; i < n; i++)
                {
                    x[i] = -1;
                }

                break;
            }
            default:
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = rank / base_[i] % m + 1;
                }

                break;
            }
        }
    }

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
        switch (rank)
        {
            //
            case < 0 when m <= 0:
                Console.WriteLine("");
                Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                Console.WriteLine("  The value M <= 0 is not legal.");
                Console.WriteLine("  M = " + m + "");
                return;
            case < 0 when n <= 0:
                Console.WriteLine("");
                Console.WriteLine("TUPLE_NEXT_FAST - Fatal error!");
                Console.WriteLine("  The value N <= 0 is not legal.");
                Console.WriteLine("  N = " + n + "");
                return;
            case < 0:
            {
                data.base_[n - 1] = 1;
                for (i = n - 2; 0 <= i; i--)
                {
                    data.base_[i] = data.base_[i + 1] * m;
                }

                for (i = 0; i < n; i++)
                {
                    x[i] = -1;
                }

                break;
            }
            default:
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = rank / data.base_[i] % m + 1;
                }

                break;
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

        if (m2 < m1)
        {
            rank = 0;
            return;
        }

        switch (rank)
        {
            case <= 0:
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = m1;
                }

                rank = 1;
                break;
            }
            default:
            {
                rank += 1;
                i = n - 1;

                for (;;)
                {

                    if (x[i] < m2)
                    {
                        x[i] += 1;
                        break;
                    }

                    x[i] = m1;

                    if (i == 0)
                    {
                        rank = 0;
                        int j;
                        for (j = 0; j < n; j++)
                        {
                            x[j] = m1;
                        }

                        break;
                    }

                    i -= 1;
                }

                break;
            }
        }
    }

    public static void tuple_next_ge(int m, int n, ref int k, ref int[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT_GE computes the next "nondecreasing" element of a tuple space.
        //
        //  Discussion:
        //
        //    The elements are N vectors.  Each element is constrained to lie
        //    between 1 and M, and to have components that are nondecreasing.
        //    That is, for an element X, and any positive K,
        //      X(I) <= X(I+K)
        //
        //    The elements are produced one at a time.
        //    The first element is
        //      (1,1,...,1)
        //    and the last element is
        //      (M,M,...,M)
        //    Intermediate elements are produced in lexicographic order.
        //
        //  Example:
        //
        //    N = 3, M = 3
        //
        //    K   X
        //    --  -----
        //     1  1 1 1
        //     2  1 1 2
        //     3  1 1 3
        //     4  1 2 2
        //     5  1 2 3
        //     6  1 3 3
        //     7  2 2 2
        //     8  2 2 3
        //     9  2 3 3
        //    10  3 3 3
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
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int M, the maximum entry.
        //
        //    Input, int N, the number of components.
        //
        //    Input/output, int &K, counts the elements.
        //    On first call, set K to 0.  Thereafter, K will indicate the
        //    order of the element returned.  When there are no more elements,
        //    K will be returned as 0.
        //
        //    Input/output, int X[N], on input the previous tuple.
        //    On output, the next tuple.
        //
    {
        int i;
        int j;

        switch (m)
        {
            case < 1:
                return;
        }

        switch (k)
        {
            case <= 0:
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = 1;
                }

                k = 1;
                return;
            }
        }

        for (i = n - 1; 0 <= i; i--)
        {
            if (x[i] >= m)
            {
                continue;
            }

            x[i] += 1;
            for (j = i + 1; j < n; j++)
            {
                x[j] = x[i];
            }

            k += 1;
            return;
        }

        k = 0;
        for (j = 0; j < n; j++)
        {
            x[j] = 0;
        }
    }

    public static void tuple_next2(int n, int[] xmin, int[] xmax, ref int[] x, ref int rank)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUPLE_NEXT2 computes the next element of an integer tuple space.
        //
        //  Discussion:
        //
        //    The elements X are N vectors.
        //
        //    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
        //
        //    The elements are produced one at a time.
        //
        //    The first element is
        //      (XMIN(1), XMIN(2), ..., XMIN(N)),
        //    the second is (probably)
        //      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
        //    and the last element is
        //      (XMAX(1), XMAX(2), ..., XMAX(N))
        //
        //    Intermediate elements are produced in a lexicographic order, with
        //    the first index more important than the last, and the ordering of
        //    values at a fixed index implicitly defined by the sign of
        //    XMAX(I) - XMIN(I).
        //
        //  Example:
        //
        //    N = 2,
        //    XMIN = (/ 1, 10 /)
        //    XMAX = (/ 3,  8 /)
        //
        //    RANK    X
        //    ----  -----
        //      1   1 10
        //      2   1  9
        //      3   1  8
        //      4   2 10
        //      5   2  9
        //      6   2  8
        //      7   3 10
        //      8   3  9
        //      9   3  8
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
        //    Input, int N, the number of components.
        //
        //    Input, int XMIN[N], XMAX[N], the "minimum" and "maximum" entry values.
        //    These values are minimum and maximum only in the sense of the lexicographic
        //    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
        //    than XMAX(I).
        //
        //    Input/output, int X[N], on input the previous tuple.
        //    On output, the next tuple.
        //
        //    Input/output, int &RANK, the rank of the item.  On first call,
        //    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
        //    there are no more items in the sequence.
        //
    {
        int i;

        switch (rank)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("TUPLE_NEXT2 - Fatal error!");
                Console.WriteLine("  Illegal value of RANK = " + rank + "");
                return;
        }

        int test = 1;
        for (i = 0; i < n; i++)
        {
            test *= 1 + Math.Abs(xmax[i] - xmin[i]);
        }

        if (test < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("TUPLE_NEXT2 - Fatal error!");
            Console.WriteLine("  Illegal value of RANK = " + rank + "");
            return;
        }

        switch (rank)
        {
            case 0:
            {
                for (i = 0; i < n; i++)
                {
                    x[i] = xmin[i];
                }

                rank = 1;
                return;
            }
        }

        rank += 1;
        i = n - 1;

        for (;;)
        {
            if (x[i] != xmax[i])
            {
                if (xmin[i] < xmax[i])
                {
                    x[i] += 1;
                }
                else
                {
                    x[i] -= 1;
                }

                break;
            }

            x[i] = xmin[i];

            if (i == 0)
            {
                rank = 0;
                break;
            }

            i -= 1;
        }
    }
}