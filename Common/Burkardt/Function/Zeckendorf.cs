using Burkardt.Sequence;

namespace Burkardt.Function;

public class Zeckendorf
{
    public static void zeckendorf(int n, int m_max, ref int m, ref int[] i_list, ref int[] f_list )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZECKENDORF produces the Zeckendorf decomposition of a positive integer.
        //
        //  Discussion:
        //
        //    Zeckendorf proved that every positive integer can be represented
        //    uniquely as the sum of non-consecutive Fibonacci numbers.
        //
        //    N = sum ( 1 <= I <= M ) F_LIST(I)
        //
        //  Example:
        //
        //     N    Decomposition
        //
        //    50    34 + 13 + 3
        //    51    34 + 13 + 3 + 1
        //    52    34 + 13 + 5
        //    53    34 + 13 + 5 + 1
        //    54    34 + 13 + 5 + 2
        //    55    55
        //    56    55 + 1
        //    57    55 + 2
        //    58    55 + 3
        //    59    55 + 3 + 1
        //    60    55 + 5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2000
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the positive integer to be decomposed.
        //
        //    Input, int M_MAX, the maximum number of parts in the decomposition.
        //    Set M_MAX = 100 to be safe.
        //
        //    Output, int M, the number of parts in the decomposition.
        //
        //    Output, int I_LIST[M], the index of the Fibonacci numbers
        //    in the decomposition.
        //
        //    Output, int F_LIST[M], the value of the Fibonacci numbers
        //    in the decomposition.
        //
    {
        int f = 0;
        int i = 0;
        int j = 0;

        m = 0;
        //
        //  Extract a sequence of Fibonacci numbers.
        //
        while (0 < n && m < m_max)
        {
            Fibonacci.fibonacci_floor(n, ref f, ref i);

            i_list[m] = i;
            m += 1;
            n -= f;
        }

        //
        //  Replace any pair of consecutive indices ( I, I-1 ) by I+1.
        //
        for (i = m; 2 <= i; i--)
        {
            if (i_list[i - 2] != i_list[i - 1] + 1)
            {
                continue;
            }

            i_list[i - 2] += 1;
            for (j = i; j <= m - 1; j++)
            {
                i_list[j - 1] = i_list[j];
            }

            i_list[m - 1] = 0;
            m -= 1;

        }

        //
        //  Fill in the actual values of the Fibonacci numbers.
        //
        for (i = 0; i < m; i++)
        {
            f_list[i] = Fibonacci.fibonacci_direct(i_list[i]);
        }

    }

}