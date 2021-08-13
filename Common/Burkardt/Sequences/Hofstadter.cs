namespace Burkardt.Sequence
{
    public static class Hofstadter
    {
        public static int f_hofstadter(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_HOFSTADTER computes the Hofstadter F sequence.
            //
            //  Discussion:
            //
            //    F(N) = 0                if N = 0
            //         = N - F ( N - 1 ), otherwise.
            //
            //    F(N) is defined for all nonnegative integers, and turns out
            //    to be equal to int ( ( N + 1 ) / 2 ).
            //
            //  Table:
            //
            //     N  F(N)
            //    --  ----
            //
            //     0     0
            //     1     1
            //     2     1
            //     3     2
            //     4     2
            //     5     3
            //     6     3
            //     7     4
            //     8     4
            //     9     5
            //    10     5
            //    11     6
            //    12     6
            //    13     7
            //    14     7
            //    15     8
            //    16     8
            //    17     9
            //    18     9
            //    19    10
            //    20    10
            //    21    11
            //    22    11
            //    23    12
            //    24    12
            //    25    13
            //    26    13
            //    27    14
            //    28    14
            //    29    15
            //    30    15
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Douglas Hofstadter,
            //    Goedel, Escher, Bach,
            //    Basic Books, 1979.
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the function.
            //
            //    Output, int F_HOFSTADTER, the value of the function.
            //
        {
            if (n <= 0)
            {
                return 0;
            }
            else
            {
                return (n - f_hofstadter(n - 1));
            }
        }

        public static int g_hofstadter(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    G_HOFSTADTER computes the Hofstadter G sequence.
            //
            //  Discussion:
            //
            //    G(N) = 0                      if N = 0
            //         = N - G ( G ( N - 1 ) ), otherwise.
            //
            //    G(N) is defined for all nonnegative integers.
            //
            //    The value of G(N) turns out to be related to the Zeckendorf
            //    representation of N as a sum of non-consecutive Fibonacci numbers.
            //    To compute G(N), determine the Zeckendorf representation:
            //
            //      N = sum ( 1 <= I <= M ) F(I)
            //
            //    and reduce the index of each Fibonacci number by 1:
            //
            //      G(N) = sum ( 1 <= I <= M ) F(I-1)
            //
            //    However, this is NOT how the computation is done in this routine.
            //    Instead, a straightforward recursive function call is defined
            //    to correspond to the definition of the mathematical function.
            //
            //  Table:
            //
            //     N  G(N)  Zeckendorf   Decremented
            //    --  ----  ----------   -----------
            //
            //     1   1    1            1
            //     2   1    2            1
            //     3   2    3            2
            //     4   3    3 + 1        2 + 1
            //     5   3    5            3
            //     6   4    5 + 1        3 + 1
            //     7   4    5 + 2        3 + 1
            //     8   5    8            5
            //     9   6    8 + 1        5 + 1
            //    10   6    8 + 2        5 + 1
            //    11   7    8 + 3        5 + 2
            //    12   8    8 + 3 + 1    5 + 2 + 1
            //    13   8    13           8
            //    14   9    13 + 1       8 + 1
            //    15   9    13 + 2       8 + 1
            //    16  10    13 + 3       8 + 2
            //    17  11    13 + 3 + 1   8 + 2 + 1
            //    18  11    13 + 5       8 + 3
            //    19  12    13 + 5 + 1   8 + 3 + 1
            //    20  12    13 + 5 + 2   8 + 3 + 1
            //    21  13    21           13
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Douglas Hofstadter,
            //    Goedel, Escher, Bach,
            //    Basic Books, 1979.
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the function.
            //
            //    Output, int G_HOFSTADTER, the value of the function.
            //
        {
            if (n <= 0)
            {
                return 0;
            }
            else
            {
                return (n - g_hofstadter(g_hofstadter(n - 1)));
            }

        }
        
        public static int h_hofstadter ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    H_HOFSTADTER computes the Hofstadter H sequence.
            //
            //  Discussion:
            //
            //    H(N) = 0                          if N = 0
            //         = N - H ( H ( H ( N - 1 ) ), otherwise.
            //
            //    H(N) is defined for all nonnegative integers.
            //
            //  Table:
            //
            //     N  H(N)
            //    --  ----
            //
            //     0     0
            //     1     1
            //     2     1
            //     3     2
            //     4     3
            //     5     4
            //     6     4
            //     7     5
            //     8     5
            //     9     6
            //    10     7
            //    11     7
            //    12     8
            //    13     9
            //    14    10
            //    15    10
            //    16    11
            //    17    12
            //    18    13
            //    19    13
            //    20    14
            //    21    14
            //    22    15
            //    23    16
            //    24    17
            //    25    17
            //    26    18
            //    27    18
            //    28    19
            //    29    20
            //    30    20
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Douglas Hofstadter,
            //    Goedel, Escher, Bach,
            //    Basic Books, 1979.
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the function.
            //
            //    Output, int H_HOFSTADTER, the value of the function.
            //
        {
            if ( n <= 0 )
            {
                return 0;
            }
            else
            {
                return ( n - h_hofstadter ( h_hofstadter ( h_hofstadter ( n-1 ) ) ) );
            }

        }

        public static int v_hofstadter ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_HOFSTADTER computes the Hofstadter V sequence.
            //
            //  Discussion:
            //
            //    V(N) = 0                          if N = 0
            //         = 1                          if 1 <= N <= 4
            //         = V(N-V(N-1)) + V(N-V(N-4)), otherwise.
            //
            //    V(N) is defined for all nonnegative integers.
            //
            //  Table:
            //
            //     N  V(N)
            //    --  ----
            //
            //     0     0
            //     1     1
            //     2     1
            //     3     1
            //     4     1
            //     5     2
            //     6     3
            //     7     4
            //     8     5
            //     9     5
            //    10     6
            //    11     6
            //    12     7
            //    13     8
            //    14     8
            //    15     9
            //    16     9
            //    17    10
            //    18    11
            //    19    11
            //    20    11
            //    21    12
            //    22    12
            //    23    13
            //    24    14
            //    25    14
            //    26    15
            //    27    15
            //    28    16
            //    29    17
            //    30    17
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the function.
            //
            //    Output, int V_HOFSTADTER, the value of the function.
            //
        {
            if ( n <= 0 )
            {
                return 0;
            }
            else if ( n <= 4 )
            {
                return 1;
            }
            else
            {
                return (  v_hofstadter ( n - v_hofstadter(n-1) ) 
                          + v_hofstadter ( n - v_hofstadter(n-4) ) );
            }

        }

    }
}