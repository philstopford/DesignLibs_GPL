namespace Burkardt.Treepack
{
    public static class Catalan
    {
        public static int[] catalan(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CATALAN computes the Catalan numbers, from C(0) to C(N).
            //
            //  Discussion:
            //
            //    The Catalan number C(N) counts:
            //
            //    1) the number of binary trees on N vertices;
            //    2) the number of ordered trees on N+1 vertices;
            //    3) the number of full binary trees on 2N+1 vertices;
            //    4) the number of well formed sequences of 2N parentheses;
            //    5) the number of ways 2N ballots can be counted, in order,
            //       with N positive and N negative, so that the running sum
            //       is never negative;
            //    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
            //    7) the number of monotone functions from [1..N} to [1..N} which
            //       satisfy f(i) <= i for all i;
            //    8) the number of ways to triangulate a polygon with N+2 vertices.
            //
            //    The formula is:
            //
            //      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
            //           = 1 / (N+1) * COMB ( 2N, N )
            //           = 1 / (2N+1) * COMB ( 2N+1, N+1).
            //
            //  First values:
            //
            //     C(0)     1
            //     C(1)     1
            //     C(2)     2
            //     C(3)     5
            //     C(4)    14
            //     C(5)    42
            //     C(6)   132
            //     C(7)   429
            //     C(8)  1430
            //     C(9)  4862
            //    C(10) 16796
            //
            //  Recursion:
            //
            //    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
            //    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
            //
            //  Example:
            //
            //    N = 3
            //
            //    ()()()
            //    ()(())
            //    (()())
            //    (())()
            //    ((()))
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dennis Stanton, Dennis White,
            //    Constructive Combinatorics,
            //    Springer, 1986,
            //    ISBN: 0387963472,
            //    LC: QA164.S79.
            //
            //  Parameters:
            //
            //    Input, int N, the number of Catalan numbers desired.
            //
            //    Output, int CATALAN[N+1], the Catalan numbers from C(0) to C(N).
            //
        {
            int[] c;
            int i;

            if (n < 0)
            {
                return null;
            }

            c = new int[n + 1];

            c[0] = 1;
            //
            //  The extra parentheses ensure that the integer division is
            //  done AFTER the integer multiplication.
            //
            for (i = 1; i <= n; i++)
            {
                c[i] = (c[i - 1] * 2 * (2 * i - 1)) / (i + 1);
            }

            return c;
        }

        public static void catalan_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CATALAN_VALUES returns some values of the Catalan numbers for testing.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    US Department of Commerce, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input/output, int &N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and N_DATA
            //    is set to 1.  On each subsequent call, the input value of N_DATA is
            //    incremented and that test data item is returned, if available.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, int &N, the order of the Catalan number.
            //
            //    Output, int &C, the value of the Catalan number.
            //
        {
            int N_MAX = 11;

            int[] c_vec = {1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796};
            int[] n_vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

            if (n_data < 0)
            {
                n_data = 0;
            }

            if (N_MAX <= n_data)
            {
                n_data = 0;
                n = 0;
                c = 0;
            }
            else
            {
                n = n_vec[n_data];
                c = c_vec[n_data];
                n_data = n_data + 1;
            }

        }
    }
}