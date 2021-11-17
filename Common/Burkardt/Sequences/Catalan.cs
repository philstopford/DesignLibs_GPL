namespace Burkardt.Sequence;

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

        switch (n)
        {
            case < 0:
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
            c[i] = c[i - 1] * 2 * (2 * i - 1) / (i + 1);
        }

        return c;
    }

    public static void catalan_row_next(bool next, int n, ref int[] irow)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CATALAN_ROW computes row N of Catalan's triangle.
        //
        //  Example:
        //
        //    I\J 0   1   2   3   4   5   6
        //
        //    0   1
        //    1   1   1
        //    2   1   2   2
        //    3   1   3   5   5
        //    4   1   4   9  14  14
        //    5   1   5  14  28  42  42
        //    6   1   6  20  48  90 132 132
        //
        //  Recursion:
        //
        //    C(0,0) = 1
        //    C(I,0) = 1
        //    C(I,J) = 0 for I < J
        //    C(I,J) = C(I,J-1) + C(I-1,J)
        //    C(I,I) is the I-th Catalan number.
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
        //  Parameters:
        //
        //    Input, bool NEXT, indicates whether this is a call for
        //    the 'next' row of the triangle.
        //    NEXT = FALSE, this is a startup call.  Row N is desired, but
        //    presumably this is a first call, or row N-1 was not computed
        //    on the previous call.
        //    NEXT = TRUE, this is not the first call, and row N-1 was computed
        //    on the previous call.  In this case, much work can be saved
        //    by using the information from the previous values of IROW
        //    to build the next values.
        //
        //    Input, int N, the index of the row of the triangle desired.
        //
        //    Input/output, int IROW[N+1], the row of coefficients.
        //    If NEXT = FALSE, then IROW is not required to be set on input.
        //    If NEXT = TRUE, then IROW must be set on input to the value of
        //    row N-1.
        //
    {
        int i;
        int j;
        switch (n)
        {
            //
            case < 0:
                return;
        }

        switch (next)
        {
            case false:
            {
                irow[0] = 1;
                for (i = 1; i <= n; i++)
                {
                    irow[i] = 0;
                }

                for (i = 1; i <= n; i++)
                {
                    irow[0] = 1;

                    for (j = 1; j <= i - 1; j++)
                    {
                        irow[j] += irow[j - 1];
                    }

                    irow[i] = irow[i - 1];

                }

                break;
            }
            default:
            {
                irow[0] = 1;

                for (j = 1; j <= n - 1; j++)
                {
                    irow[j] += irow[j - 1];
                }

                irow[n] = n switch
                {
                    >= 1 => irow[n - 1],
                    _ => irow[n]
                };

                break;
            }
        }
    }
        
}