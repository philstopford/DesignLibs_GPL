namespace Burkardt.Sequence
{
    public static class Stirling
    {
        public static int[] stirling1(int n, int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STIRLING1 computes the Stirling numbers of the first kind.
            //
            //  Discussion:
            //
            //    The absolute value of the Stirling number S1(N,M) gives the number
            //    of permutations on N objects having exactly M cycles, while the
            //    sign of the Stirling number records the sign (odd or even) of
            //    the permutations.  For example, there are six permutations on 3 objects:
            //
            //      A B C   3 cycles (A) (B) (C)
            //      A C B   2 cycles (A) (BC)
            //      B A C   2 cycles (AB) (C)
            //      B C A   1 cycle  (ABC)
            //      C A B   1 cycle  (ABC)
            //      C B A   2 cycles (AC) (B)
            //
            //    There are
            //
            //      2 permutations with 1 cycle, and S1(3,1) = 2
            //      3 permutations with 2 cycles, and S1(3,2) = -3,
            //      1 permutation with 3 cycles, and S1(3,3) = 1.
            //
            //    Since there are N! permutations of N objects, the sum of the absolute
            //    values of the Stirling numbers in a given row,
            //
            //      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
            //
            //  First terms:
            //
            //    N/M:  1     2      3     4     5    6    7    8
            //
            //    1     1     0      0     0     0    0    0    0
            //    2    -1     1      0     0     0    0    0    0
            //    3     2    -3      1     0     0    0    0    0
            //    4    -6    11     -6     1     0    0    0    0
            //    5    24   -50     35   -10     1    0    0    0
            //    6  -120   274   -225    85   -15    1    0    0
            //    7   720 -1764   1624  -735   175  -21    1    0
            //    8 -5040 13068 -13132  6769 -1960  322  -28    1
            //
            //  Recursion:
            //
            //    S1(N,1) = (-1)^(N-1) * (N-1)! for all N.
            //    S1(I,I) = 1 for all I.
            //    S1(I,J) = 0 if I < J.
            //
            //    S1(N,M) = S1(N-1,M-1) - (N-1) * S1(N-1,M)
            //
            //  Properties:
            //
            //    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
            //
            //    X_N = sum ( 0 <= K <= N ) S1(N,K) X^K
            //    where X_N is the falling factorial function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows of the table.
            //
            //    Input, int M, the number of columns of the table.
            //
            //    Output, int STIRLING1[N*M], the Stirling numbers of the first kind.
            //
        {
            int i;
            int j;
            int[] s1;

            if (n <= 0)
            {
                return null;
            }

            if (m <= 0)
            {
                return null;
            }

            s1 = new int[n * m];

            s1[0 + 0 * n] = 1;
            for (j = 2; j <= m; j++)
            {
                s1[0 + (j - 1) * n] = 0;
            }

            for (i = 2; i <= n; i++)
            {
                s1[i - 1 + 0 * n] = -(i - 1) * s1[i - 2 + 0 * n];
                for (j = 2; j <= m; j++)
                {
                    s1[i - 1 + (j - 1) * n] = s1[i - 2 + (j - 2) * n] - (i - 1) * s1[i - 2 + (j - 1) * n];
                }

            }

            return s1;
        }

        public static int[] stirling2(int n, int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STIRLING2 computes the Stirling numbers of the second kind.
            //
            //  Discussion:
            //
            //    S2(N,M) represents the number of distinct partitions of N elements
            //    into M nonempty sets.  For a fixed N, the sum of the Stirling
            //    numbers S2(N,M) is represented by B(N), called "Bell's number",
            //    and represents the number of distinct partitions of N elements.
            //
            //    For example, with 4 objects, there are:
            //
            //    1 partition into 1 set:
            //
            //      (A,B,C,D)
            //
            //    7 partitions into 2 sets:
            //
            //      (A,B,C) (D)
            //      (A,B,D) (C)
            //      (A,C,D) (B)
            //      (A) (B,C,D)
            //      (A,B) (C,D)
            //      (A,C) (B,D)
            //      (A,D) (B,C)
            //
            //    6 partitions into 3 sets:
            //
            //      (A,B) (C) (D)
            //      (A) (B,C) (D)
            //      (A) (B) (C,D)
            //      (A,C) (B) (D)
            //      (A,D) (B) (C)
            //      (A) (B,D) (C)
            //
            //    1 partition into 4 sets:
            //
            //      (A) (B) (C) (D)
            //
            //    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
            //
            //
            //  First terms:
            //
            //    N/M: 1    2    3    4    5    6    7    8
            //
            //    1    1    0    0    0    0    0    0    0
            //    2    1    1    0    0    0    0    0    0
            //    3    1    3    1    0    0    0    0    0
            //    4    1    7    6    1    0    0    0    0
            //    5    1   15   25   10    1    0    0    0
            //    6    1   31   90   65   15    1    0    0
            //    7    1   63  301  350  140   21    1    0
            //    8    1  127  966 1701 1050  266   28    1
            //
            //  Recursion:
            //
            //    S2(N,1) = 1 for all N.
            //    S2(I,I) = 1 for all I.
            //    S2(I,J) = 0 if I < J.
            //
            //    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
            //
            //  Properties:
            //
            //    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
            //
            //    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
            //    where X_K is the falling factorial function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows of the table.
            //
            //    Input, int M, the number of columns of the table.
            //
            //    Output, int STIRLING2[N*M], the Stirling numbers of the second kind.
            //
        {
            int i;
            int j;
            int[] s2;

            if (n <= 0)
            {
                return null;
            }

            if (m <= 0)
            {
                return null;
            }

            s2 = new int[n * m];

            s2[0 + (0) * n] = 1;
            for (j = 2; j <= m; j++)
            {
                s2[0 + (j - 1) * n] = 0;
            }

            for (i = 2; i <= n; i++)
            {
                s2[i - 1 + (0) * n] = 1;

                for (j = 2; j <= m; j++)
                {
                    s2[i - 1 + (j - 1) * n] = j * s2[i - 2 + (j - 1) * n] + s2[i - 2 + (j - 2) * n];
                }

            }

            return s2;
        }
    }
}