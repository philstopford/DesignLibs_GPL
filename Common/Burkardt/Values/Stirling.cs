namespace Burkardt.Values;

public static class Stirling
{

    public static void stirling1_values(ref int n_data, ref int n, ref int m, ref int fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING1_VALUES returns some values of the Stirling numbers, kind 1.
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
        //    In Mathematica, the function can be evaluated by:
        //
        //      StirlingS1[n,m]
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
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N, &M, the arguments of the function.
        //
        //    Output, ref int FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        int[] fx_vec =
        {
            0,
            1,
            -3,
            11,
            -50,
            274,
            -1764,
            13068,
            -109584,
            1026576,
            -13132,
            6769,
            -1960,
            322,
            -28,
            1
        };

        int[] m_vec =
        {
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8
        };

        int[] n_vec =
        {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            m = 0;
            fx = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            m = m_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static void stirling2_values(ref int n_data, ref int n, ref int m, ref int fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING2_VALUES returns some values of the Stirling numbers, kind 2.
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
        //    In Mathematica, the function can be evaluated by:
        //
        //      StirlingS2[n,m]
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
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int N, &M, the arguments of the function.
        //
        //    Output, ref int FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        int[] fx_vec =
        {
            0,
            1,
            3,
            7,
            15,
            31,
            63,
            127,
            255,
            511,
            966,
            1701,
            1050,
            266,
            28,
            1
        };

        int[] m_vec =
        {
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8
        };

        int[] n_vec =
        {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            m = 0;
            fx = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            m = m_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

}