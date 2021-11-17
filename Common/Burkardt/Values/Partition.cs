namespace Burkardt.Values;

public static class Partition
{

    public static void partition_count_values(ref int n_data, ref int n, ref int c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_COUNT_VALUES returns some values of the int *partition count.
        //
        //  Discussion:
        //
        //    A partition of an int *N is a representation of the integer
        //    as the sum of nonzero positive integers.  The order of the summands
        //    does not matter.  The number of partitions of N is symbolized
        //    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
        //    following partitions:
        //
        //    5 = 5
        //      = 4 + 1
        //      = 3 + 2
        //      = 3 + 1 + 1
        //      = 2 + 2 + 1
        //      = 2 + 1 + 1 + 1
        //      = 1 + 1 + 1 + 1 + 1
        //
        //    In Mathematica, the function can be evaluated by
        //
        //      PartitionsP[n]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 February 2003
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
        //    Output, ref int N, the integer.
        //
        //    Output, ref int C, the number of partitions of the integer.
        //
    {
        const int N_MAX = 21;

        int[] c_vec =
        {
            1,
            1, 2, 3, 5, 7, 11, 15, 22, 30, 42,
            56, 77, 101, 135, 176, 231, 297, 385, 490, 627
        };

        int[] n_vec =
        {
            0,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20
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
            c = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            c = c_vec[n_data - 1];
        }
    }

    public static void partition_distinct_count_values(ref int n_data, ref int n, ref int c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
        //
        //  Discussion:
        //
        //    A partition of an int *N is a representation of the integer
        //    as the sum of nonzero positive integers.  The order of the summands
        //    does not matter.  The number of partitions of N is symbolized
        //    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
        //    following partitions:
        //
        //    5 = 5
        //      = 4 + 1
        //      = 3 + 2
        //      = 3 + 1 + 1
        //      = 2 + 2 + 1
        //      = 2 + 1 + 1 + 1
        //      = 1 + 1 + 1 + 1 + 1
        //
        //    However, if we require that each member of the partition
        //    be distinct, so that no nonzero summand occurs more than once,
        //    we are computing something symbolized by Q(N).
        //    The number 5 has Q(N) = 3, because it has the following partitions
        //    into distinct parts:
        //
        //    5 = 5
        //      = 4 + 1
        //      = 3 + 2
        //
        //    In Mathematica, the function can be evaluated by
        //
        //      PartitionsQ[n]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 2003
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
        //    Output, ref int N, the integer.
        //
        //    Output, ref int C, the number of partitions of the integer
        //    into distinct parts.
        //
    {
        const int N_MAX = 21;

        int[] c_vec =
        {
            1,
            1, 1, 2, 2, 3, 4, 5, 6, 8, 10,
            12, 15, 18, 22, 27, 32, 38, 46, 54, 64
        };

        int[] n_vec =
        {
            0,
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20
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
            c = 0;
        }
        else
        {
            n = n_vec[n_data - 1];
            c = c_vec[n_data - 1];
        }
    }

}