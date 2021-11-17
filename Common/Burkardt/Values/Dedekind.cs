﻿namespace Burkardt.Values;

public static class Dedekind
{
    public static void dedekind_sum_values(ref int n_data, ref int p, ref int q, ref int n, ref int d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEDEKIND_SUM_VALUES returns some values of the Dedekind sum.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Hans Rademacher, Emil Grosswald,
        //    Dedekind Sums,
        //    Mathematics Association of America, 1972,
        //    LC: QA241.R2.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
        //    before the first call.  On each call, the routine increments N_DATA
        //    by 1, and returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref int P, &Q, the arguments of the function.
        //
        //    Output, ref int N, &D, the numerator and denominator of the
        //    function value.
        //
    {
        const int N_MAX = 95;

        int[] d_vec =
        {
            1, 1, 18, 8, 5, 18, 14, 16, 27, 5,
            22, 72, 13, 14, 90, 32, 17, 27, 38, 40,
            1, 18, 1, 14, 27, 22, 13, 18, 17, 38,
            1, 1, 8, 1, 14, 16, 1, 22, 13, 14,
            32, 17, 38, 8, 1, 18, 5, 14, 27, 22,
            13, 90, 1, 38, 1, 1, 18, 8, 18, 14,
            16, 27, 22, 72, 1, 14, 32, 17, 27, 38,
            1, 5, 14, 22, 13, 17, 38, 1, 1, 18,
            8, 1, 18, 16, 27, 1, 22, 72, 13, 18,
            32, 17, 27, 38, 8
        };
        int[] n_vec =
        {
            0, 0, 1, 1, 1, 5, 5, 7, 14, 3,
            15, 55, 11, 13, 91, 35, 20, 34, 51, 57,
            0, -1, 0, 1, 4, 5, 4, 7, 8, 21,
            0, 0, -1, 0, -1, 1, 0, 3, 1, 3,
            5, 5, 9, 3, 0, 1, -1, 1, -4, 3,
            -1, 19, 0, 11, 0, 0, -1, 1, -5, -1,
            -1, 4, -5, -1, 0, 3, -5, 1, 2, 11,
            0, 1, -5, 5, -4, 5, -9, 0, 0, 1,
            -1, 0, 5, -7, -4, 0, -3, 1, 4, -7,
            -3, 1, -2, 3, 3
        };
        int[] p_vec =
        {
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
            4, 4, 4, 4, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            6, 6, 6, 6, 6, 6, 6, 7, 7, 7,
            7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
            7, 7, 7, 7, 7
        };
        int[] q_vec =
        {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            1, 3, 5, 7, 9, 11, 13, 15, 17, 19,
            1, 2, 4, 5, 7, 8, 10, 11, 13, 14,
            16, 17, 19, 20, 1, 3, 5, 7, 9, 11,
            13, 15, 17, 19, 1, 2, 3, 4, 6, 7,
            8, 9, 11, 12, 13, 14, 16, 17, 18, 19,
            1, 5, 7, 11, 13, 17, 19, 1, 2, 3,
            4, 5, 6, 8, 9, 10, 11, 12, 13, 15,
            16, 17, 18, 19, 20
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
            p = 0;
            q = 0;
            n = 0;
            d = 0;
        }
        else
        {
            p = p_vec[n_data - 1];
            q = q_vec[n_data - 1];
            n = n_vec[n_data - 1];
            d = d_vec[n_data - 1];
        }
    }
}