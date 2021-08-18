namespace Burkardt.Values
{
    public static class Prime
    {
        public static void prime_values(ref int n_data, ref int n, ref int p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PRIME_VALUES returns values of the prime function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Prime[n]
            //
            //    Thanks to Morten Welinder for pointing out that the index of 145253029
            //    is 8192000, 12 April 2013.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //    Output, ref int N, the index of the prime.
            //
            //    Output, ref int P, the value of the prime.
            //
        {
            int N_MAX = 24;

            int[] n_vec =
            {
                1,
                2,
                4,
                8,
                16,
                32,
                64,
                128,
                256,
                512,
                1000,
                2000,
                4000,
                8000,
                16000,
                32000,
                64000,
                128000,
                256000,
                512000,
                1024000,
                2048000,
                4096000,
                8192000
            };

            int[] p_vec =
            {
                2,
                3,
                7,
                19,
                53,
                131,
                311,
                719,
                1619,
                3671,
                7919,
                17389,
                37813,
                81799,
                176081,
                376127,
                800573,
                1698077,
                3588941,
                7559173,
                15881419,
                33283031,
                69600977,
                145253029
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                p = 0;
            }
            else
            {
                n = n_vec[n_data - 1];
                p = p_vec[n_data - 1];
            }
        }
    }
}