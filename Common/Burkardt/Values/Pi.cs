namespace Burkardt.Values
{
    public static class Pi
    {

        public static void pi_values(ref int n_data, ref int n, ref int p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PI_VALUES returns values of the Pi function.
            //
            //  Discussion:
            //
            //    Pi[n] is the number of primes less than or equal to n.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      PrimePi[n]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2004
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
            //    Output, ref int N, the argument.
            //
            //    Output, ref int P, the value of the function.
            //
        {
            int N_MAX = 17;

            int[] n_vec =
            {
                10,
                20,
                30,
                40,
                50,
                60,
                70,
                80,
                90,
                100,
                1000,
                10000,
                100000,
                1000000,
                10000000,
                100000000,
                1000000000
            };

            int[] p_vec =
            {
                4,
                8,
                10,
                12,
                15,
                17,
                19,
                22,
                24,
                25,
                168,
                1229,
                9592,
                78498,
                664579,
                5761455,
                50847534
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