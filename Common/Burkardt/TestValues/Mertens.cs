namespace TestValues
{
    public static class Mertens
    {
        
        public static void mertens_values ( ref int n_data, ref int n, ref int c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MERTENS_VALUES returns some values of the Mertens function.
            //
            //  Discussion:
            //
            //    The Mertens function M(N) is the sum from 1 to N of the Moebius
            //    function MU.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 October 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    M Deleglise, J Rivat,
            //    Computing the Summation of the Moebius Function,
            //    Experimental Mathematics,
            //    Volume 5, 1996, pages 291-295.
            //
            //    Eric Weisstein,
            //    CRC Concise Encyclopedia of Mathematics,
            //    CRC Press, 2002,
            //    Second edition,
            //    ISBN: 1584883472,
            //    LC: QA5.W45.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and N_DATA
            //    is set to 1.  On each subsequent call, the input value of N_DATA is
            //    incremented and that test data item is returned, if available.  When
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, ref int N, the argument of the Mertens function.
            //
            //    Output, ref int C, the value of the Mertens function.
            //
        {
            int N_MAX = 15;

            int[] c_vec = {
                1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1,
                -2,  -2,   1,    2, -23 };
            int[] n_vec = {
                1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
                11,  12,  100, 1000, 10000 };

            if ( n_data < 0 )
            {
                n_data = 0;
            }

            if ( N_MAX <= n_data )
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