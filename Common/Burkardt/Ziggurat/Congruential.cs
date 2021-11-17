namespace Burkardt.Ziggurat;

public static class Congruential
{
    public static int cong_seeded ( ref int jcong )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONG_SEEDED evaluates the CONG congruential random number generator.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Marsaglia, Wai Wan Tsang,
        //    The Ziggurat Method for Generating Random Variables,
        //    Journal of Statistical Software,
        //    Volume 5, Number 8, October 2000, seven pages.
        //
        //  Parameters:
        //
        //    Input/output, uint32_t &JCONG, the seed, which is updated 
        //    on each call.
        //
        //    Output, uint32_t CONG_SEEDED, the new value.
        //
    {
        jcong = 69069 * jcong + 1234567;

        int value = jcong;

        return value;
    }
}