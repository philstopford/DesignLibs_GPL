namespace Burkardt.Ziggurat;

public static class SHR3
{
    public static int shr3_seeded ( ref int jsr )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHR3_SEEDED evaluates the SHR3 generator for integers.
        //
        //  Discussion:
        //
        //    Thanks to Dirk Eddelbuettel for pointing out that this code needed to
        //    use the uint32_t data type in order to execute properly in 64 bit mode,
        //    03 October 2013.
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
        //    Input/output, uint32_t &JSR, the seed, which is updated 
        //    on each call.
        //
        //    Output, uint32_t SHR3_SEEDED, the new value.
        //
    {
        int jsr_input = jsr;

        jsr ^= jsr <<   13;
        jsr ^= jsr >>   17;
        jsr ^= jsr <<    5;

        int value = jsr_input + jsr;

        return value;
    }
}