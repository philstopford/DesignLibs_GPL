namespace Burkardt.Ziggurat;

public static class KISS
{
    public static int kiss_seeded ( ref int jcong, ref int jsr, ref int w, ref int z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KISS_SEEDED evaluates the KISS random number generator.
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
        //    Input/output, uint32_t &JCONG, uint32_t &JSR, uint32_t &W, uint32_t &Z, 
        //    the seeds, which are updated on each call.
        //
        //    Output, uint32_t KISS_SEEDED, the new value.
        //
    {
        int value = ( MultiplyWithCarry.mwc_seeded ( ref w, ref z ) ^ Congruential.cong_seeded ( ref jcong ) ) + SHR3.shr3_seeded ( ref jsr );

        return value;
    }
}