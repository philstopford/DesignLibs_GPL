namespace Burkardt.Ziggurat;

public static class MultiplyWithCarry
{
    public static int mwc_seeded ( ref int w, ref int z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MWC_SEEDED evaluates the MWC multiply-with-carry random number generator.
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
        //    Input/output, uint32_t &W, uint32_t &Z, the seeds, which are updated 
        //    on each call.
        //
        //    Output, uint32_t MWC_SEEDED, the new value.
        //
    {
        z = 36969 * ( z & 65535 ) + ( z >> 16 );
        w = 18000 * ( w & 65535 ) + ( w >> 16 );

        int value = ( z << 16 ) + w;

        return value;
    }
}