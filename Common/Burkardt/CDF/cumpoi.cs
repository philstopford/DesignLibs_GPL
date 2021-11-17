namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumpoi(double s, double xlam, ref double cum, ref double ccum)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMPOI evaluates the cumulative Poisson distribution.
        //
        //  Discussion:
        //
        //    CUMPOI returns the probability of S or fewer events in a Poisson
        //    distribution with mean XLAM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    Formula 26.4.21.
        //
        //  Parameters:
        //
        //    Input, double *S, the upper limit of cumulation of the
        //    Poisson density function.
        //
        //    Input, double *XLAM, the mean of the Poisson distribution.
        //
        //    Output, double *CUM, *CCUM, the Poisson density CDF and
        //    complementary CDF.
        //
    {
        double df = 2.0e0 * (s + 1.0e0);
        double chi = 2.0e0 * xlam;
        cumchi(chi, df, ref ccum, ref cum);
    }
}