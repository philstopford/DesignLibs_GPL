namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumf(double f, double dfn, double dfd, ref double cum, ref double ccum)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMF evaluates the cumulative F distribution.
        //
        //  Discussion:
        //
        //    CUMF computes the integral from 0 to F of the F density with DFN
        //    numerator and DFD denominator degrees of freedom.
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
        //    Handbook of Mathematical Functions
        //    1966, Formula 26.5.28.
        //
        //  Parameters:
        //
        //    Input, double *F, the upper limit of integration.
        //
        //    Input, double *DFN, *DFD, the number of degrees of
        //    freedom for the numerator and denominator.
        //
        //    Output, double *CUM, *CCUM, the value of the F CDF and
        //    the complementary F CDF.
        //
    {
        const double half = 0.5e0;
        const double done = 1.0e0;

        int ierr = 0;
        double yy;

        switch (f)
        {
            case <= 0.0e0:
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
        }

        double prod = dfn * f;
        //
        //  XX is such that the incomplete beta with parameters
        //  DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
        //  YY is 1 - XX
        //  Calculate the smaller of XX and YY accurately
        //
        double dsum = dfd + prod;
        double xx = dfd / dsum;

        if (xx > half)
        {
            yy = prod / dsum;
            xx = done - yy;
        }
        else
        {
            yy = done - xx;
        }

        double T1 = dfd * half;
        double T2 = dfn * half;
        beta_inc(T1, T2, xx, yy, ref ccum, ref cum, ref ierr);
    }
}