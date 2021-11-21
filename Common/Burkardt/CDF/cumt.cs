namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumt ( double t, double df, ref double cum, ref double ccum )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMT evaluates the cumulative T distribution.
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
        //    Formula 26.5.27.
        //
        //  Parameters:
        //
        //    Input, double *T, the upper limit of integration.
        //
        //    Input, double *DF, the number of degrees of freedom of
        //    the T distribution.
        //
        //    Output, double *CUM, *CCUM, the T distribution CDF and
        //    complementary CDF.
        //
    {
        double a = 0;
        const double K2 = 0.5e0;
        double oma = 0;

        double tt = t * t;
        double dfptt = df + tt;
        double xx = df / dfptt;
        double yy = tt / dfptt;
        double T1 = 0.5e0 * df;
        cumbet ( xx, yy, T1, K2, ref a, ref oma );

        switch (t)
        {
            case <= 0.0e0:
                cum = 0.5e0 * a;
                ccum = oma + cum;
                break;
            default:
                ccum = 0.5e0 * a;
                cum = oma + ccum;
                break;
        }
    }
}