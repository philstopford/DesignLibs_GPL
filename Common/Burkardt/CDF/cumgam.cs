namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumgam ( double x, double a, ref double cum, ref double ccum )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMGAM evaluates the cumulative incomplete gamma distribution.
        //
        //  Discussion:
        //
        //    This routine computes the cumulative distribution function of the
        //    incomplete gamma distribution, i.e., the integral from 0 to X of
        //
        //      (1/GAM(A))*EXP(-T)*T^(A-1) DT
        //
        //    where GAM(A) is the complete gamma function of A, i.e.,
        //
        //      GAM(A) = integral from 0 to infinity of EXP(-T)*T^(A-1) DT
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
        //  Parameters:
        //
        //    Input, double *X, the upper limit of integration.
        //
        //    Input, double *A, the shape parameter of the incomplete
        //    Gamma distribution.
        //
        //    Output, double *CUM, *CCUM, the incomplete Gamma CDF and
        //    complementary CDF.
        //
    {
        int K1 = 0;

        switch (x)
        {
            case <= 0.0e0:
                cum = 0.0e0;
                ccum = 1.0e0;
                break;
            default:
                gamma_inc ( a, x, ref cum, ref ccum, K1 );
                break;
        }
    }
}