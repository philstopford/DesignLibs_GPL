namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumbin ( double s, double xn, double pr, double ompr,
            ref double cum, ref double ccum )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMBIN evaluates the cumulative binomial distribution.
        //
        //  Discussion:
        //
        //    This routine returns the probability of 0 to S successes in XN binomial
        //    trials, each of which has a probability of success, PR.
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
        //    1966, Formula 26.5.24.
        //
        //  Parameters:
        //
        //    Input, double *S, the upper limit of summation.
        //
        //    Input, double *XN, the number of trials.
        //
        //    Input, double *PR, the probability of success in one trial.
        //
        //    Input, double *OMPR, equals ( 1 - PR ).
        //
        //    Output, double *CUM, the cumulative binomial distribution.
        //
        //    Output, double *CCUM, the complement of the cumulative
        //    binomial distribution.
        //
    {
        if ( s < xn )
        {
            double T1 = s + 1.0;
            double T2 = xn - s;
            cumbet ( pr, ompr, T1, T2, ref ccum, ref cum );
        }
        else
        {
            cum = 1.0;
            ccum = 0.0;
        }
    }
}