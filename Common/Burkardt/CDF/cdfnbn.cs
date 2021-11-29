using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdfnbn(int which, ref double p, ref double q, ref double s, ref double xn,
            ref double pr, ref double ompr, ref int status_, ref double bound)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CDFNBN evaluates the CDF of the Negative Binomial distribution
        //
        //  Discussion:
        //
        //    This routine calculates any one parameter of the negative binomial
        //    distribution given values for the others.
        //
        //    The cumulative negative binomial distribution returns the
        //    probability that there will be F or fewer failures before the
        //    S-th success in binomial trials each of which has probability of
        //    success PR.
        //
        //    The individual term of the negative binomial is the probability of
        //    F failures before S successes and is
        //    Choose( F, S+F-1 ) * PR^(S) * (1-PR)^F
        //
        //    Computation of other parameters involve a seach for a value that
        //    produces the desired value of P.  The search relies on the
        //    monotonicity of P with respect to the other parameters.
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
        //    1966, Formula 26.5.26.
        //
        //  Parameters:
        //
        //    Input, int WHICH, indicates which argument is to be calculated
        //    from the others.
        //    1: Calculate P and Q from F, S, PR and OMPR;
        //    2: Calculate F from P, Q, S, PR and OMPR;
        //    3: Calculate S from P, Q, F, PR and OMPR;
        //    4: Calculate PR and OMPR from P, Q, F and S.
        //
        //    Input/output, double P, the cumulation from 0 to F of
        //    the negative binomial distribution.  If P is an input value, it
        //    should lie in the range [0,1].
        //
        //    Input/output, double Q, equal to 1-P.  If Q is an input
        //    value, it should lie in the range [0,1].  If Q is an output value,
        //    it will lie in the range [0,1].
        //
        //    Input/output, double F, the upper limit of cumulation of
        //    the binomial distribution.  There are F or fewer failures before
        //    the S-th success.  If this is an input value, it may lie in the
        //    range [0,+infinity), and if it is an output value, it will be searched
        //    for in the range [0,1.0D+300].
        //
        //    Input/output, double S, the number of successes.
        //    If this is an input value, it should lie in the range: [0, +infinity).
        //    If it is an output value, it will be searched for in the range:
        //    [0, 1.0D+300].
        //
        //    Input/output, double PR, the probability of success in each
        //    binomial trial.  Whether an input or output value, it should lie in the
        //    range [0,1].
        //
        //    Input/output, double OMPR, the value of (1-PR).  Whether an
        //    input or output value, it should lie in the range [0,1].
        //
        //    Output, int STATUS, reports the status of the computation.
        //     0, if the calculation completed correctly;
        //    -I, if the input parameter number I is out of range;
        //    +1, if the answer appears to be lower than lowest search bound;
        //    +2, if the answer appears to be higher than greatest search bound;
        //    +3, if P + Q /= 1;
        //    +4, if PR + OMPR /= 1.
        //
        //    Output, double BOUND, is only defined if STATUS is nonzero.
        //    If STATUS is negative, then this is the value exceeded by parameter I.
        //    if STATUS is 1 or 2, this is the search bound that was exceeded.
        //
    {
        const double tol = 1.0e-8;
        const double atol = 1.0e-50;
        const double inf = 1.0e300;
        const double one = 1.0e0;

        double ccum = 0;
        double cum = 0;
        const int K1 = 1;
        const double K2 = 0.0e0;
        const double K4 = 0.5e0;
        const double K5 = 5.0e0;
        const double K11 = 1.0e0;
        bool qporq = false;

        E0000_E0001_Data data = new()
        {
            status = 0
        };

        bound = 0.0;
        switch (which is < 1 or > 4)
        {
            //
            //  Check arguments
            //
            case false:
                goto S30;
        }

        switch (which < 1)
        {
            case false:
                goto S10;
        }

        bound = 1.0e0;
        goto S20;
        S10:
        bound = 4.0e0;
        S20:
        data.status = -1;
        return;
        S30:
        switch (which)
        {
            case 1:
                goto S70;
        }

        switch (p is < 0.0e0 or > 1.0e0)
        {
            //
            //     P
            //
            case false:
                goto S60;
        }

        switch (p < 0.0e0)
        {
            case false:
                goto S40;
        }

        bound = 0.0e0;
        goto S50;
        S40:
        bound = 1.0e0;
        S50:
        data.status = -2;
        return;
        S70:
        S60:
        switch (which)
        {
            case 1:
                goto S110;
        }

        switch (q is <= 0.0e0 or > 1.0e0)
        {
            //
            //     Q
            //
            case false:
                goto S100;
        }

        switch (q <= 0.0e0)
        {
            case false:
                goto S80;
        }

        bound = 0.0e0;
        goto S90;
        S80:
        bound = 1.0e0;
        S90:
        data.status = -3;
        return;
        S110:
        S100:
        switch (which)
        {
            case 2:
                goto S130;
        }

        switch (s < 0.0e0)
        {
            //
            //     S
            //
            case false:
                goto S120;
        }

        bound = 0.0e0;
        data.status = -4;
        return;
        S130:
        S120:
        switch (which)
        {
            case 3:
                goto S150;
        }

        switch (xn < 0.0e0)
        {
            //
            //     XN
            //
            case false:
                goto S140;
        }

        bound = 0.0e0;
        data.status = -5;
        return;
        S150:
        S140:
        switch (which)
        {
            case 4:
                goto S190;
        }

        switch (pr is < 0.0e0 or > 1.0e0)
        {
            //
            //     PR
            //
            case false:
                goto S180;
        }

        switch (pr < 0.0e0)
        {
            case false:
                goto S160;
        }

        bound = 0.0e0;
        goto S170;
        S160:
        bound = 1.0e0;
        S170:
        data.status = -6;
        return;
        S190:
        S180:
        switch (which)
        {
            case 4:
                goto S230;
        }

        switch (ompr is < 0.0e0 or > 1.0e0)
        {
            //
            //     OMPR
            //
            case false:
                goto S220;
        }

        switch (ompr < 0.0e0)
        {
            case false:
                goto S200;
        }

        bound = 0.0e0;
        goto S210;
        S200:
        bound = 1.0e0;
        S210:
        data.status = -7;
        return;
        S230:
        S220:
        switch (which)
        {
            case 1:
                goto S270;
        }

        //
        //     P + Q
        //
        double pq = p + q;
        if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S260;
        }

        switch (pq < 0.0e0)
        {
            case false:
                goto S240;
        }

        bound = 0.0e0;
        goto S250;
        S240:
        bound = 1.0e0;
        S250:
        data.status = 3;
        return;
        S270:
        S260:
        switch (which)
        {
            case 4:
                goto S310;
        }

        //
        //     PR + OMPR
        //
        double prompr = pr + ompr;
        if (!(Math.Abs(prompr - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S300;
        }

        switch (prompr < 0.0e0)
        {
            case false:
                goto S280;
        }

        bound = 0.0e0;
        goto S290;
        S280:
        bound = 1.0e0;
        S290:
        data.status = 4;
        return;
        S310:
        S300:
        qporq = (which == 1) switch
        {
            false => p <= q,
            _ => qporq
        };

        switch (which)
        {
            //
            //     Select the minimum of P or Q
            //     Calculate ANSWERS
            //
            case 1:
                //
                //  Calculating P
                //
                cumnbn(s, xn, pr, ompr, ref p, ref q);
                data.status = 0;
                break;
            case 2:
            {
                //
                //     Calculating S
                //
                s = 5.0e0;
                E0000E0001.dstinv(ref data, K2, inf, K4, K4, K5, atol, tol);
                data.status = 0;
                data.x = s;
                E0000E0001.dinvr(ref data);
                S320:
                switch (data.status == 1)
                {
                    case false:
                        goto S350;
                }

                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S330;
                }

                data.fx = cum - p;
                goto S340;
                S330:
                data.fx = ccum - q;
                S340:
                data.x = s;
                E0000E0001.dinvr(ref data);
                goto S320;
                S350:
                switch (data.status == -1)
                {
                    case false:
                        goto S380;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S360;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S370;
                S360:
                data.status = 2;
                bound = inf;
                S380:
                S370: ;
                break;
            }
            case 3:
            {
                //
                //     Calculating XN
                //
                xn = 5.0e0;
                E0000E0001.dstinv(ref data, K2, inf, K4, K4, K5, atol, tol);
                data.status = 0;
                data.x = xn;
                E0000E0001.dinvr(ref data);
                S390:
                switch (data.status == 1)
                {
                    case false:
                        goto S420;
                }

                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S400;
                }

                data.fx = cum - p;
                goto S410;
                S400:
                data.fx = ccum - q;
                S410:
                data.x = xn;
                E0000E0001.dinvr(ref data);
                goto S390;
                S420:
                switch (data.status == -1)
                {
                    case false:
                        goto S450;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S430;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S440;
                S430:
                data.status = 2;
                bound = inf;
                S450:
                S440: ;
                break;
            }
            case 4:
            {
                //
                //     Calculating PR and OMPR
                //
                double T12 = atol;
                double T13 = tol;
                E0000E0001.dstzr(ref data, K2, K11, T12, T13);
                switch (qporq)
                {
                    case false:
                        goto S480;
                }

                data.status = 0;
                data.x = pr;
                E0000E0001.dzror(ref data);
                ompr = one - data.x;
                S460:
                switch (data.status == 1)
                {
                    case false:
                        goto S470;
                }

                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = pr;
                E0000E0001.dzror(ref data);
                ompr = one - data.x;
                goto S460;
                S470:
                goto S510;
                S480:
                data.status = 0;
                data.x = ompr;
                E0000E0001.dzror(ref data);
                pr = one - data.x;
                S490:
                switch (data.status == 1)
                {
                    case false:
                        goto S500;
                }

                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                data.fx = ccum - q;
                data.x = ompr;
                E0000E0001.dzror(ref data);
                pr = one - data.x;
                goto S490;
                S510:
                S500:
                switch (data.status == -1)
                {
                    case false:
                        goto S540;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S520;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S530;
                S520:
                data.status = 2;
                bound = 1.0e0;
                S530: ;
                break;
            }
        }

        S540:
        status_ = data.status;
    }
}