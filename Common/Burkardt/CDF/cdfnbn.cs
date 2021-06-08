using System;

namespace Burkardt.CDFLib
{
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
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double inf = 1.0e300;
            double one = 1.0e0;

            double ccum = 0;
            double cum = 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double K11 = 1.0e0;
            double pq = 0;
            double prompr = 0;
            bool qporq = false;
            double T3 = 0;
            double T6 = 0;
            double T7 = 0;
            double T8 = 0;
            double T9 = 0;
            double T10 = 0;
            double T12 = 0;
            double T13 = 0;

            E0000_E0001_Data data = new E0000_E0001_Data();

            data.status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 4)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 4.0e0;
            S20:
            data.status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //     P
            //
            if (!(p < 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = 1.0e0;
            S50:
            data.status = -2;
            return;
            S70:
            S60:
            if (which == 1) goto S110;
            //
            //     Q
            //
            if (!(q <= 0.0e0 || q > 1.0e0)) goto S100;
            if (!(q <= 0.0e0)) goto S80;
            bound = 0.0e0;
            goto S90;
            S80:
            bound = 1.0e0;
            S90:
            data.status = -3;
            return;
            S110:
            S100:
            if (which == 2) goto S130;
            //
            //     S
            //
            if (!(s < 0.0e0)) goto S120;
            bound = 0.0e0;
            data.status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     XN
            //
            if (!(xn < 0.0e0)) goto S140;
            bound = 0.0e0;
            data.status = -5;
            return;
            S150:
            S140:
            if (which == 4) goto S190;
            //
            //     PR
            //
            if (!(pr < 0.0e0 || pr > 1.0e0)) goto S180;
            if (!(pr < 0.0e0)) goto S160;
            bound = 0.0e0;
            goto S170;
            S160:
            bound = 1.0e0;
            S170:
            data.status = -6;
            return;
            S190:
            S180:
            if (which == 4) goto S230;
            //
            //     OMPR
            //
            if (!(ompr < 0.0e0 || ompr > 1.0e0)) goto S220;
            if (!(ompr < 0.0e0)) goto S200;
            bound = 0.0e0;
            goto S210;
            S200:
            bound = 1.0e0;
            S210:
            data.status = -7;
            return;
            S230:
            S220:
            if (which == 1) goto S270;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S260;
            if (!(pq < 0.0e0)) goto S240;
            bound = 0.0e0;
            goto S250;
            S240:
            bound = 1.0e0;
            S250:
            data.status = 3;
            return;
            S270:
            S260:
            if (which == 4) goto S310;
            //
            //     PR + OMPR
            //
            prompr = pr + ompr;
            if (!(Math.Abs(prompr - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S300;
            if (!(prompr < 0.0e0)) goto S280;
            bound = 0.0e0;
            goto S290;
            S280:
            bound = 1.0e0;
            S290:
            data.status = 4;
            return;
            S310:
            S300:
            if (!(which == 1)) qporq = p <= q;
            //
            //     Select the minimum of P or Q
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P
                //
                cumnbn(s, xn, pr, ompr, p, q);
                data.status = 0;
            }
            else if (2 == which)
            {
                //
                //     Calculating S
                //
                s = 5.0e0;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                E0000E0001.dstinv(ref data, K2, T3, K4, K4, K5, T6, T7);
                data.status = 0;
                data.x = s;
                E0000E0001.dinvr(ref data);
                S320:
                if (!(data.status == 1)) goto S350;
                cumnbn(s, xn, pr, ompr, cum, ccum);
                if (!qporq) goto S330;
                data.fx = cum - p;
                goto S340;
                S330:
                data.fx = ccum - q;
                S340:
                data.x = s;
                E0000E0001.dinvr(ref data);
                goto S320;
                S350:
                if (!(data.status == -1)) goto S380;
                if (!data.qleft) goto S360;
                data.status = 1;
                bound = 0.0e0;
                goto S370;
                S360:
                data.status = 2;
                bound = inf;
                S380:
                S370: ;
            }
            else if (3 == which)
            {
                //
                //     Calculating XN
                //
                xn = 5.0e0;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                E0000E0001.dstinv(ref data, K2, T8, K4, K4, K5, T9, T10);
                data.status = 0;
                data.x = xn;
                E0000E0001.dinvr(ref data);
                S390:
                if (!(data.status == 1)) goto S420;
                cumnbn(s, xn, pr, ompr, cum, ccum);
                if (!qporq) goto S400;
                data.fx = cum - p;
                goto S410;
                S400:
                data.fx = ccum - q;
                S410:
                data.x = xn;
                E0000E0001.dinvr(ref data);
                goto S390;
                S420:
                if (!(data.status == -1)) goto S450;
                if (!data.qleft) goto S430;
                data.status = 1;
                bound = 0.0e0;
                goto S440;
                S430:
                data.status = 2;
                bound = inf;
                S450:
                S440: ;
            }
            else if (4 == which)
            {
                //
                //     Calculating PR and OMPR
                //
                T12 = atol;
                T13 = tol;
                E0000E0001.dstzr(ref data, K2, K11, T12, T13);
                if (!qporq) goto S480;
                data.status = 0;
                data.x = pr;
                E0000E0001.dzror(ref data);
                ompr = one - data.x;
                S460:
                if (!(data.status == 1)) goto S470;
                cumnbn(s, xn, pr, ompr, cum, ccum);
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
                if (!(data.status == 1)) goto S500;
                cumnbn(s, xn, pr, ompr, cum, ccum);
                data.fx = ccum - q;
                data.x = ompr;
                E0000E0001.dzror(ref data);
                pr = one - data.x;
                goto S490;
                S510:
                S500:
                if (!(data.status == -1)) goto S540;
                if (!data.qleft) goto S520;
                data.status = 1;
                bound = 0.0e0;
                goto S530;
                S520:
                data.status = 2;
                bound = 1.0e0;
                S530: ;
            }

            S540:
            status_ = data.status;
            return;
        }
    }
}