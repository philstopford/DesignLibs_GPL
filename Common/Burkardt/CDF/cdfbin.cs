using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfbin(int which, ref double p, ref double q, ref double s, ref double xn,
                ref double pr, ref double ompr, ref int status_, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFBIN evaluates the CDF of the Binomial distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the binomial distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    P is the probablility of S or fewer successes in XN binomial trials,
            //    each trial having an individual probability of success of PR.
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
            //    Input, int *WHICH, indicates which of argument values is to
            //    be calculated from the others.
            //    1: Calculate P and Q from S, XN, PR and OMPR;
            //    2: Calculate S from P, Q, XN, PR and OMPR;
            //    3: Calculate XN from P, Q, S, PR and OMPR;
            //    4: Calculate PR and OMPR from P, Q, S and XN.
            //
            //    Input/output, double *P, the cumulation, from 0 to S,
            //    of the binomial distribution.  If P is an input value, it should
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *S, the number of successes observed.
            //    Whether this is an input or output value, it should lie in the
            //    range [0,XN].
            //
            //    Input/output, double *XN, the number of binomial trials.
            //    If this is an input value it should lie in the range: (0, +infinity).
            //    If it is an output value it will be searched for in the
            //    range [1.0D-300, 1.0D+300].
            //
            //    Input/output, double *PR, the probability of success in each
            //    binomial trial.  Whether this is an input or output value, it should
            //    lie in the range: [0,1].
            //
            //    Input/output, double *OMPR, equal to 1-PR.  Whether this is an
            //    input or output value, it should lie in the range [0,1].  Also, it should
            //    be the case that PR + OMPR = 1.
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +4, if PR + OMPR /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double atol = (1.0e-50);
            double tol = (1.0e-8);
            double zero = (1.0e-300);
            double inf = 1.0e300;
            double one = 1.0e0;

            double ccum= 0;
            double cum= 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K3 = 0.5e0;
            double K4 = 5.0e0;
            double K11 = 1.0e0;
            double pq= 0;
            double prompr= 0;
            bool qporq = false;
            double T5= 0;
            double T6= 0;
            double T7= 0;
            double T8= 0;
            double T9= 0;
            double T10= 0;
            double T12= 0;
            double T13= 0;

            E0000_E0001_Data data = new E0000_E0001_Data();

            bound = 0.0;
            //
            //  Check arguments
            //
            if (!(which < 1 && which > 4)) goto S30;
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
            if (!(q < 0.0e0 || q > 1.0e0)) goto S100;
            if (!(q < 0.0e0)) goto S80;
            bound = 0.0e0;
            goto S90;
            S80:
            bound = 1.0e0;
            S90:
            data.status = -3;
            return;
            S110:
            S100:
            if (which == 3) goto S130;
            //
            //     XN
            //
            if (!(xn <= 0.0e0)) goto S120;
            bound = 0.0e0;
            data.status = -5;
            return;
            S130:
            S120:
            if (which == 2) goto S170;
            //
            //     S
            //
            if (!(s < 0.0e0 || (which != 3 && s > xn))) goto S160;
            if (!(s < 0.0e0)) goto S140;
            bound = 0.0e0;
            goto S150;
            S140:
            bound = xn;
            S150:
            data.status = -4;
            return;
            S170:
            S160:
            if (which == 4) goto S210;
            //
            //     PR
            //
            if (!(pr < 0.0e0 || pr > 1.0e0)) goto S200;
            if (!(pr < 0.0e0)) goto S180;
            bound = 0.0e0;
            goto S190;
            S180:
            bound = 1.0e0;
            S190:
            data.status = -6;
            return;
            S210:
            S200:
            if (which == 4) goto S250;
            //
            //     OMPR
            //
            if (!(ompr < 0.0e0 || ompr > 1.0e0)) goto S240;
            if (!(ompr < 0.0e0)) goto S220;
            bound = 0.0e0;
            goto S230;
            S220:
            bound = 1.0e0;
            S230:
            data.status = -7;
            return;
            S250:
            S240:
            if (which == 1) goto S290;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S280;
            if (!(pq < 0.0e0)) goto S260;
            bound = 0.0e0;
            goto S270;
            S260:
            bound = 1.0e0;
            S270:
            data.status = 3;
            return;
            S290:
            S280:
            if (which == 4) goto S330;
            //
            //     PR + OMPR
            //
            prompr = pr + ompr;
            if (!(Math.Abs(prompr - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S320;
            if (!(prompr < 0.0e0)) goto S300;
            bound = 0.0e0;
            goto S310;
            S300:
            bound = 1.0e0;
            S310:
            data.status = 4;
            return;
            S330:
            S320:
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
                cumbin(s, xn, pr, ompr, ref p, ref q);
                data.status = 0;
            }
            else if (2 == which)
            {
                //
                //  Calculating S
                //
                s = 5.0e0;
                T5 = atol;
                T6 = tol;
                E0000E0001.dstinv(ref data, K2, xn, K3, K3, K4, T5, T6);
                data.status = 0;
                data.x = s;
                E0000E0001.dinvr(ref data);
                S340:
                if (!(data.status == 1)) goto S370;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S350;
                data.fx = cum - p;
                goto S360;
                S350:
                data.fx = ccum - q;
                S360:
                data.x = s;
                E0000E0001.dinvr(ref data);
                goto S340;
                S370:
                if (!(data.status == -1)) goto S400;
                if (!data.qleft) goto S380;
                data.status = 1;
                bound = 0.0e0;
                goto S390;
                S380:
                data.status = 2;
                bound = xn;
                S400:
                S390: ;
            }
            else if (3 == which)
            {
                //
                //  Calculating XN
                //
                xn = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                E0000E0001.dstinv(ref data, T7, T8, K3, K3, K4, T9, T10);
                data.status = 0;
                data.x = xn;
                E0000E0001.dinvr(ref data);
                S410:
                if (!(data.status == 1)) goto S440;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S420;
                data.fx = cum - p;
                goto S430;
                S420:
                data.fx = ccum - q;
                S430:
                data.x = xn;
                E0000E0001.dinvr(ref data);
                goto S410;
                S440:
                if (!(data.status == -1)) goto S470;
                if (!data.qleft) goto S450;
                data.status = 1;
                bound = zero;
                goto S460;
                S450:
                data.status = 2;
                bound = inf;
                S470:
                S460: ;
            }
            else if (4 == which)
            {
                //
                //  Calculating PR and OMPR
                //
                T12 = atol;
                T13 = tol;
                E0000E0001.dstzr(ref data, K2, K11, T12, T13);
                if (!qporq) goto S500;
                data.status = 0;
                data.x = pr;
                E0000E0001.dzror(ref data);
                ompr = one - pr;
                S480:
                if (!(data.status == 1)) goto S490;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = pr;
                E0000E0001.dzror(ref data);
                ompr = one - pr;
                goto S480;
                S490:
                goto S530;
                S500:
                data.status = 0;
                data.x = ompr;
                E0000E0001.dzror(ref data);
                pr = one - ompr;
                S510:
                if (!(data.status == 1)) goto S520;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                data.fx = ccum - q;
                data.x = ompr;
                E0000E0001.dzror(ref data);
                pr = one - ompr;
                goto S510;
                S530:
                S520:
                if (!(data.status == -1)) goto S560;
                if (!data.qleft) goto S540;
                data.status = 1;
                bound = 0.0e0;
                goto S550;
                S540:
                data.status = 2;
                bound = 1.0e0;
                S550: ;
            }

            S560:
            status_ = data.status;
        }
    }
}