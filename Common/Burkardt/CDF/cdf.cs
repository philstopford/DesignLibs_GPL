using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfgam(int which, ref double p, ref double q, ref double x, ref double shape,
                ref double scale, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFGAM evaluates the CDF of the Gamma Distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the Gamma distribution
            //    given the others.
            //
            //    The cumulative distribution function P is calculated directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The gamma density is proportional to T^(SHAPE - 1) * EXP(- SCALE * T)
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
            //    Armido DiDinato and Alfred Morris,
            //    Computation of the incomplete gamma function ratios and their inverse,
            //    ACM Transactions on Mathematical Software,
            //    Volume 12, 1986, pages 377-393.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, SHAPE and SCALE;
            //    2: Calculate X from P, Q, SHAPE and SCALE;
            //    3: Calculate SHAPE from P, Q, X and SCALE;
            //    4: Calculate SCALE from P, Q, X and SHAPE.
            //
            //    Input/output, double *P, the integral from 0 to X of the
            //    Gamma density.  If this is an input value, it should lie in the
            //    range: [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration of
            //    the Gamma density.  If this is an input value, it should lie in the
            //    range: [0, +infinity).  If it is an output value, it will lie in
            //    the range: [0,1E300].
            //
            //    Input/output, double *SHAPE, the shape parameter of the
            //    Gamma density.  If this is an input value, it should lie in the range:
            //    (0, +infinity).  If it is an output value, it will be searched for
            //    in the range: [1.0D-300,1.0D+300].
            //
            //    Input/output, double *SCALE, the scale parameter of the
            //    Gamma density.  If this is an input value, it should lie in the range
            //    (0, +infinity).  If it is an output value, it will be searched for
            //    in the range: (1.0D-300,1.0D+300].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +10, if the Gamma or inverse Gamma routine cannot compute the answer.
            //    This usually happens only for X and SHAPE very large (more than 1.0D+10.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double inf = 1.0e300;

            double ccum = 0;
            double cum = 0;
            double fx = 0;
            int ierr = 0;
            int K1 = 1;
            double K5 = 0.5e0;
            double K6 = 5.0e0;
            double porq = 0;
            double pq;
            bool qhi = false;
            bool qleft = false;
            bool qporq = false;
            double T2;
            double T3;
            double T4;
            double T7;
            double T8;
            double T9;
            double xscale;
            double xx = 0;

            status = 0;
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
            status = -1;
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
            status = -2;
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
            status = -3;
            return;
            S110:
            S100:
            if (which == 2) goto S130;
            //
            //     X
            //
            if (!(x < 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //  SHAPE
            //
            if (!(shape <= 0.0e0)) goto S140;
            bound = 0.0e0;
            status = -5;
            return;
            S150:
            S140:
            if (which == 4) goto S170;
            //
            //  SCALE
            //
            if (!(scale <= 0.0e0)) goto S160;
            bound = 0.0e0;
            status = -6;
            return;
            S170:
            S160:
            if (which == 1) goto S210;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S200;
            if (!(pq < 0.0e0)) goto S180;
            bound = 0.0e0;
            goto S190;
            S180:
            bound = 1.0e0;
            S190:
            status = 3;
            return;
            S210:
            S200:
            if (which == 1) goto S240;
            //
            //     Select the minimum of P or Q
            //
            qporq = p <= q;
            if (!qporq) goto S220;
            porq = p;
            goto S230;
            S220:
            porq = q;
            S240:
            S230:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //     Calculating P
                //
                status = 0;
                xscale = x * scale;
                cumgam(xscale, shape, ref p, ref q);
                if (porq > 1.5e0) status = 10;
            }
            else if (2 == which)
            {
                //
                //     Computing X
                //
                T2 = -1.0e0;
                gamma_inc_inv(shape, ref xx, T2, p, q, ref ierr);
                if (ierr < 0.0e0)
                {
                    status = 10;
                    return;
                }
                else
                {
                    x = xx / scale;
                    status = 0;
                }
            }
            else if (3 == which)
            {
                //
                //     Computing SHAPE
                //
                shape = 5.0e0;
                xscale = x * scale;
                T3 = zero;
                T4 = inf;
                T7 = atol;
                T8 = tol;
                dstinv(T3, T4, K5, K5, K6, T7, T8);
                status = 0;
                dinvr(status, ref shape, fx, ref qleft, ref qhi);
                S250:
                if (!(status == 1)) goto S290;
                cumgam(xscale, shape, ref cum, ref ccum);
                if (!qporq) goto S260;
                fx = cum - p;
                goto S270;
                S260:
                fx = ccum - q;
                S270:
                if (!((qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))) goto S280;
                status = 10;
                return;
                S280:
                dinvr(status, ref shape, fx, ref qleft, ref qhi);
                goto S250;
                S290:
                if (!(status == -1)) goto S320;
                if (!qleft) goto S300;
                status = 1;
                bound = zero;
                goto S310;
                S300:
                status = 2;
                bound = inf;
                S320:
                S310: ;
            }
            else if (4 == which)
            {
                //
                //  Computing SCALE
                //
                T9 = -1.0e0;
                gamma_inc_inv(shape, ref xx, T9, p, q, ref ierr);
                if (ierr < 0.0e0)
                {
                    status = 10;
                    return;
                }
                else
                {
                    scale = xx / x;
                    status = 0;
                }
            }
        }

        public static void cdfnbn(int which, ref double p, ref double q, ref double s, ref double xn,
                ref double pr, ref double ompr, ref int status, ref double bound)

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
            double fx = 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double K11 = 1.0e0;
            double pq;
            double prompr = 0;
            bool qhi = false;
            bool qleft = false;
            bool qporq = false;
            double T3;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T12;
            double T13;
            double xhi = 0;
            double xlo = 0;

            status = 0;
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
            status = -1;
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
            status = -2;
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
            status = -3;
            return;
            S110:
            S100:
            if (which == 2) goto S130;
            //
            //     S
            //
            if (!(s < 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     XN
            //
            if (!(xn < 0.0e0)) goto S140;
            bound = 0.0e0;
            status = -5;
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
            status = -6;
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
            status = -7;
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
            status = 3;
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
            status = 4;
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
                cumnbn(s, xn, pr, ompr, ref p, ref q);
                status = 0;
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
                dstinv(K2, T3, K4, K4, K5, T6, T7);
                status = 0;
                dinvr(status, ref s, fx, ref qleft, ref qhi);
                S320:
                if (!(status == 1)) goto S350;
                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S330;
                fx = cum - p;
                goto S340;
                S330:
                fx = ccum - q;
                S340:
                dinvr(status, ref s, fx, ref qleft, ref qhi);
                goto S320;
                S350:
                if (!(status == -1)) goto S380;
                if (!qleft) goto S360;
                status = 1;
                bound = 0.0e0;
                goto S370;
                S360:
                status = 2;
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
                dstinv(K2, T8, K4, K4, K5, T9, T10);
                status = 0;
                dinvr(status, ref xn, fx, ref qleft, ref qhi);
                S390:
                if (!(status == 1)) goto S420;
                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S400;
                fx = cum - p;
                goto S410;
                S400:
                fx = ccum - q;
                S410:
                dinvr(status, ref xn, fx, ref qleft, ref qhi);
                goto S390;
                S420:
                if (!(status == -1)) goto S450;
                if (!qleft) goto S430;
                status = 1;
                bound = 0.0e0;
                goto S440;
                S430:
                status = 2;
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
                dstzr(K2, K11, T12, T13);
                if (!qporq) goto S480;
                status = 0;
                dzror(ref status, ref pr, ref fx, ref xlo, ref xhi, ref qleft, ref qhi);
                ompr = one - pr;
                S460:
                if (!(status == 1)) goto S470;
                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                fx = cum - p;
                dzror(ref status, ref pr, ref fx, ref xlo, ref xhi, ref qleft, ref qhi);
                ompr = one - pr;
                goto S460;
                S470:
                goto S510;
                S480:
                status = 0;
                dzror(ref status, ref ompr, ref fx, ref xlo, ref xhi, ref qleft, ref qhi);
                pr = one - ompr;
                S490:
                if (!(status == 1)) goto S500;
                cumnbn(s, xn, pr, ompr, ref cum, ref ccum);
                fx = ccum - q;
                dzror(ref status, ref ompr, ref fx, ref xlo, ref xhi, ref qleft, ref qhi);
                pr = one - ompr;
                goto S490;
                S510:
                S500:
                if (!(status == -1)) goto S540;
                if (!qleft) goto S520;
                status = 1;
                bound = 0.0e0;
                goto S530;
                S520:
                status = 2;
                bound = 1.0e0;
                S530: ;
            }

            S540:
            return;
        }

        public static void cdfnor(int which, ref double p, ref double q, ref double x, ref double mean,
                ref double sd, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFNOR evaluates the CDF of the Normal distribution.
            //
            //  Discussion:
            //
            //    A slightly modified version of ANORM from SPECFUN
            //    is used to calculate the cumulative standard normal distribution.
            //
            //    The rational functions from pages 90-95 of Kennedy and Gentle
            //    are used as starting values to Newton's Iterations which
            //    compute the inverse standard normal.  Therefore no searches are
            //    necessary for any parameter.
            //
            //    For X < -15, the asymptotic expansion for the normal is used  as
            //    the starting value in finding the inverse standard normal.
            //
            //    The normal density is proportional to
            //    exp( - 0.5 * (( X - MEAN)/SD)^2)
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
            //    1966, Formula 26.2.12.
            //
            //    William Cody,
            //    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
            //      Special Function Routines and Test Drivers,
            //    ACM Transactions on Mathematical Software,
            //    Volume 19, pages 22-32, 1993.
            //
            //    Kennedy and Gentle,
            //    Statistical Computing,
            //    Marcel Dekker, NY, 1980,
            //    QA276.4  K46
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, MEAN and SD;
            //    2: Calculate X from P, Q, MEAN and SD;
            //    3: Calculate MEAN from P, Q, X and SD;
            //    4: Calculate SD from P, Q, X and MEAN.
            //
            //    Input/output, double *P, the integral from -infinity to X
            //    of the Normal density.  If this is an input or output value, it will
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration of
            //    the Normal density.
            //
            //    Input/output, double *MEAN, the mean of the Normal density.
            //
            //    Input/output, double *SD, the standard deviation of the
            //    Normal density.  If this is an input value, it should lie in the
            //    range (0,+infinity).
            //
            //    Output, int *STATUS, the status of the calculation.
            //    0, if calculation completed correctly;
            //    -I, if input parameter number I is out of range;
            //    1, if answer appears to be lower than lowest search bound;
            //    2, if answer appears to be higher than greatest search bound;
            //    3, if P + Q /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            int K1 = 1;
            double pq;
            double z;

            status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            status = 0;
            if (!(which < 1 || which > 4)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 4.0e0;
            S20:
            status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //     P
            //
            if (!(p <= 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p <= 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = 1.0e0;
            S50:
            status = -2;
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
            status = -3;
            return;
            S110:
            S100:
            if (which == 1) goto S150;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S140;
            if (!(pq < 0.0e0)) goto S120;
            bound = 0.0e0;
            goto S130;
            S120:
            bound = 1.0e0;
            S130:
            status = 3;
            return;
            S150:
            S140:
            if (which == 4) goto S170;
            //
            //     SD
            //
            if (!(sd <= 0.0e0)) goto S160;
            bound = 0.0e0;
            status = -6;
            return;
            S170:
            S160:
            //
            //  Computing P
            //
            if (1 == which)
            {
                z = (x - mean) / sd;
                cumnor(z, ref p, ref q);
            }
            //
            //  Computing X
            //
            else if (2 == which)
            {
                z = dinvnr(p, q);
                x = sd * z + mean;
            }
            //
            //  Computing the MEAN
            //
            else if (3 == which)
            {
                z = dinvnr(p, q);
                mean = x - sd * z;
            }
            //
            //  Computing SD
            //
            else if (4 == which)
            {
                z = dinvnr(p, q);
                sd = (x - mean) / z;
            }

            return;
        }

        public static void cdfpoi(int which, ref double p, ref double q, ref double s, ref double xlam,
                ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFPOI evaluates the CDF of the Poisson distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the Poisson distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
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
            //    1966, Formula 26.4.21.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from S and XLAM;
            //    2: Calculate A from P, Q and XLAM;
            //    3: Calculate XLAM from P, Q and S.
            //
            //    Input/output, double *P, the cumulation from 0 to S of the
            //    Poisson density.  Whether this is an input or output value, it will
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *S, the upper limit of cumulation of
            //    the Poisson CDF.  If this is an input value, it should lie in
            //    the range: [0, +infinity).  If it is an output value, it will be
            //    searched for in the range: [0,1.0D+300].
            //
            //    Input/output, double *XLAM, the mean of the Poisson
            //    distribution.  If this is an input value, it should lie in the range
            //    [0, +infinity).  If it is an output value, it will be searched for
            //    in the range: [0,1E300].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double inf = 1.0e300;

            double ccum = 0;
            double cum = 0;
            double fx = 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double pq;
            bool qhi = false;
            bool qleft = false;
            bool qporq = false;
            double T3;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;

            status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 3)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 3.0e0;
            S20:
            status = -1;
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
            status = -2;
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
            status = -3;
            return;
            S110:
            S100:
            if (which == 2) goto S130;
            //
            //     S
            //
            if (!(s < 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     XLAM
            //
            if (!(xlam < 0.0e0)) goto S140;
            bound = 0.0e0;
            status = -5;
            return;
            S150:
            S140:
            if (which == 1) goto S190;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S180;
            if (!(pq < 0.0e0)) goto S160;
            bound = 0.0e0;
            goto S170;
            S160:
            bound = 1.0e0;
            S170:
            status = 3;
            return;
            S190:
            S180:
            if (!(which == 1)) qporq = p <= q;
            //
            //  Select the minimum of P or Q
            //  Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P
                //
                cumpoi(s, xlam, ref p, ref q);
                status = 0;
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
                dstinv(K2, T3, K4, K4, K5, T6, T7);
                status = 0;
                dinvr(status, ref s, fx, ref qleft, ref qhi);
                S200:
                if (!(status == 1)) goto S230;
                cumpoi(s, xlam, ref cum, ref ccum);
                if (!qporq) goto S210;
                fx = cum - p;
                goto S220;
                S210:
                fx = ccum - q;
                S220:
                dinvr(status, ref s, fx, ref qleft, ref qhi);
                goto S200;
                S230:
                if (!(status == -1)) goto S260;
                if (!qleft) goto S240;
                status = 1;
                bound = 0.0e0;
                goto S250;
                S240:
                status = 2;
                bound = inf;
                S260:
                S250: ;
            }
            else if (3 == which)
            {
                //
                //     Calculating XLAM
                //
                xlam = 5.0e0;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                dstinv(K2, T8, K4, K4, K5, T9, T10);
                status = 0;
                dinvr(status, ref xlam, fx, ref qleft, ref qhi);
                S270:
                if (!(status == 1)) goto S300;
                cumpoi(s, xlam, ref cum, ref ccum);
                if (!qporq) goto S280;
                fx = cum - p;
                goto S290;
                S280:
                fx = ccum - q;
                S290:
                dinvr(status, ref xlam, fx, ref qleft, ref qhi);
                goto S270;
                S300:
                if (!(status == -1)) goto S330;
                if (!qleft) goto S310;
                status = 1;
                bound = 0.0e0;
                goto S320;
                S310:
                status = 2;
                bound = inf;
                S320: ;
            }

            S330:
            return;
        }

        public static void cdft(int which, ref double p, ref double q, ref double t, ref double df,
                ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFT evaluates the CDF of the T distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the T distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of other parameters involve a seach for a value that
            //    produces the desired value of P.   The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The original version of this routine allowed the search interval
            //    to extend from -1.0E+300 to +1.0E+300, which is fine until you
            //    try to evaluate a function at such a point!
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
            //    1966, Formula 26.5.27.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1 : Calculate P and Q from T and DF;
            //    2 : Calculate T from P, Q and DF;
            //    3 : Calculate DF from P, Q and T.
            //
            //    Input/output, double *P, the integral from -infinity to T of
            //    the T-density.  Whether an input or output value, this will lie in the
            //    range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *T, the upper limit of integration of
            //    the T-density.  If this is an input value, it may have any value.
            //    It it is an output value, it will be searched for in the range
            //    [ -1.0D+30, 1.0D+30 ].
            //
            //    Input/output, double *DF, the number of degrees of freedom
            //    of the T distribution.  If this is an input value, it should lie
            //    in the range: (0 , +infinity).  If it is an output value, it will be
            //    searched for in the range: [1, 1.0D+10].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double inf = 1.0e30;
            double maxdf = 1.0e10;

            double ccum = 0;
            double cum = 0;
            double fx = 0;
            int K1 = 1;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double pq;
            bool qhi = false;
            bool qleft = false;
            bool qporq = false;
            double T2;
            double T3;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11;

            status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 3)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 3.0e0;
            S20:
            status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //     P
            //
            if (!(p <= 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p <= 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = 1.0e0;
            S50:
            status = -2;
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
            status = -3;
            return;
            S110:
            S100:
            if (which == 3) goto S130;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -5;
            return;
            S130:
            S120:
            if (which == 1) goto S170;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S160;
            if (!(pq < 0.0e0)) goto S140;
            bound = 0.0e0;
            goto S150;
            S140:
            bound = 1.0e0;
            S150:
            status = 3;
            return;
            S170:
            S160:
            if (!(which == 1)) qporq = p <= q;
            //
            //  Select the minimum of P or Q.  Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Computing P and Q
                //
                cumt(t, df, ref p, ref q);
                status = 0;
            }
            else if (2 == which)
            {
                //
                //  Computing T
                //  Get initial approximation for T
                //
                t = dt1(p, q, df);
                T2 = -inf;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                dstinv(T2, T3, K4, K4, K5, T6, T7);
                status = 0;
                dinvr(status, ref t, fx, ref qleft, ref qhi);
                S180:
                if (!(status == 1)) goto S210;
                cumt(t, df, ref cum, ref ccum);
                if (!qporq) goto S190;
                fx = cum - p;
                goto S200;
                S190:
                fx = ccum - q;
                S200:
                dinvr(status, ref t, fx, ref qleft, ref qhi);
                goto S180;
                S210:
                if (!(status == -1)) goto S240;
                if (!qleft) goto S220;
                status = 1;
                bound = -inf;
                goto S230;
                S220:
                status = 2;
                bound = inf;
                S240:
                S230: ;
            }
            else if (3 == which)
            {
                //
                //  Computing DF
                //
                df = 5.0e0;
                T8 = zero;
                T9 = maxdf;
                T10 = atol;
                T11 = tol;
                dstinv(T8, T9, K4, K4, K5, T10, T11);
                status = 0;
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                S250:
                if (!(status == 1)) goto S280;
                cumt(t, df, ref cum, ref ccum);
                if (!qporq) goto S260;
                fx = cum - p;
                goto S270;
                S260:
                fx = ccum - q;
                S270:
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                goto S250;
                S280:
                if (!(status == -1)) goto S310;
                if (!qleft) goto S290;
                status = 1;
                bound = zero;
                goto S300;
                S290:
                status = 2;
                bound = maxdf;
                S300: ;
            }

            S310:
            return;
        }


    }
}