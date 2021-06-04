using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfchi(int which, ref double p, ref double q, ref double x, ref double df,
                ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFCHI evaluates the CDF of the chi square distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the chi square distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The CDF of the chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF [ ChiSquareDistribution [ DF ], X ]
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
            //    1966, Formula 26.4.19.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X and DF;
            //    2: Calculate X from P, Q and DF;
            //    3: Calculate DF from P, Q and X.
            //
            //    Input/output, double *P, the integral from 0 to X of
            //    the chi-square distribution.  If this is an input value, it should
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration
            //    of the chi-square distribution.  If this is an input
            //    value, it should lie in the range: [0, +infinity).  If it is an output
            //    value, it will be searched for in the range: [0,1.0D+300].
            //
            //    Input/output, double *DF, the degrees of freedom of the
            //    chi-square distribution.  If this is an input value, it should lie
            //    in the range: (0, +infinity).  If it is an output value, it will be
            //    searched for in the range: [ 1.0D-300, 1.0D+300].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +10, an error was returned from CUMGAM.
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
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double porq = 0;
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
            //  P
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
            //     DF
            //
            if (!(df <= 0.0e0)) goto S140;
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
            if (which == 1) goto S220;
            //
            //  Select the minimum of P or Q
            //
            qporq = p <= q;
            if (!qporq) goto S200;
            porq = p;
            goto S210;
            S200:
            porq = q;
            S220:
            S210:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P and Q
                //
                status = 0;
                cumchi(x, df, ref p, ref q);
                if (porq > 1.5e0)
                {
                    status = 10;
                    return;
                }
            }
            else if (2 == which)
            {
                //
                //  Calculating X
                //
                x = 5.0e0;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                dstinv(K2, T3, K4, K4, K5, T6, T7);
                status = 0;
                dinvr(status, ref x, fx, ref qleft, ref qhi);
                S230:
                if (!(status == 1)) goto S270;
                cumchi(x, df, ref cum, ref ccum);
                if (!qporq) goto S240;
                fx = cum - p;
                goto S250;
                S240:
                fx = ccum - q;
                S250:
                if (!(fx + porq > 1.5e0)) goto S260;
                status = 10;
                return;
                S260:
                dinvr(status, ref x, fx, ref qleft, ref qhi);
                goto S230;
                S270:
                if (!(status == -1)) goto S300;
                if (!qleft) goto S280;
                status = 1;
                bound = 0.0e0;
                goto S290;
                S280:
                status = 2;
                bound = inf;
                S300:
                S290: ;
            }
            else if (3 == which)
            {
                //
                //  Calculating DF
                //
                df = 5.0e0;
                T8 = zero;
                T9 = inf;
                T10 = atol;
                T11 = tol;
                dstinv(T8, T9, K4, K4, K5, T10, T11);
                status = 0;
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                S310:
                if (!(status == 1)) goto S350;
                cumchi(x, df, ref cum, ref ccum);
                if (!qporq) goto S320;
                fx = cum - p;
                goto S330;
                S320:
                fx = ccum - q;
                S330:
                if (!(fx + porq > 1.5e0)) goto S340;
                status = 10;
                return;
                S340:
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                goto S310;
                S350:
                if (!(status == -1)) goto S380;
                if (!qleft) goto S360;
                status = 1;
                bound = zero;
                goto S370;
                S360:
                status = 2;
                bound = inf;
                S370: ;
            }

            S380:
            return;
        }

        public static void cdfchn(int which, ref double p, ref double q, double x, ref double df,
                ref double pnonc, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the noncentral chi-square
            //    distribution given values for the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The computation time required for this routine is proportional
            //    to the noncentrality parameter (PNONC).  Very large values of
            //    this parameter can consume immense computer resources.  This is
            //    why the search range is bounded by 10,000.
            //
            //    The CDF of the noncentral chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
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
            //    1966, Formula 26.5.25.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, DF and PNONC;
            //    2: Calculate X from P, DF and PNONC;
            //    3: Calculate DF from P, X and PNONC;
            //    4: Calculate PNONC from P, X and DF.
            //
            //    Input/output, double *P, the integral from 0 to X of
            //    the noncentral chi-square distribution.  If this is an input
            //    value, it should lie in the range: [0, 1.0-1.0D-16).
            //
            //    Input/output, double *Q, is generally not used by this
            //    subroutine and is only included for similarity with other routines.
            //    However, if P is to be computed, then a value will also be computed
            //    for Q.
            //
            //    Input, double *X, the upper limit of integration of the
            //    noncentral chi-square distribution.  If this is an input value, it
            //    should lie in the range: [0, +infinity).  If it is an output value,
            //    it will be sought in the range: [0,1.0D+300].
            //
            //    Input/output, double *DF, the number of degrees of freedom
            //    of the noncentral chi-square distribution.  If this is an input value,
            //    it should lie in the range: (0, +infinity).  If it is an output value,
            //    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
            //
            //    Input/output, double *PNONC, the noncentrality parameter of
            //    the noncentral chi-square distribution.  If this is an input value, it
            //    should lie in the range: [0, +infinity).  If it is an output value,
            //    it will be searched for in the range: [0,1.0D+4]
            //
            //    Output, int *STATUS, reports on the calculation.
            //    0, if calculation completed correctly;
            //    -I, if input parameter number I is out of range;
            //    1, if the answer appears to be lower than the lowest search bound;
            //    2, if the answer appears to be higher than the greatest search bound.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tent4 = 1.0e4;
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double one = (1.0e0 - 1.0e-16);
            double inf = 1.0e300;

            double cum = 0;
            double ccum = 0;
            double fx = 0;
            double K1 = 0.0e0;
            double K3 = 0.5e0;
            double K4 = 5.0e0;
            bool qhi = false;
            bool qleft = false;
            double T2;
            double T5;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11;
            double T12;
            double T13;

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
            //  P
            //
            if (!(p < 0.0e0 || p > one)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = one;
            S50:
            status = -2;
            return;
            S70:
            S60:
            if (which == 2) goto S90;
            //
            //     X
            //
            if (!(x < 0.0e0)) goto S80;
            bound = 0.0e0;
            status = -4;
            return;
            S90:
            S80:
            if (which == 3) goto S110;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S100;
            bound = 0.0e0;
            status = -5;
            return;
            S110:
            S100:
            if (which == 4) goto S130;
            //
            //  PNONC
            //
            if (!(pnonc < 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -6;
            return;
            S130:
            S120:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //     Calculating P and Q
                //
                cumchn(x, df, pnonc, ref p, ref q);
                status = 0;
            }
            else if (2 == which)
            {
                //
                //     Calculating X
                //
                x = 5.0e0;
                T2 = inf;
                T5 = atol;
                T6 = tol;
                dstinv(K1, T2, K3, K3, K4, T5, T6);
                status = 0;
                dinvr(status, ref x, fx, ref qleft, ref qhi);
                S140:
                if (!(status == 1)) goto S150;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref x, fx, ref qleft, ref qhi);
                goto S140;
                S150:
                if (!(status == -1)) goto S180;
                if (!qleft) goto S160;
                status = 1;
                bound = 0.0e0;
                goto S170;
                S160:
                status = 2;
                bound = inf;
                S180:
                S170: ;
            }
            else if (3 == which)
            {
                //
                //     Calculating DF
                //
                df = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                dstinv(T7, T8, K3, K3, K4, T9, T10);
                status = 0;
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                S190:
                if (!(status == 1)) goto S200;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref df, fx, ref qleft, ref qhi);
                goto S190;
                S200:
                if (!(status == -1)) goto S230;
                if (!qleft) goto S210;
                status = 1;
                bound = zero;
                goto S220;
                S210:
                status = 2;
                bound = inf;
                S230:
                S220: ;
            }
            else if (4 == which)
            {
                //
                //     Calculating PNONC
                //
                pnonc = 5.0e0;
                T11 = tent4;
                T12 = atol;
                T13 = tol;
                dstinv(K1, T11, K3, K3, K4, T12, T13);
                status = 0;
                dinvr(status, ref pnonc, fx, ref qleft, ref qhi);
                S240:
                if (!(status == 1)) goto S250;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref pnonc, fx, ref qleft, ref qhi);
                goto S240;
                S250:
                if (!(status == -1)) goto S280;
                if (!qleft) goto S260;
                status = 1;
                bound = zero;
                goto S270;
                S260:
                status = 2;
                bound = tent4;
                S270: ;
            }

            S280:
            return;
        }
    }
}