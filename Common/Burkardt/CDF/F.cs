using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdff(int which, ref double p, ref double q, ref double f, ref double dfn,
                ref double dfd, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFF evaluates the CDF of the F distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the F distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The value of the cumulative F distribution is not necessarily
            //    monotone in either degree of freedom.  There thus may be two
            //    values that provide a given CDF value.  This routine assumes
            //    monotonicity and will find an arbitrary one of the two values.
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
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions
            //    1966, Formula 26.6.2.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from F, DFN and DFD;
            //    2: Calculate F from P, Q, DFN and DFD;
            //    3: Calculate DFN from P, Q, F and DFD;
            //    4: Calculate DFD from P, Q, F and DFN.
            //
            //    Input/output, double *P, the integral from 0 to F of
            //    the F-density.  If it is an input value, it should lie in the
            //    range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *F, the upper limit of integration
            //    of the F-density.  If this is an input value, it should lie in the
            //    range [0, +infinity).  If it is an output value, it will be searched
            //    for in the range [0,1.0D+300].
            //
            //    Input/output, double *DFN, the number of degrees of
            //    freedom of the numerator sum of squares.  If this is an input value,
            //    it should lie in the range: (0, +infinity).  If it is an output value,
            //    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
            //
            //    Input/output, double *DFD, the number of degrees of freedom
            //    of the denominator sum of squares.  If this is an input value, it should
            //    lie in the range: (0, +infinity).  If it is an output value, it will
            //    be searched for in the  range: [ 1.0D-300, 1.0D+300].
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
            double T11;
            double T12;
            double T13;
            double T14;
            double T15;

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
            //     F
            //
            if (!(f < 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     DFN
            //
            if (!(dfn <= 0.0e0)) goto S140;
            bound = 0.0e0;
            status = -5;
            return;
            S150:
            S140:
            if (which == 4) goto S170;
            //
            //     DFD
            //
            if (!(dfd <= 0.0e0)) goto S160;
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
            if (!(which == 1)) qporq = p <= q;
            //
            //     Select the minimum of P or Q
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //     Calculating P
                //
                cumf(f, dfn, dfd, ref p, ref q);
                status = 0;
            }
            else if (2 == which)
            {
                //
                //  Calculating F
                //
                f = 5.0e0;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                dstinv(K2, T3, K4, K4, K5, T6, T7);
                status = 0;
                dinvr(status, ref f, fx, ref qleft, ref qhi);
                S220:
                if (!(status == 1)) goto S250;
                cumf(f, dfn, dfd, ref cum, ref ccum);
                if (!qporq) goto S230;
                fx = cum - p;
                goto S240;
                S230:
                fx = ccum - q;
                S240:
                dinvr(status, ref f, fx, ref qleft, ref qhi);
                goto S220;
                S250:
                if (!(status == -1)) goto S280;
                if (!qleft) goto S260;
                status = 1;
                bound = 0.0e0;
                goto S270;
                S260:
                status = 2;
                bound = inf;
                S280:
                S270: ;
            }
            //
            //  Calculate DFN.
            //
            //  Note that, in the original calculation, the lower bound for DFN was 0.
            //  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
            //  The lower bound was set to the more reasonable value of 1.
            //  JVB, 14 April 2007.
            //
            else if (3 == which)
            {

                T8 = 1.0;
                T9 = inf;
                T10 = atol;
                T11 = tol;
                dstinv(T8, T9, K4, K4, K5, T10, T11);

                status = 0;
                dfn = 5.0;
                fx = 0.0;

                dinvr(status, ref dfn, fx, ref qleft, ref qhi);

                while (status == 1)
                {
                    cumf(f, dfn, dfd, ref cum, ref ccum);

                    if (p <= q)
                    {
                        fx = cum - p;
                    }
                    else
                    {
                        fx = ccum - q;
                    }

                    dinvr(status, ref dfn, fx, ref qleft, ref qhi);
                }

                if (status == -1)
                {
                    if (qleft)
                    {
                        status = 1;
                        bound = 1.0;
                    }
                    else
                    {
                        status = 2;
                        bound = inf;
                    }
                }
            }
            //
            //  Calculate DFD.
            //
            //  Note that, in the original calculation, the lower bound for DFD was 0.
            //  Using DFD = 0 causes an error in CUMF when it calls BETA_INC.
            //  The lower bound was set to the more reasonable value of 1.
            //  JVB, 14 April 2007.
            //
            //
            else if (4 == which)
            {
                T12 = 1.0;
                T13 = inf;
                T14 = atol;
                T15 = tol;
                dstinv(T12, T13, K4, K4, K5, T14, T15);

                status = 0;
                dfd = 5.0;
                fx = 0.0;
                dinvr(status, ref dfd, fx, ref qleft, ref qhi);

                while (status == 1)
                {
                    cumf(f, dfn, dfd, ref cum, ref ccum);

                    if (p <= q)
                    {
                        fx = cum - p;
                    }
                    else
                    {
                        fx = ccum - q;
                    }

                    dinvr(status, ref dfd, fx, ref qleft, ref qhi);
                }

                if (status == -1)
                {
                    if (qleft)
                    {
                        status = 1;
                        bound = 1.0;
                    }
                    else
                    {
                        status = 2;
                        bound = inf;
                    }
                }
            }

            return;
        }

        public static void cdffnc(int which, ref double p, ref double q, ref double f, ref double dfn,
                ref double dfd, ref double phonc, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFFNC evaluates the CDF of the Noncentral F distribution.
            //
            //  Discussion:
            //
            //    This routine originally used 1.0E+300 as the upper bound for the
            //    interval in which many of the missing parameters are to be sought.
            //    Since the underlying rootfinder routine needs to evaluate the
            //    function at this point, it is no surprise that the program was
            //    experiencing overflows.  A less extravagant upper bound
            //    is being tried for now!
            //
            //    This routine calculates any one parameter of the Noncentral F distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The computation time required for this routine is proportional
            //    to the noncentrality parameter PNONC.  Very large values of
            //    this parameter can consume immense computer resources.  This is
            //    why the search range is bounded by 10,000.
            //
            //    The value of the cumulative noncentral F distribution is not
            //    necessarily monotone in either degree of freedom.  There thus
            //    may be two values that provide a given CDF value.  This routine
            //    assumes monotonicity and will find an arbitrary one of the two
            //    values.
            //
            //    The CDF of the noncentral F distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF [ NoncentralFRatioDistribution [ DFN, DFD, PNONC ], X ]
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
            //    1966, Formula 26.6.20.
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
            //    1: Calculate P and Q from F, DFN, DFD and PNONC;
            //    2: Calculate F from P, Q, DFN, DFD and PNONC;
            //    3: Calculate DFN from P, Q, F, DFD and PNONC;
            //    4: Calculate DFD from P, Q, F, DFN and PNONC;
            //    5: Calculate PNONC from P, Q, F, DFN and DFD.
            //
            //    Input/output, double *P, the integral from 0 to F of
            //    the noncentral F-density.  If P is an input value it should
            //    lie in the range [0,1) (Not including 1!).
            //
            //    Dummy, double *Q, is not used by this subroutine,
            //    and is only included for similarity with the other routines.
            //    Its input value is not checked.  If P is to be computed, the
            //    Q is set to 1 - P.
            //
            //    Input/output, double *F, the upper limit of integration
            //    of the noncentral F-density.  If this is an input value, it should
            //    lie in the range: [0, +infinity).  If it is an output value, it
            //    will be searched for in the range: [0,1.0D+30].
            //
            //    Input/output, double *DFN, the number of degrees of freedom
            //    of the numerator sum of squares.  If this is an input value, it should
            //    lie in the range: (0, +infinity).  If it is an output value, it will
            //    be searched for in the range: [ 1.0, 1.0D+30].
            //
            //    Input/output, double *DFD, the number of degrees of freedom
            //    of the denominator sum of squares.  If this is an input value, it should
            //    be in range: (0, +infinity).  If it is an output value, it will be
            //    searched for in the range [1.0, 1.0D+30].
            //
            //    Input/output, double *PNONC, the noncentrality parameter
            //    If this is an input value, it should be nonnegative.
            //    If it is an output value, it will be searched for in the range: [0,1.0D+4].
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
            double tent4 = 1.0e4;
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double one = (1.0e0 - 1.0e-16);
            double inf = 1.0e300;

            double ccum = 0;
            double cum = 0;
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
            double T14;
            double T15;
            double T16;
            double T17;

            status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 5)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 5.0e0;
            S20:
            status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //     P
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
            //     F
            //
            if (!(f < 0.0e0)) goto S80;
            bound = 0.0e0;
            status = -4;
            return;
            S90:
            S80:
            if (which == 3) goto S110;
            //
            //  DFN
            //
            if (!(dfn <= 0.0e0)) goto S100;
            bound = 0.0e0;
            status = -5;
            return;
            S110:
            S100:
            if (which == 4) goto S130;
            //
            //     DFD
            //
            if (!(dfd <= 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -6;
            return;
            S130:
            S120:
            if (which == 5) goto S150;
            //
            //     PHONC
            //
            if (!(phonc < 0.0e0)) goto S140;
            bound = 0.0e0;
            status = -7;
            return;
            S150:
            S140:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P
                //
                cumfnc(f, dfn, dfd, phonc, ref p, ref q);
                status = 0;
            }
            else if (2 == which)
            {
                //
                //  Calculating F
                //
                f = 5.0e0;
                T2 = inf;
                T5 = atol;
                T6 = tol;
                dstinv(K1, T2, K3, K3, K4, T5, T6);
                status = 0;
                dinvr(status, ref f, fx, ref qleft, ref qhi);
                S160:
                if (!(status == 1)) goto S170;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref f, fx, ref qleft, ref qhi);
                goto S160;
                S170:
                if (!(status == -1)) goto S200;
                if (!qleft) goto S180;
                status = 1;
                bound = 0.0e0;
                goto S190;
                S180:
                status = 2;
                bound = inf;
                S200:
                S190: ;
            }
            else if (3 == which)
            {
                //
                //  Calculating DFN
                //
                dfn = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                dstinv(T7, T8, K3, K3, K4, T9, T10);
                status = 0;
                dinvr(status, ref dfn, fx, ref qleft, ref qhi);
                S210:
                if (!(status == 1)) goto S220;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref dfn, fx, ref qleft, ref qhi);
                goto S210;
                S220:
                if (!(status == -1)) goto S250;
                if (!qleft) goto S230;
                status = 1;
                bound = zero;
                goto S240;
                S230:
                status = 2;
                bound = inf;
                S250:
                S240: ;
            }
            else if (4 == which)
            {
                //
                //     Calculating DFD
                //
                dfd = 5.0e0;
                T11 = zero;
                T12 = inf;
                T13 = atol;
                T14 = tol;
                dstinv(T11, T12, K3, K3, K4, T13, T14);
                status = 0;
                dinvr(status, ref dfd, fx, ref qleft, ref qhi);
                S260:
                if (!(status == 1)) goto S270;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref dfd, fx, ref qleft, ref qhi);
                goto S260;
                S270:
                if (!(status == -1)) goto S300;
                if (!qleft) goto S280;
                status = 1;
                bound = zero;
                goto S290;
                S280:
                status = 2;
                bound = inf;
                S300:
                S290: ;
            }
            else if (5 == which)
            {
                //
                //     Calculating PHONC
                //
                phonc = 5.0e0;
                T15 = tent4;
                T16 = atol;
                T17 = tol;
                dstinv(K1, T15, K3, K3, K4, T16, T17);
                status = 0;
                dinvr(status, ref phonc, fx, ref qleft, ref qhi);
                S310:
                if (!(status == 1)) goto S320;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                dinvr(status, ref phonc, fx, ref qleft, ref qhi);
                goto S310;
                S320:
                if (!(status == -1)) goto S350;
                if (!qleft) goto S330;
                status = 1;
                bound = 0.0e0;
                goto S340;
                S330:
                status = 2;
                bound = tent4;
                S340: ;
            }

            S350:
            return;
        }

    }
}