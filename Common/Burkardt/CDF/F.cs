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

            E0000Data e0000Data = new E0000Data();
            E0001Data e0001Data = new E0001Data();
            
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
                e0000Data.zsmall = K2;
                e0000Data.zbig = T3;
                e0000Data.zabsst = K4;
                e0000Data.zrelst = K4;
                e0000Data.zstpmu = K5;
                e0000Data.zabsto = T6;
                e0000Data.zrelto = T7;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S220:
                if (!(status == 1)) goto S250;
                cumf(f, dfn, dfd, ref cum, ref ccum);
                if (!qporq) goto S230;
                fx = cum - p;
                goto S240;
                S230:
                fx = ccum - q;
                S240:
                e0000Data.status = status;
                dinvr(ref e0000Data);
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
                e0000Data.zsmall = T8;
                e0000Data.zbig = T9;
                e0000Data.zabsst = K4;
                e0000Data.zrelst = K4;
                e0000Data.zstpmu = K5;
                e0000Data.zabsto = T10;
                e0000Data.zrelto = T11;
                dstinv(ref e0000Data);

                status = 0;
                dfn = 5.0;
                fx = 0.0;

                e0000Data.status = status;
                dinvr(ref e0000Data);

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

                    e0000Data.status = status;
                    dinvr(ref e0000Data);
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
                e0000Data.zsmall = T12;
                e0000Data.zbig = T13;
                e0000Data.zabsst = K4;
                e0000Data.zrelst = K4;
                e0000Data.zstpmu = K5;
                e0000Data.zabsto = T14;
                e0000Data.zrelto = T15;
                dstinv(ref e0000Data);

                status = 0;
                dfd = 5.0;
                fx = 0.0;
                e0000Data.status = status;
                dinvr(ref e0000Data);

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

                    e0000Data.status = status;
                    dinvr(ref e0000Data);
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

            E0000Data e0000Data = new E0000Data();
            E0001Data e0001Data = new E0001Data();
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
                e0000Data.zsmall = K1;
                e0000Data.zbig = T2;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T5;
                e0000Data.zrelto = T6;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S160:
                if (!(status == 1)) goto S170;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                e0000Data.status = status;
                dinvr(ref e0000Data);
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
                e0000Data.zsmall = T7;
                e0000Data.zbig = T8;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T9;
                e0000Data.zrelto = T10;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S210:
                if (!(status == 1)) goto S220;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                e0000Data.status = status;
                dinvr(ref e0000Data);
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
                e0000Data.zsmall = T11;
                e0000Data.zbig = T12;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T13;
                e0000Data.zrelto = T14;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S260:
                if (!(status == 1)) goto S270;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                e0000Data.status = status;
                dinvr(ref e0000Data);
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
                e0000Data.zsmall = K1;
                e0000Data.zbig = T15;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T16;
                e0000Data.zrelto = T17;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S310:
                if (!(status == 1)) goto S320;
                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                fx = cum - p;
                e0000Data.status = status;
                dinvr(ref e0000Data);
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

        public static void f_noncentral_cdf_values(ref int n_data, ref int a, ref int b, ref double lambda,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
            //
            //  Discussion:
            //
            //    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
            //    in Mathematica by commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
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
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz and Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    US Department of Commerce, 1964.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int *A, int *B, double *LAMBDA, the
            //    parameters of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 22;

            int[] a_vec =
                {
                    1, 1, 1, 1,
                    1, 1, 1, 1,
                    1, 1, 2, 2,
                    3, 3, 4, 4,
                    5, 5, 6, 6,
                    8, 16
                }
                ;
            int[] b_vec =
                {
                    1, 5, 5, 5,
                    5, 5, 5, 5,
                    5, 5, 5, 10,
                    5, 5, 5, 5,
                    1, 5, 6, 12,
                    16, 8
                }
                ;
            double[] fx_vec =
                {
                    0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
                    0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
                    0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
                    0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
                    0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
                    0.844589E+00, 0.816368E+00
                }
                ;
            double[] lambda_vec =
                {
                    0.00E+00, 0.000E+00, 0.25E+00, 1.00E+00,
                    1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                    1.00E+00, 2.00E+00, 1.00E+00, 1.00E+00,
                    1.00E+00, 2.00E+00, 1.00E+00, 1.00E+00,
                    0.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                    1.00E+00, 1.00E+00
                }
                ;
            double[] x_vec =
                {
                    1.00E+00, 1.00E+00, 1.00E+00, 0.50E+00,
                    1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
                    5.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                    1.00E+00, 1.00E+00, 1.00E+00, 2.00E+00,
                    1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                    2.00E+00, 2.00E+00
                }
                ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0;
                lambda = 0.0E+00;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void f_cdf_values(ref int n_data, ref int a, ref int b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_CDF_VALUES returns some values of the F CDF test function.
            //
            //  Discussion:
            //
            //    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
            //    commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[FRatioDistribution[ DFN, DFD ], X ]
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
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz and Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    US Department of Commerce, 1964.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int *A, int *B, the parameters of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 20;

            int[] a_vec =  {
                1, 1, 5, 1,
                2, 4, 1, 6,
                8, 1, 3, 6,
                1, 1, 1, 1,
                2, 3, 4, 5
            }
            ;
            int[] b_vec =  {
                1, 5, 1, 5,
                10, 20, 5, 6,
                16, 5, 10, 12,
                5, 5, 5, 5,
                5, 5, 5, 5
            }
            ;
            double[] fx_vec =  {
                0.500000E+00, 0.499971E+00, 0.499603E+00, 0.749699E+00,
                0.750466E+00, 0.751416E+00, 0.899987E+00, 0.899713E+00,
                0.900285E+00, 0.950025E+00, 0.950057E+00, 0.950193E+00,
                0.975013E+00, 0.990002E+00, 0.994998E+00, 0.999000E+00,
                0.568799E+00, 0.535145E+00, 0.514343E+00, 0.500000E+00
            }
            ;
            double[] x_vec =  {
                1.00E+00, 0.528E+00, 1.89E+00, 1.69E+00,
                1.60E+00, 1.47E+00, 4.06E+00, 3.05E+00,
                2.09E+00, 6.61E+00, 3.71E+00, 3.00E+00,
                10.01E+00, 16.26E+00, 22.78E+00, 47.18E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }


    }
}