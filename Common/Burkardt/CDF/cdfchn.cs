namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfchn(int which, ref double p, ref double q, ref double x_, ref double df,
                ref double pnonc, ref int status_, ref double bound)

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
            double one = (1.0e0-1.0e-16);
            double inf = 1.0e300;

            double cum = 0;
            double ccum = 0;
            double K1 = 0.0e0;
            double K3 = 0.5e0;
            double K4 = 5.0e0;
            double T2 = 0;
            double T5 = 0;
            double T6 = 0;
            double T7 = 0;
            double T8 = 0;
            double T9 = 0;
            double T10 = 0;
            double T11 = 0;
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
            //  P
            //
            if (!(p < 0.0e0 || p > one)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = one;
            S50:
            data.status = -2;
            return;
            S70:
            S60:
            if (which == 2) goto S90;
            //
            //     X
            //
            if (!(data.x < 0.0e0)) goto S80;
            bound = 0.0e0;
            data.status = -4;
            return;
            S90:
            S80:
            if (which == 3) goto S110;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S100;
            bound = 0.0e0;
            data.status = -5;
            return;
            S110:
            S100:
            if (which == 4) goto S130;
            //
            //  PNONC
            //
            if (!(pnonc < 0.0e0)) goto S120;
            bound = 0.0e0;
            data.status = -6;
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
                cumchn(data.x, df, pnonc, ref p, ref q);
                data.status = 0;
            }
            else if (2 == which)
            {
                //
                //     Calculating X
                //
                data.x = 5.0e0;
                T2 = inf;
                T5 = atol;
                T6 = tol;
                E0000E0001.dstinv(ref data, K1, T2, K3, K3, K4, T5, T6);
                data.status = 0;
                E0000E0001.dinvr(ref data);
                S140:
                if (!(data.status == 1)) goto S150;
                cumchn(data.x, df, pnonc, ref cum, ref ccum);
                data.fx = cum - p;
                E0000E0001.dinvr(ref data);
                goto S140;
                S150:
                if (!(data.status == -1)) goto S180;
                if (!data.qleft) goto S160;
                data.status = 1;
                bound = 0.0e0;
                goto S170;
                S160:
                data.status = 2;
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
                E0000E0001.dstinv(ref data, T7, T8, K3, K3, K4, T9, T10);
                data.status = 0;
                data.x = df;
                E0000E0001.dinvr(ref data);
                S190:
                if (!(data.status == 1)) goto S200;
                cumchn(data.x, df, pnonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = df;
                E0000E0001.dinvr(ref data);
                goto S190;
                S200:
                if (!(data.status == -1)) goto S230;
                if (!data.qleft) goto S210;
                data.status = 1;
                bound = zero;
                goto S220;
                S210:
                data.status = 2;
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
                E0000E0001.dstinv(ref data, K1, T11, K3, K3, K4, T12, T13);
                data.status = 0;
                data.x = pnonc;
                E0000E0001.dinvr(ref data);
                S240:
                if (!(data.status == 1)) goto S250;
                cumchn(data.x, df, pnonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = pnonc;
                E0000E0001.dinvr(ref data);
                goto S240;
                S250:
                if (!(data.status == -1)) goto S280;
                if (!data.qleft) goto S260;
                data.status = 1;
                bound = zero;
                goto S270;
                S260:
                data.status = 2;
                bound = tent4;
                S270: ;
            }

            S280:
            x_ = data.x;
            status_ = data.status;
        }
    }
}