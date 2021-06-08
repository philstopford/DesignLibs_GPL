using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfchi(int which, ref double p, ref double q, ref double x_, ref double df,
                ref int status_, ref double bound)

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
            double cum= 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double porq= 0;
            double pq= 0;
            bool qporq = false;
            double T3= 0;
            double T6= 0;
            double T7= 0;
            double T8= 0;
            double T9= 0;
            double T10= 0;
            double T11= 0;

            E0000_E0001_Data data = new E0000_E0001_Data();

            data.status = 0;
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
            data.status = -1;
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
            //     X
            //
            if (!(data.x < 0.0e0)) goto S120;
            bound = 0.0e0;
            data.status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S140;
            bound = 0.0e0;
            data.status = -5;
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
            data.status = 3;
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
                data.status = 0;
                cumchi(data.x, df, p, q);
                if (porq > 1.5e0)
                {
                    data.status = 10;
                    return;
                }
            }
            else if (2 == which)
            {
                //
                //  Calculating X
                //
                data.x = 5.0e0;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                E0000E0001.dstinv(ref data, K2, T3, K4, K4, K5, T6, T7);
                data.status = 0;
                E0000E0001.dinvr(ref data);
                S230:
                if (!(data.status == 1)) goto S270;
                cumchi(data.x, df, cum, ccum);
                if (!qporq) goto S240;
                data.fx = cum - p;
                goto S250;
                S240:
                data.fx = ccum - q;
                S250:
                if (!(data.fx + porq > 1.5e0)) goto S260;
                data.status = 10;
                return;
                S260:
                E0000E0001.dinvr(ref data);
                goto S230;
                S270:
                if (!(data.status == -1)) goto S300;
                if (!data.qleft) goto S280;
                data.status = 1;
                bound = 0.0e0;
                goto S290;
                S280:
                data.status = 2;
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
                E0000E0001.dstinv(ref data, T8, T9, K4, K4, K5, T10, T11);
                data.status = 0;
                data.x = df;
                E0000E0001.dinvr(ref data);
                S310:
                if (!(data.status == 1)) goto S350;
                cumchi(data.x, df, cum, ccum);
                if (!qporq) goto S320;
                data.fx = cum - p;
                goto S330;
                S320:
                data.fx = ccum - q;
                S330:
                if (!(data.fx + porq > 1.5e0)) goto S340;
                data.status = 10;
                return;
                S340:
                data.x = df;
                E0000E0001.dinvr(ref data);
                goto S310;
                S350:
                if (!(data.status == -1)) goto S380;
                if (!data.qleft) goto S360;
                data.status = 1;
                bound = zero;
                goto S370;
                S360:
                data.status = 2;
                bound = inf;
                S370: ;
            }

            S380:
            x_ = data.x;
            status_ = data.status;
            return;
        }

    }
}