using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdft(int which, ref double p, ref double q, ref double t, ref double df,
            ref int status_, ref double bound)

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
        const double tol = 1.0e-8;
        const double atol = 1.0e-50;
        const double zero = 1.0e-300;
        const double inf = 1.0e30;
        const double maxdf = 1.0e10;

        double ccum = 0;
        double cum = 0;
        const int K1 = 1;
        const double K4 = 0.5e0;
        const double K5 = 5.0e0;
        bool qporq = false;

        E0000_E0001_Data data = new()
        {
            status = 0
        };

        bound = 0.0;
        switch (which is < 1 or > 3)
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
        bound = 3.0e0;
        S20:
        data.status = -1;
        return;
        S30:
        switch (which)
        {
            case 1:
                goto S70;
        }

        switch (p is <= 0.0e0 or > 1.0e0)
        {
            //
            //     P
            //
            case false:
                goto S60;
        }

        switch (p <= 0.0e0)
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
            case 3:
                goto S130;
        }

        switch (df <= 0.0e0)
        {
            //
            //     DF
            //
            case false:
                goto S120;
        }

        bound = 0.0e0;
        data.status = -5;
        return;
        S130:
        S120:
        switch (which)
        {
            case 1:
                goto S170;
        }

        //
        //     P + Q
        //
        double pq = p + q;
        if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S160;
        }

        switch (pq < 0.0e0)
        {
            case false:
                goto S140;
        }

        bound = 0.0e0;
        goto S150;
        S140:
        bound = 1.0e0;
        S150:
        data.status = 3;
        return;
        S170:
        S160:
        qporq = (which == 1) switch
        {
            false => p <= q,
            _ => qporq
        };

        switch (which)
        {
            //
            //  Select the minimum of P or Q.  Calculate ANSWERS
            //
            case 1:
                //
                //  Computing P and Q
                //
                cumt(t, df, ref p, ref q);
                data.status = 0;
                break;
            case 2:
            {
                //
                //  Computing T
                //  Get initial approximation for T
                //
                t = dt1(p, q, df);
                double T2 = -inf;
                double T3 = inf;
                double T6 = atol;
                double T7 = tol;
                E0000E0001.dstinv(ref data, T2, T3, K4, K4, K5, T6, T7);
                data.status = 0;
                data.x = t;
                E0000E0001.dinvr(ref data);
                S180:
                switch (data.status == 1)
                {
                    case false:
                        goto S210;
                }

                cumt(t, df, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S190;
                }

                data.fx = cum - p;
                goto S200;
                S190:
                data.fx = ccum - q;
                S200:
                data.x = t;
                E0000E0001.dinvr(ref data);
                goto S180;
                S210:
                switch (data.status == -1)
                {
                    case false:
                        goto S240;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S220;
                }

                data.status = 1;
                bound = -inf;
                goto S230;
                S220:
                data.status = 2;
                bound = inf;
                S240:
                S230: ;
                break;
            }
            case 3:
            {
                //
                //  Computing DF
                //
                df = 5.0e0;
                double T8 = zero;
                double T9 = maxdf;
                double T10 = atol;
                double T11 = tol;
                E0000E0001.dstinv(ref data, T8, T9, K4, K4, K5, T10, T11);
                data.status = 0;
                data.x = df;
                E0000E0001.dinvr(ref data);
                S250:
                switch (data.status == 1)
                {
                    case false:
                        goto S280;
                }

                cumt(t, df, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S260;
                }

                data.fx = cum - p;
                goto S270;
                S260:
                data.fx = ccum - q;
                S270:
                data.x = df;
                E0000E0001.dinvr(ref data);
                goto S250;
                S280:
                switch (data.status == -1)
                {
                    case false:
                        goto S310;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S290;
                }

                data.status = 1;
                bound = zero;
                goto S300;
                S290:
                data.status = 2;
                bound = maxdf;
                S300: ;
                break;
            }
        }

        S310:
        status_ = data.status;
    }
}