using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdfpoi(int which, ref double p, ref double q, ref double s, ref double xlam,
            ref int status_, ref double bound)

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
        const double tol = 1.0e-8;
        const double atol = 1.0e-50;
        const double inf = 1.0e300;

        double ccum = 0;
        double cum = 0;
        const int K1 = 1;
        const double K2 = 0.0e0;
        const double K4 = 0.5e0;
        const double K5 = 5.0e0;
        bool qporq = false;

        E0000_E0001_Data data = new()
        {
            status = 0
        };

        bound = 0.0;
        switch (which is >= 1 and <= 3)
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

        switch (xlam < 0.0e0)
        {
            //
            //     XLAM
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
            case 1:
                goto S190;
        }

        //
        //     P + Q
        //
        double pq = p + q;
        if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S180;
        }

        switch (pq < 0.0e0)
        {
            case false:
                goto S160;
        }

        bound = 0.0e0;
        goto S170;
        S160:
        bound = 1.0e0;
        S170:
        data.status = 3;
        return;
        S190:
        S180:
        qporq = (which == 1) switch
        {
            false => p <= q,
            _ => qporq
        };

        switch (which)
        {
            //
            //  Select the minimum of P or Q
            //  Calculate ANSWERS
            //
            case 1:
                //
                //  Calculating P
                //
                cumpoi(s, xlam, ref p, ref q);
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
                S200:
                switch (data.status == 1)
                {
                    case false:
                        goto S230;
                }

                cumpoi(s, xlam, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S210;
                }

                data.fx = cum - p;
                goto S220;
                S210:
                data.fx = ccum - q;
                S220:
                data.x = s;
                E0000E0001.dinvr(ref data);
                goto S200;
                S230:
                switch (data.status == -1)
                {
                    case false:
                        goto S260;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S240;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S250;
                S240:
                data.status = 2;
                bound = inf;
                S260:
                S250: ;
                break;
            }
            case 3:
            {
                //
                //     Calculating XLAM
                //
                xlam = 5.0e0;
                E0000E0001.dstinv(ref data, K2, inf, K4, K4, K5, atol, tol);
                data.status = 0;
                data.x = xlam;
                E0000E0001.dinvr(ref data);
                S270:
                switch (data.status == 1)
                {
                    case false:
                        goto S300;
                }

                cumpoi(s, xlam, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S280;
                }

                data.fx = cum - p;
                goto S290;
                S280:
                data.fx = ccum - q;
                S290:
                data.x = xlam;
                E0000E0001.dinvr(ref data);
                goto S270;
                S300:
                switch (data.status == -1)
                {
                    case false:
                        goto S330;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S310;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S320;
                S310:
                data.status = 2;
                bound = inf;
                S320: ;
                break;
            }
        }

        S330:
        status_ = data.status;
    }

}