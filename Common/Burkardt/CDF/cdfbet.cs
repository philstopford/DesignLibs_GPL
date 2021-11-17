using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdfbet(int which, ref double p, ref double q, ref double x_, ref double y,
            ref double a, ref double b, ref int status_, ref double bound)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CDFBET evaluates the CDF of the Beta Distribution.
        //
        //  Discussion:
        //
        //    This routine calculates any one parameter of the beta distribution
        //    given the others.
        //
        //    The value P of the cumulative distribution function is calculated
        //    directly by code associated with the reference.
        //
        //    Computation of the other parameters involves a seach for a value that
        //    produces the desired value of P.  The search relies on the
        //    monotonicity of P with respect to the other parameters.
        //
        //    The beta density is proportional to t^(A-1) * (1-t)^(B-1).
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
        //    Algorithm 708:
        //    Significant Digit Computation of the Incomplete Beta Function Ratios,
        //    ACM Transactions on Mathematical Software,
        //    Volume 18, 1993, pages 360-373.
        //
        //  Parameters:
        //
        //    Input, int *WHICH, indicates which of the next four argument
        //    values is to be calculated from the others.
        //    1: Calculate P and Q from X, Y, A and B;
        //    2: Calculate X and Y from P, Q, A and B;
        //    3: Calculate A from P, Q, X, Y and B;
        //    4: Calculate B from P, Q, X, Y and A.
        //
        //    Input/output, double *P, the integral from 0 to X of the
        //    chi-square distribution.  Input range: [0, 1].
        //
        //    Input/output, double *Q, equals 1-P.  Input range: [0, 1].
        //
        //    Input/output, double *X, the upper limit of integration
        //    of the beta density.  If it is an input value, it should lie in
        //    the range [0,1].  If it is an output value, it will be searched for
        //    in the range [0,1].
        //
        //    Input/output, double *Y, equal to 1-X.  If it is an input
        //    value, it should lie in the range [0,1].  If it is an output value,
        //    it will be searched for in the range [0,1].
        //
        //    Input/output, double *A, the first parameter of the beta
        //    density.  If it is an input value, it should lie in the range
        //    (0, +infinity).  If it is an output value, it will be searched
        //    for in the range [1D-300,1D300].
        //
        //    Input/output, double *B, the second parameter of the beta
        //    density.  If it is an input value, it should lie in the range
        //    (0, +infinity).  If it is an output value, it will be searched
        //    for in the range [1D-300,1D300].
        //
        //    Output, int *STATUS, reports the status of the computation.
        //     0, if the calculation completed correctly;
        //    -I, if the input parameter number I is out of range;
        //    +1, if the answer appears to be lower than lowest search bound;
        //    +2, if the answer appears to be higher than greatest search bound;
        //    +3, if P + Q /= 1;
        //    +4, if X + Y /= 1.
        //
        //    Output, double *BOUND, is only defined if STATUS is nonzero.
        //    If STATUS is negative, then this is the value exceeded by parameter I.
        //    if STATUS is 1 or 2, this is the search bound that was exceeded.
        //
    {
        double tol = 1.0e-8;
        double atol = 1.0e-50;
        double zero = 1.0e-300;
        double inf = 1.0e300;
        double one = 1.0e0;

        double ccum = 0;
        double cum = 0;
        int K1 = 1;
        double K2 = 0.0e0;
        double K3 = 1.0e0;
        double K8 = 0.5e0;
        double K9 = 5.0e0;
        double pq = 0;
        bool qporq = false;
        double xy = 0;

        E0000_E0001_Data data = new();

        double T4 = 0, T5 = 0, T6 = 0, T7 = 0, T10 = 0, T11 = 0, T12 = 0, T13 = 0, T14 = 0, T15 = 0;

        data.status = 0;
        data.x = x_;
        bound = 0.0;
        switch ((which < 1 || which > 4))
        {
            //
            //  Check arguments
            //
            case false:
                goto S30;
        }

        switch ((which < 1))
        {
            case false:
                goto S10;
        }

        bound = 1.0e0;
        goto S20;
        S10:
        bound = 4.0e0;
        S20:
        data.status = -1;
        return;
        S30:
        switch (which)
        {
            case 1:
                goto S70;
        }

        switch ((p < 0.0e0 || p > 1.0e0))
        {
            //
            //  P
            //
            case false:
                goto S60;
        }

        switch ((p < 0.0e0))
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

        switch ((q < 0.0e0 || q > 1.0e0))
        {
            //
            //  Q
            //
            case false:
                goto S100;
        }

        switch ((q < 0.0e0))
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
                goto S150;
        }

        switch ((data.x < 0.0e0 || data.x > 1.0e0))
        {
            //
            //  X
            //
            case false:
                goto S140;
        }

        switch ((data.x < 0.0e0))
        {
            case false:
                goto S120;
        }

        bound = 0.0e0;
        goto S130;
        S120:
        bound = 1.0e0;
        S130:
        data.status = -4;
        return;
        S150:
        S140:
        switch (which)
        {
            case 2:
                goto S190;
        }

        switch ((y < 0.0e0 || y > 1.0e0))
        {
            //
            //  Y
            //
            case false:
                goto S180;
        }

        switch ((y < 0.0e0))
        {
            case false:
                goto S160;
        }

        bound = 0.0e0;
        goto S170;
        S160:
        bound = 1.0e0;
        S170:
        data.status = -5;
        return;
        S190:
        S180:
        switch (which)
        {
            case 3:
                goto S210;
        }

        switch ((a <= 0.0e0))
        {
            //
            //  A
            //
            case false:
                goto S200;
        }

        bound = 0.0e0;
        data.status = -6;
        return;
        S210:
        S200:
        switch (which)
        {
            case 4:
                goto S230;
        }

        switch ((b <= 0.0e0))
        {
            //
            //  B
            //
            case false:
                goto S220;
        }

        bound = 0.0e0;
        data.status = -7;
        return;
        S230:
        S220:
        switch (which)
        {
            case 1:
                goto S270;
        }

        //
        //  P + Q
        //
        pq = p + q;
        if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S260;
        }

        switch ((pq < 0.0e0))
        {
            case false:
                goto S240;
        }

        bound = 0.0e0;
        goto S250;
        S240:
        bound = 1.0e0;
        S250:
        data.status = 3;
        return;
        S270:
        S260:
        switch (which)
        {
            case 2:
                goto S310;
        }

        //
        //  X + Y
        //
        xy = data.x + y;
        if (!(Math.Abs(xy - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S300;
        }

        switch ((xy < 0.0e0))
        {
            case false:
                goto S280;
        }

        bound = 0.0e0;
        goto S290;
        S280:
        bound = 1.0e0;
        S290:
        data.status = 4;
        return;
        S310:
        S300:
        qporq = (which == 1) switch
        {
            false => p <= q,
            _ => qporq
        };

        switch (which)
        {
            //
            //     Select the minimum of P or Q
            //     Calculate ANSWERS
            //
            case 1:
                //
                //  Calculating P and Q
                //
                cumbet(data.x, y, a, b, ref p, ref q);
                data.status = 0;
                break;
            case 2:
            {
                //
                //  Calculating X and Y
                //
                T4 = atol;
                T5 = tol;
                E0000E0001.dstzr(ref data, K2, K3, T4, T5);
                switch (qporq)
                {
                    case false:
                        goto S340;
                }

                data.status = 0;
                E0000E0001.dzror(ref data);
                y = one - data.x;
                S320:
                switch ((data.status == 1))
                {
                    case false:
                        goto S330;
                }

                cumbet(data.x, y, a, b, ref cum, ref ccum);
                data.fx = cum - p;
                E0000E0001.dzror(ref data);
                y = one - data.x;
                goto S320;
                S330:
                goto S370;
                S340:
                data.status = 0;
                data.x = y;
                E0000E0001.dzror(ref data);
                data.x = one - data.x;
                S350:
                switch ((data.status == 1))
                {
                    case false:
                        goto S360;
                }

                cumbet(data.x, y, a, b, ref cum, ref ccum);
                data.fx = ccum - q;
                data.x = y;
                E0000E0001.dzror(ref data);
                data.x = one - data.x;
                goto S350;
                S370:
                S360:
                switch ((data.status == -1))
                {
                    case false:
                        goto S400;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S380;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S390;
                S380:
                data.status = 2;
                bound = 1.0e0;
                S400:
                S390: ;
                break;
            }
            case 3:
            {
                //
                //  Computing A
                //
                a = 5.0e0;
                T6 = zero;
                T7 = inf;
                T10 = atol;
                T11 = tol;
                E0000E0001.dstinv(ref data, T6, T7, K8, K8, K9, T10, T11);
                data.status = 0;
                data.x = a;
                E0000E0001.dinvr(ref data);
                S410:
                switch ((data.status == 1))
                {
                    case false:
                        goto S440;
                }

                cumbet(data.x, y, a, b, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S420;
                }

                data.fx = cum - p;
                goto S430;
                S420:
                data.fx = ccum - q;
                S430:
                data.x = a;
                E0000E0001.dinvr(ref data);
                goto S410;
                S440:
                switch ((data.status == -1))
                {
                    case false:
                        goto S470;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S450;
                }

                data.status = 1;
                bound = zero;
                goto S460;
                S450:
                data.status = 2;
                bound = inf;
                S470:
                S460: ;
                break;
            }
            case 4:
            {
                //
                //  Computing B
                //
                b = 5.0e0;
                T12 = zero;
                T13 = inf;
                T14 = atol;
                T15 = tol;
                E0000E0001.dstinv(ref data, T12, T13, K8, K8, K9, T14, T15);
                data.status = 0;
                data.x = b;
                E0000E0001.dinvr(ref data);
                S480:
                switch ((data.status == 1))
                {
                    case false:
                        goto S510;
                }

                cumbet(data.x, y, a, b, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S490;
                }

                data.fx = cum - p;
                goto S500;
                S490:
                data.fx = ccum - q;
                S500:
                data.x = b;
                E0000E0001.dinvr(ref data);
                goto S480;
                S510:
                switch ((data.status == -1))
                {
                    case false:
                        goto S540;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S520;
                }

                data.status = 1;
                bound = zero;
                goto S530;
                S520:
                data.status = 2;
                bound = inf;
                S530: ;
                break;
            }
        }

        S540:
        x_ = data.x;
        y = 1 - x_;
        status_ = data.status;
    }
}