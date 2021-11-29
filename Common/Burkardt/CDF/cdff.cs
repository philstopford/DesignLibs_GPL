﻿using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdff(int which, ref double p, ref double q, ref double f, ref double dfn,
            ref double dfd, ref int status_, ref double bound)

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
        const double tol = 1.0e-8;
        const double atol = 1.0e-50;
        const double inf = 1.0e300;

        double ccum = 0;
        double cum = 0;
        const int K1 = 1;
        const double K2 = 0.0e0;
        const double K4 = 0.5e0;
        const double K5 = 5.0e0;
        double pq = 0;
        bool qporq = false;

        E0000_E0001_Data data = new()
        {
            status = 0
        };

        bound = 0.0;
        switch (which is < 1 or > 4)
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

        switch (p is < 0.0e0 or > 1.0e0)
        {
            //
            //  P
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

        switch (f < 0.0e0)
        {
            //
            //     F
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

        switch (dfn <= 0.0e0)
        {
            //
            //     DFN
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
            case 4:
                goto S170;
        }

        switch (dfd <= 0.0e0)
        {
            //
            //     DFD
            //
            case false:
                goto S160;
        }

        bound = 0.0e0;
        data.status = -6;
        return;
        S170:
        S160:
        switch (which)
        {
            case 1:
                goto S210;
        }

        //
        //     P + Q
        //
        pq = p + q;
        if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1)))
        {
            goto S200;
        }

        switch (pq < 0.0e0)
        {
            case false:
                goto S180;
        }

        bound = 0.0e0;
        goto S190;
        S180:
        bound = 1.0e0;
        S190:
        data.status = 3;
        return;
        S210:
        S200:
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
                //     Calculating P
                //
                cumf(f, dfn, dfd, ref p, ref q);
                data.status = 0;
                break;
            case 2:
            {
                //
                //  Calculating F
                //
                f = 5.0e0;
                E0000E0001.dstinv(ref data, K2, inf, K4, K4, K5, atol, tol);
                data.status = 0;
                data.x = f;
                E0000E0001.dinvr(ref data);
                S220:
                switch (data.status == 1)
                {
                    case false:
                        goto S250;
                }

                cumf(f, dfn, dfd, ref cum, ref ccum);
                switch (qporq)
                {
                    case false:
                        goto S230;
                }

                data.fx = cum - p;
                goto S240;
                S230:
                data.fx = ccum - q;
                S240:
                data.x = f;
                E0000E0001.dinvr(ref data);
                goto S220;
                S250:
                switch (data.status == -1)
                {
                    case false:
                        goto S280;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S260;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S270;
                S260:
                data.status = 2;
                bound = inf;
                S280:
                S270: ;
                break;
            }
            //
            //  Calculate DFN.
            //
            //  Note that, in the original calculation, the lower bound for DFN was 0.
            //  Using DFN = 0 causes an error in CUMF when it calls BETA_INC.
            //  The lower bound was set to the more reasonable value of 1.
            //  JVB, 14 April 2007.
            //
            case 3:
            {
                const double T8 = 1.0;
                E0000E0001.dstinv(ref data, T8, inf, K4, K4, K5, atol, tol);

                data.status = 0;
                dfn = 5.0;
                data.fx = 0.0;
                data.x = dfn;

                E0000E0001.dinvr(ref data);

                while (data.status == 1)
                {
                    cumf(f, dfn, dfd, ref cum, ref ccum);

                    if (p <= q)
                    {
                        data.fx = cum - p;
                    }
                    else
                    {
                        data.fx = ccum - q;
                    }

                    data.x = dfn;
                    E0000E0001.dinvr(ref data);
                }

                switch (data.status)
                {
                    case -1 when data.qleft:
                        data.status = 1;
                        bound = 1.0;
                        break;
                    case -1:
                        data.status = 2;
                        bound = inf;
                        break;
                }

                break;
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
            case 4:
            {
                const double T12 = 1.0;
                E0000E0001.dstinv(ref data, T12, inf, K4, K4, K5, atol, tol);

                data.status = 0;
                dfd = 5.0;
                data.fx = 0.0;
                data.x = dfd;
                E0000E0001.dinvr(ref data);

                while (data.status == 1)
                {
                    cumf(f, dfn, dfd, ref cum, ref ccum);

                    if (p <= q)
                    {
                        data.fx = cum - p;
                    }
                    else
                    {
                        data.fx = ccum - q;
                    }

                    data.x = dfd;
                    E0000E0001.dinvr(ref data);
                }

                switch (data.status)
                {
                    case -1 when data.qleft:
                        data.status = 1;
                        bound = 1.0;
                        break;
                    case -1:
                        data.status = 2;
                        bound = inf;
                        break;
                }

                break;
            }
        }

        status_ = data.status;
    }
}