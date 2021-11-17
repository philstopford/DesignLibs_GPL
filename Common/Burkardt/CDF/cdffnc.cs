namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cdffnc(int which, ref double p, ref double q, ref double f, ref double dfn,
            ref double dfd, ref double phonc, ref int status_, ref double bound)

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
        double tol = 1.0e-8;
        double atol = 1.0e-50;
        double zero = 1.0e-300;
        double one = 1.0e0-1.0e-16;
        double inf = 1.0e300;

        double ccum = 0;
        double cum = 0;
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
        double T14 = 0;
        double T15 = 0;
        double T16 = 0;
        double T17 = 0;
            
        E0000_E0001_Data data = new()
        {
            status = 0
        };

        bound = 0.0;
        switch ((which < 1 || which > 5))
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
        bound = 5.0e0;
        S20:
        data.status = -1;
        return;
        S30:
        switch (which)
        {
            case 1:
                goto S70;
        }

        switch ((p < 0.0e0 || p > one))
        {
            //
            //     P
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
        bound = one;
        S50:
        data.status = -2;
        return;
        S70:
        S60:
        switch (which)
        {
            case 2:
                goto S90;
        }

        switch ((f < 0.0e0))
        {
            //
            //     F
            //
            case false:
                goto S80;
        }

        bound = 0.0e0;
        data.status = -4;
        return;
        S90:
        S80:
        switch (which)
        {
            case 3:
                goto S110;
        }

        switch ((dfn <= 0.0e0))
        {
            //
            //  DFN
            //
            case false:
                goto S100;
        }

        bound = 0.0e0;
        data.status = -5;
        return;
        S110:
        S100:
        switch (which)
        {
            case 4:
                goto S130;
        }

        switch ((dfd <= 0.0e0))
        {
            //
            //     DFD
            //
            case false:
                goto S120;
        }

        bound = 0.0e0;
        data.status = -6;
        return;
        S130:
        S120:
        switch (which)
        {
            case 5:
                goto S150;
        }

        switch ((phonc < 0.0e0))
        {
            //
            //     PHONC
            //
            case false:
                goto S140;
        }

        bound = 0.0e0;
        data.status = -7;
        return;
        S150:
        S140:
        switch (which)
        {
            //
            //     Calculate ANSWERS
            //
            case 1:
                //
                //  Calculating P
                //
                cumfnc(f, dfn, dfd, phonc, ref p, ref q);
                data.status = 0;
                break;
            case 2:
            {
                //
                //  Calculating F
                //
                f = 5.0e0;
                T2 = inf;
                T5 = atol;
                T6 = tol;
                E0000E0001.dstinv(ref data, K1, T2, K3, K3, K4, T5, T6);
                data.status = 0;
                data.x = f;
                E0000E0001.dinvr(ref data);
                S160:
                switch ((data.status == 1))
                {
                    case false:
                        goto S170;
                }

                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = f;
                E0000E0001.dinvr(ref data);
                goto S160;
                S170:
                switch ((data.status == -1))
                {
                    case false:
                        goto S200;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S180;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S190;
                S180:
                data.status = 2;
                bound = inf;
                S200:
                S190: ;
                break;
            }
            case 3:
            {
                //
                //  Calculating DFN
                //
                dfn = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                E0000E0001.dstinv(ref data, T7, T8, K3, K3, K4, T9, T10);
                data.status = 0;
                data.x = dfn;
                E0000E0001.dinvr(ref data);
                S210:
                switch ((data.status == 1))
                {
                    case false:
                        goto S220;
                }

                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = dfn;
                E0000E0001.dinvr(ref data);
                goto S210;
                S220:
                switch ((data.status == -1))
                {
                    case false:
                        goto S250;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S230;
                }

                data.status = 1;
                bound = zero;
                goto S240;
                S230:
                data.status = 2;
                bound = inf;
                S250:
                S240: ;
                break;
            }
            case 4:
            {
                //
                //     Calculating DFD
                //
                dfd = 5.0e0;
                T11 = zero;
                T12 = inf;
                T13 = atol;
                T14 = tol;
                E0000E0001.dstinv(ref data, T11, T12, K3, K3, K4, T13, T14);
                data.status = 0;
                data.x = dfd;
                E0000E0001.dinvr(ref data);
                S260:
                switch ((data.status == 1))
                {
                    case false:
                        goto S270;
                }

                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = dfd;
                E0000E0001.dinvr(ref data);
                goto S260;
                S270:
                switch ((data.status == -1))
                {
                    case false:
                        goto S300;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S280;
                }

                data.status = 1;
                bound = zero;
                goto S290;
                S280:
                data.status = 2;
                bound = inf;
                S300:
                S290: ;
                break;
            }
            case 5:
            {
                //
                //     Calculating PHONC
                //
                phonc = 5.0e0;
                T15 = tent4;
                T16 = atol;
                T17 = tol;
                E0000E0001.dstinv(ref data, K1, T15, K3, K3, K4, T16, T17);
                data.status = 0;
                data.x = phonc;
                E0000E0001.dinvr(ref data);
                S310:
                switch ((data.status == 1))
                {
                    case false:
                        goto S320;
                }

                cumfnc(f, dfn, dfd, phonc, ref cum, ref ccum);
                data.fx = cum - p;
                data.x = phonc;
                E0000E0001.dinvr(ref data);
                goto S310;
                S320:
                switch ((data.status == -1))
                {
                    case false:
                        goto S350;
                }

                switch (data.qleft)
                {
                    case false:
                        goto S330;
                }

                data.status = 1;
                bound = 0.0e0;
                goto S340;
                S330:
                data.status = 2;
                bound = tent4;
                S340: ;
                break;
            }
        }

        S350:
        status_ = data.status;
    }
}