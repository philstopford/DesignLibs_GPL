using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfgam(int which, ref double p, ref double q, ref double x_, ref double shape,
                ref double scale, ref int status_, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFGAM evaluates the CDF of the Gamma Distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the Gamma distribution
            //    given the others.
            //
            //    The cumulative distribution function P is calculated directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The gamma density is proportional to T^(SHAPE - 1) * EXP(- SCALE * T)
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
            //    Computation of the incomplete gamma function ratios and their inverse,
            //    ACM Transactions on Mathematical Software,
            //    Volume 12, 1986, pages 377-393.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, SHAPE and SCALE;
            //    2: Calculate X from P, Q, SHAPE and SCALE;
            //    3: Calculate SHAPE from P, Q, X and SCALE;
            //    4: Calculate SCALE from P, Q, X and SHAPE.
            //
            //    Input/output, double *P, the integral from 0 to X of the
            //    Gamma density.  If this is an input value, it should lie in the
            //    range: [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration of
            //    the Gamma density.  If this is an input value, it should lie in the
            //    range: [0, +infinity).  If it is an output value, it will lie in
            //    the range: [0,1E300].
            //
            //    Input/output, double *SHAPE, the shape parameter of the
            //    Gamma density.  If this is an input value, it should lie in the range:
            //    (0, +infinity).  If it is an output value, it will be searched for
            //    in the range: [1.0D-300,1.0D+300].
            //
            //    Input/output, double *SCALE, the scale parameter of the
            //    Gamma density.  If this is an input value, it should lie in the range
            //    (0, +infinity).  If it is an output value, it will be searched for
            //    in the range: (1.0D-300,1.0D+300].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +10, if the Gamma or inverse Gamma routine cannot compute the answer.
            //    This usually happens only for X and SHAPE very large (more than 1.0D+10.
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
            int ierr = 0;
            int K1 = 1;
            double K5 = 0.5e0;
            double K6 = 5.0e0;
            double porq = 0;
            double pq = 0;
            bool qporq = false;
            double T2 = 0;
            double T3 = 0;
            double T4 = 0;
            double T7 = 0;
            double T8 = 0;
            double T9 = 0;
            double xscale = 0;
            double xx = 0;

            E0000_E0001_Data data = new E0000_E0001_Data();

            data.status = 0;
            data.x = x_;
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
            //     P
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
            //  SHAPE
            //
            if (!(shape <= 0.0e0)) goto S140;
            bound = 0.0e0;
            data.status = -5;
            return;
            S150:
            S140:
            if (which == 4) goto S170;
            //
            //  SCALE
            //
            if (!(scale <= 0.0e0)) goto S160;
            bound = 0.0e0;
            data.status = -6;
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
            data.status = 3;
            return;
            S210:
            S200:
            if (which == 1) goto S240;
            //
            //     Select the minimum of P or Q
            //
            qporq = p <= q;
            if (!qporq) goto S220;
            porq = p;
            goto S230;
            S220:
            porq = q;
            S240:
            S230:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //     Calculating P
                //
                data.status = 0;
                xscale = data.x * scale;
                cumgam(xscale, shape, ref p, ref q);
                if (porq > 1.5e0) data.status = 10;
            }
            else if (2 == which)
            {
                //
                //     Computing X
                //
                T2 = -1.0e0;
                gamma_inc_inv(shape, xx, T2, p, q, ierr);
                if (ierr < 0.0e0)
                {
                    data.status = 10;
                    return;
                }
                else
                {
                    data.x = xx / scale;
                    data.status = 0;
                }
            }
            else if (3 == which)
            {
                //
                //     Computing SHAPE
                //
                shape = 5.0e0;
                xscale = data.x * scale;
                T3 = zero;
                T4 = inf;
                T7 = atol;
                T8 = tol;
                E0000E0001.dstinv(ref data, T3, T4, K5, K5, K6, T7, T8);
                data.status = 0;
                data.x = shape;
                E0000E0001.dinvr(ref data);
                S250:
                if (!(data.status == 1)) goto S290;
                cumgam(xscale, shape, ref cum, ref ccum);
                if (!qporq) goto S260;
                data.fx = cum - p;
                goto S270;
                S260:
                data.fx = ccum - q;
                S270:
                if (!((qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))) goto S280;
                data.status = 10;
                return;
                S280:
                data.x = shape;
                E0000E0001.dinvr(ref data);
                goto S250;
                S290:
                if (!(data.status == -1)) goto S320;
                if (!data.qleft) goto S300;
                data.status = 1;
                bound = zero;
                goto S310;
                S300:
                data.status = 2;
                bound = inf;
                S320:
                S310: ;
            }
            else if (4 == which)
            {
                //
                //  Computing SCALE
                //
                T9 = -1.0e0;
                gamma_inc_inv(shape, xx, T9, p, q, ierr);
                if (ierr < 0.0e0)
                {
                    data.status = 10;
                    return;
                }
                else
                {
                    scale = xx / data.x;
                    data.status = 0;
                }
            }

            status_ = data.status;
            x_ = data.x;
            return;
        }
    }
}