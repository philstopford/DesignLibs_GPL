using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfnor(int which, ref double p, ref double q, ref double x, ref double mean,
                ref double sd, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFNOR evaluates the CDF of the Normal distribution.
            //
            //  Discussion:
            //
            //    A slightly modified version of ANORM from SPECFUN
            //    is used to calculate the cumulative standard normal distribution.
            //
            //    The rational functions from pages 90-95 of Kennedy and Gentle
            //    are used as starting values to Newton's Iterations which
            //    compute the inverse standard normal.  Therefore no searches are
            //    necessary for any parameter.
            //
            //    For X < -15, the asymptotic expansion for the normal is used  as
            //    the starting value in finding the inverse standard normal.
            //
            //    The normal density is proportional to
            //    exp( - 0.5 * (( X - MEAN)/SD)^2)
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
            //    1966, Formula 26.2.12.
            //
            //    William Cody,
            //    Algorithm 715: SPECFUN - A Portable FORTRAN Package of
            //      Special Function Routines and Test Drivers,
            //    ACM Transactions on Mathematical Software,
            //    Volume 19, pages 22-32, 1993.
            //
            //    Kennedy and Gentle,
            //    Statistical Computing,
            //    Marcel Dekker, NY, 1980,
            //    QA276.4  K46
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, MEAN and SD;
            //    2: Calculate X from P, Q, MEAN and SD;
            //    3: Calculate MEAN from P, Q, X and SD;
            //    4: Calculate SD from P, Q, X and MEAN.
            //
            //    Input/output, double *P, the integral from -infinity to X
            //    of the Normal density.  If this is an input or output value, it will
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration of
            //    the Normal density.
            //
            //    Input/output, double *MEAN, the mean of the Normal density.
            //
            //    Input/output, double *SD, the standard deviation of the
            //    Normal density.  If this is an input value, it should lie in the
            //    range (0,+infinity).
            //
            //    Output, int *STATUS, the status of the calculation.
            //    0, if calculation completed correctly;
            //    -I, if input parameter number I is out of range;
            //    1, if answer appears to be lower than lowest search bound;
            //    2, if answer appears to be higher than greatest search bound;
            //    3, if P + Q /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            int K1 = 1;
            double pq;
            double z;

            status = 0;

            status = 0;
            bound = 0.0;
            //
            //  Check arguments
            //
            status = 0;
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
            //     P
            //
            if (!(p <= 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p <= 0.0e0)) goto S40;
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
            if (which == 1) goto S150;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S140;
            if (!(pq < 0.0e0)) goto S120;
            bound = 0.0e0;
            goto S130;
            S120:
            bound = 1.0e0;
            S130:
            status = 3;
            return;
            S150:
            S140:
            if (which == 4) goto S170;
            //
            //     SD
            //
            if (!(sd <= 0.0e0)) goto S160;
            bound = 0.0e0;
            status = -6;
            return;
            S170:
            S160:
            //
            //  Computing P
            //
            if (1 == which)
            {
                z = (x - mean) / sd;
                cumnor(z, ref p, ref q);
            }
            //
            //  Computing X
            //
            else if (2 == which)
            {
                z = dinvnr(p, q);
                x = sd * z + mean;
            }
            //
            //  Computing the MEAN
            //
            else if (3 == which)
            {
                z = dinvnr(p, q);
                mean = x - sd * z;
            }
            //
            //  Computing SD
            //
            else if (4 == which)
            {
                z = dinvnr(p, q);
                sd = (x - mean) / z;
            }
        }
    }
}