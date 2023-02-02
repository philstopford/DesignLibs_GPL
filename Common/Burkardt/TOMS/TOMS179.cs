using System;
using Burkardt.AppliedStatistics;
using Burkardt.Types;

namespace Burkardt.TOMSNS;

public static partial class TOMS
{
    public static double mdbeta(double x, double p, double q, ref int ier)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MDBETA evaluates the incomplete beta function.
        //
        //  Modified:
        //
        //    30 January 2008
        //
        //  Author:
        //
        //    Oliver Ludwig
        //    Modifications by John Burkardt
        //
        //  Reference:
        //
        //    Oliver Ludwig,
        //    Algorithm 179:
        //    Incomplete Beta Ratio,
        //    Communications of the ACM,
        //    Volume 6, Number 6, June 1963, page 314.
        //
        //  Parameters:
        //
        //    Input, double X, the value to which function is to be
        //    integrated.  X must be in the range [0,1] inclusive.
        //
        //    Input, double P, the first parameter.  P must be greater
        //    than 0.0.
        //
        //    Input, double Q, the second parameter.  Q must be greater
        //    than 0.0.
        //
        //    Output, int *IER, error parameter.
        //    0, normal exit.
        //    1, X is not in the range [0,1] inclusive.
        //    2, P or Q is less than or equal to 0.
        //
        //    Output, double MDBETA.  The probability that a random variable
        //    from a Beta distribution having parameters P and Q will be less than
        //    or equal to X.
        //
        //  Local parameters:
        //
        //    Local, double ALEPS, the logarithm of EPS1.
        //
        //    Local, double EPS, the machine precision.
        //
        //    Local, double EPS1, the smallest representable number.
        //
    {
        const double aleps = -179.6016;
        double cnt;
        const double eps = 2.2E-16;
        const double eps1 = 1.0E-78;
        int ifault = 0;
        int interval;
        double temp;
        double wh;
        //
        //  Check ranges of the arguments.
        //
        double prob = 0.0;
        double y = x;

        switch (x)
        {
            case < 0.0:
            case > 1.0:
                ier = 1;
                return prob;
        }

        if (p <= 0.0 || q <= 0.0)
        {
            ier = 2;
            return prob;
        }

        ier = 0;

        switch (x)
        {
            case <= 0.5:
                interval = 0;
                break;
            default:
                interval = 1;
                temp = p;
                p = q;
                q = temp;
                y = 1.0 - y;
                break;
        }

        switch (x)
        {
            case 0.0:
            case 1.0:
            {
                prob = 0.0;

                if (interval == 0)
                {
                    return prob;
                }

                prob = 1.0 - prob;

                return prob;
            }
        }

        int ib = (int)q;
        temp = ib;
        double ps = q - ib;

        if (Math.Abs(q - temp) <= typeMethods.r8_epsilon())
        {
            ps = 1.0;
        }

        double dp = p;
        double dq = q;
        double px = dp * Math.Log(y);
        double pq = Algorithms.alogam(dp + dq, ref ifault);
        double p1 = Algorithms.alogam(dp, ref ifault);
        double c = Algorithms.alogam(dq, ref ifault);
        double d4 = Math.Log(dp);
        double xb = px + Algorithms.alogam(ps + dp, ref ifault) - Algorithms.alogam(ps, ref ifault) - d4 - p1;
        //
        //  Scaling
        //
        ib = (int) (xb / aleps);
        double infsum = 0.0;
        switch (ib)
        {
            //
            //  First term of a decreasing series will underflow.
            //
            case 0:
            {
                infsum = Math.Exp(xb);
                cnt = infsum * dp;
                //
                //  CNT will equal dexp ( temp ) * ( 1.d0 - ps ) * i * p * y**i / factorial ( i ).
                //
                wh = 0.0;

                for (;;)
                {
                    wh += 1.0;
                    cnt = cnt * (wh - ps) * y / wh;
                    xb = cnt / (dp + wh);
                    infsum += xb;

                    if (xb / eps < infsum)
                    {
                        break;
                    }
                }

                break;
            }
        }

        double finsum = 0.0;

        switch (dq)
        {
            case <= 1.0:
            {
                prob = finsum + infsum;

                if (interval == 0)
                {
                    return prob;
                }

                prob = 1.0 - prob;

                return prob;
            }
        }

        xb = px + dq * Math.Log(1.0 - y) + pq - p1 - Math.Log(dq) - c;

        ib = ib switch
        {
            < 0 => 0,
            //
            //  Scaling.
            //
            _ => (int) (xb / aleps)
        };

        c = 1.0 / (1.0 - y);
        cnt = Math.Exp(xb - ib * aleps);
        ps = dq;
        wh = dq;

        for (;;)
        {
            wh -= 1.0;

            if (wh <= 0.0)
            {
                prob = finsum + infsum;

                if (interval != 0)
                {
                    prob = 1.0 - prob;
                }

                break;
            }

            px = ps * c / (dp + wh);

            if (px <= 1.0)
            {
                if (cnt / eps <= finsum || cnt <= eps1 / px)
                {
                    prob = finsum + infsum;

                    if (interval != 0)
                    {
                        prob = 1.0 - prob;
                    }

                    break;
                }
            }

            cnt *= px;
            switch (cnt)
            {
                //
                //  Rescale.
                //
                case > 1.0:
                    ib -= 1;
                    cnt *= eps1;
                    break;
            }

            ps = wh;

            switch (ib)
            {
                case 0:
                    finsum += cnt;
                    break;
            }
        }

        return prob;
    }
}