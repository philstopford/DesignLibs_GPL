using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double rlog(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RLOG computes  X - 1 - LN(X).
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
        //  Parameters:
        //
        //    Input, double *X, the argument of the function.
        //
        //    Output, double RLOG, the value of the function.
        //
    {
        const double a = .566749439387324e-01;
        const double b = .456512608815524e-01;
        const double p0 = .333333333333333e+00;
        const double p1 = -.224696413112536e+00;
        const double p2 = .620886815375787e-02;
        const double q1 = -.127408923933623e+01;
        const double q2 = .354508718369557e+00;
        double rlog, r, u, w1;

        switch (x)
        {
            case < 0.61e0:
            case > 1.57e0:
                goto S40;
        }

        switch (x)
        {
            case < 0.82e0:
                goto S10;
            case > 1.18e0:
                goto S20;
        }

        //
        //  ARGUMENT REDUCTION
        //
        u = x - 0.5e0 - 0.5e0;
        w1 = 0.0e0;
        goto S30;
        S10:
        u = x - 0.7e0;
        u /= 0.7e0;
        w1 = a - u * 0.3e0;
        goto S30;
        S20:
        u = 0.75e0 * x - 1e0;
        w1 = b + u / 3.0e0;
        S30:
        //
        //  SERIES EXPANSION
        //
        r = u / (u + 2.0e0);
        double t = r * r;
        double w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0e0);
        rlog = 2.0e0 * t * (1.0e0 / (1.0e0 - r) - r * w) + w1;
        return rlog;
        S40:
        r = x - 0.5e0 - 0.5e0;
        rlog = r - Math.Log(x);
        return rlog;
    }

    public static double rlog1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RLOG1 evaluates the function X - ln ( 1 + X ).
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
        //  Parameters:
        //
        //    Input, double *X, the argument.
        //
        //    Output, double RLOG1, the value of X - ln ( 1 + X ).
        //
    {
        const double a = .566749439387324e-01;
        const double b = .456512608815524e-01;
        double h;
        const double p0 = .333333333333333e+00;
        const double p1 = -.224696413112536e+00;
        const double p2 = .620886815375787e-02;
        const double q1 = -.127408923933623e+01;
        const double q2 = .354508718369557e+00;
        double rlog1, w1;

        switch (x)
        {
            case < -0.39e0:
            case > 0.57e0:
                goto S40;
        }

        switch (x)
        {
            case < -0.18e0:
                goto S10;
            case > 0.18e0:
                goto S20;
        }

        //
        //  ARGUMENT REDUCTION
        //
        h = x;
        w1 = 0.0e0;
        goto S30;
        S10:
        h = x + 0.3e0;
        h /= 0.7e0;
        w1 = a - h * 0.3e0;
        goto S30;
        S20:
        h = 0.75e0 * x - 0.25e0;
        w1 = b + h / 3.0e0;
        S30:
        //
        //  SERIES EXPANSION
        //
        double r = h / (h + 2.0e0);
        double t = r * r;
        double w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.0e0);
        rlog1 = 2.0e0 * t * (1.0e0 / (1.0e0 - r) - r * w) + w1;
        return rlog1;
        S40:
        w = x + 0.5e0 + 0.5e0;
        rlog1 = x - Math.Log(w);
        return rlog1;
    }
}