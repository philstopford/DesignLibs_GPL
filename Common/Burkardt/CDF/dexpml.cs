using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double dexpm1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEXPM1 evaluates the function EXP(X) - 1.
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
        //    Input, double *X, the value at which exp(X)-1 is desired.
        //
        //    Output, double DEXPM1, the value of exp(X)-1.
        //
    {
        double p1 = 0.914041914819518e-09;
        double p2 = 0.238082361044469e-01;
        double q1 = -0.499999999085958e+00;
        double q2 = 0.107141568980644e+00;
        double q3 = -0.119041179760821e-01;
        double q4 = 0.595130811860248e-03;
        double dexpm1;
        double w;

        switch (Math.Abs(x))
        {
            case <= 0.15e0:
                dexpm1 = x * (((
                                   p2 * x
                                   + p1) * x
                               + 1.0e0)
                              / ((((
                                       q4 * x
                                       + q3) * x
                                   + q2) * x
                                  + q1) * x
                                 + 1.0e0));
                break;
            default:
            {
                switch (x)
                {
                    case <= 0.0e0:
                        w = Math.Exp(x);
                        dexpm1 = w - 0.5e0 - 0.5e0;
                        break;
                    default:
                        w = Math.Exp(x);
                        dexpm1 = w * (0.5e0 + (0.5e0 - 1.0e0 / w));
                        break;
                }

                break;
            }
        }

        return dexpm1;
    }
}