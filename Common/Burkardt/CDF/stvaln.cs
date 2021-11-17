using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double stvaln(double p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STVALN provides starting values for the inverse of the normal distribution.
        //
        //  Discussion:
        //
        //    The routine returns X such that
        //      P = CUMNOR(X),
        //    that is,
        //      P = Integral from -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU.
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
        //    Kennedy and Gentle,
        //    Statistical Computing,
        //    Marcel Dekker, NY, 1980, page 95,
        //    QA276.4  K46
        //
        //  Parameters:
        //
        //    Input, double *P, the probability whose normal deviate
        //    is sought.
        //
        //    Output, double STVALN, the normal deviate whose probability
        //    is P.
        //
    {
        double[] xden =  {
                0.993484626060e-1,0.588581570495e0,0.531103462366e0,0.103537752850e0,
                0.38560700634e-2
            }
            ;
        double[] xnum =  {
                -0.322232431088e0,-1.000000000000e0,-0.342242088547e0,-0.204231210245e-1,
                -0.453642210148e-4
            }
            ;
        int K1 = 5;
        double sign;
        double stvaln, y, z;

        switch ((p <= 0.5e0))
        {
            case false:
                goto S10;
        }
        sign = -1.0e0;
        z = p;
        goto S20;
        S10:
        sign = 1.0e0;
        z = 1.0e0 - p;
        S20:
        y = Math.Sqrt(-(2.0e0 * Math.Log(z)));
        stvaln = y + eval_pol(xnum, K1, y) / eval_pol(xden, K1, y);
        stvaln = sign * stvaln;
        return stvaln;
    }
}