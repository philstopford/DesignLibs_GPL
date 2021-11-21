using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double error_fc(int ind, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERROR_FC evaluates the complementary error function ERFC.
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
        //    Input, int *IND, chooses the scaling.
        //    If IND is nonzero, then the value returned has been multiplied by
        //    EXP(X*X).
        //
        //    Input, double *X, the argument of the function.
        //
        //    Output, double ERROR_FC, the value of the complementary
        //    error function.
        //
    {
        double c = .564189583547756e0;
        double[] a =  {
                .771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
                .479137145607681e-01,.128379167095513e+00
            }
            ;
        double[] b =  {
                .301048631703895e-02,.538971687740286e-01,.375795757275549e+00
            }
            ;
        double[] p =  {
                -1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
                4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
                4.51918953711873e+02,3.00459261020162e+02
            }
            ;
        double[] q =  {
                1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
                2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
                7.90950925327898e+02,3.00459260956983e+02
            }
            ;
        double[] r =  {
                2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
                4.65807828718470e+00,2.82094791773523e-01
            }
            ;
        double[] s =  {
                9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
                1.80124575948747e+01
            }
            ;
        const int K1 = 1;
        double erfc1, bot, t, top;

        //
        //  ABS(X) <= 0.5
        //
        double ax = Math.Abs(x);
        switch (ax)
        {
            case > 0.5e0:
                goto S10;
        }
        t = x * x;
        top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.0e0;
        bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
        erfc1 = 0.5e0 + (0.5e0 - x * (top / bot));
        if (ind != 0)
        {
            erfc1 = Math.Exp(t) * erfc1;
        }

        return erfc1;
        S10:
        switch (ax)
        {
            //
            //  0.5 < ABS(X) <= 4
            //
            case > 4.0e0:
                goto S20;
        }
        top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax + p[5]) * ax + p[6]) * ax + p[
            7];
        bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax + q[5]) * ax + q[6]) * ax + q[
            7];
        erfc1 = top / bot;
        goto S40;
        S20:
        switch (x)
        {
            //
            //  4 < ABS(X)
            //
            case <= -5.6e0:
                goto S60;
        }
        if (ind != 0)
        {
            goto S30;
        }

        switch (x)
        {
            case > 100.0e0:
                goto S70;
        }
        if (x * x > -exparg(K1))
        {
            goto S70;
        }

        S30:
        t = Math.Pow(1.0e0 / x, 2.0);
        top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
        bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
        erfc1 = (c - t * top / bot) / ax;
        S40:
        switch (ind)
        {
            //
            //  FINAL ASSEMBLY
            //
            case 0:
                goto S50;
        }

        erfc1 = x switch
        {
            < 0.0e0 => 2.0e0 * Math.Exp(x * x) - erfc1,
            _ => erfc1
        };
        return erfc1;
        S50:
        double w = x * x;
        t = w;
        double e = w - t;
        erfc1 = x switch
        {
            < 0.0e0 => 2.0e0 - erfc1,
            _ => (0.5e0 + (0.5e0 - e)) * Math.Exp(-t) * erfc1
        };
        return erfc1;
        S60:
        //
        //  LIMIT VALUE FOR LARGE NEGATIVE X
        //
        erfc1 = 2.0e0;
        if (ind != 0)
        {
            erfc1 = 2.0e0 * Math.Exp(x * x);
        }

        return erfc1;
        S70:
        //
        //  LIMIT VALUE FOR LARGE POSITIVE X WHEN IND = 0
        //
        erfc1 = 0.0e0;
        return erfc1;
    }

}