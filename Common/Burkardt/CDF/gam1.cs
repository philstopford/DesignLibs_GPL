﻿namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static double gam1 ( double a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5D+00 <= A <= 1.5
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
        //    Input, double *A, forms the argument of the Gamma function.
        //
        //    Output, double GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
        //
    {
        const double s1 = .273076135303957e+00;
        const double s2 = .559398236957378e-01;
        double[] p = {
            .577215664901533e+00,-.409078193005776e+00,-.230975380857675e+00,
            .597275330452234e-01,.766968181649490e-02,-.514889771323592e-02,
            .589597428611429e-03
        };
        double[] q = {
            .100000000000000e+01,.427569613095214e+00,.158451672430138e+00,
            .261132021441447e-01,.423244297896961e-02
        };
        double[] r = {
            -.422784335098468e+00,-.771330383816272e+00,-.244757765222226e+00,
            .118378989872749e+00,.930357293360349e-03,-.118290993445146e-01,
            .223047661158249e-02,.266505979058923e-03,-.132674909766242e-03
        };
        double gam1,bot, top,w;

        double t = a;
        double d = a-0.5e0;
        t = d switch
        {
            > 0.0e0 => d - 0.5e0,
            _ => t
        };

        double T1 = t;
        switch (T1)
        {
            case < 0:
                goto S40;
            case 0:
                goto S10;
            default:
                goto S20;
        }
        S10:
        gam1 = 0.0e0;
        return gam1;
        S20:
        top = (((((p[6]*t+p[5])*t+p[4])*t+p[3])*t+p[2])*t+p[1])*t+p[0];
        bot = (((q[4]*t+q[3])*t+q[2])*t+q[1])*t+1.0e0;
        w = top/bot;
        switch (d)
        {
            case > 0.0e0:
                goto S30;
        }

        gam1 = a*w;
        return gam1;
        S30:
        gam1 = t/ a*(w-0.5e0-0.5e0);
        return gam1;
        S40:
        top = (((((((r[8]*t+r[7])*t+r[6])*t+r[5])*t+r[4])*t+r[3])*t+r[2])*t+r[1])*t+
              r[0];
        bot = (s2*t+s1)*t+1.0e0;
        w = top/bot;
        switch (d)
        {
            case > 0.0e0:
                goto S50;
        }

        gam1 = a*(w+0.5e0+0.5e0);
        return gam1;
        S50:
        gam1 = t*w/ a;
        return gam1;
    }
}