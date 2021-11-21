using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void gamma_inc(double a, double x, ref double ans, ref double qans, int ind)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
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
        //    Alfred H Morris, Jr
        //
        //  Parameters:
        //
        //    Input, double *A, *X, the arguments of the incomplete
        //    gamma ratio.  A and X must be nonnegative.  A and X cannot
        //    both be zero.
        //
        //    Output, double *ANS, *QANS.  On normal output,
        //    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
        //    A or X is negative, or both are 0, or when the answer is
        //    computationally indeterminate because A is extremely large
        //    and X is very close to A.
        //
        //    Input, int *IND, indicates the accuracy request:
        //    0, as much accuracy as possible.
        //    1, to within 1 unit of the 6-th significant digit,
        //    otherwise, to within 1 unit of the 3rd significant digit.
        //
    {
        double[] acc0 =  {
                5e-15,5e-7,5e-4
            }
            ;
        const double alog10 = 2.30258509299405e0;
        const double d10 = -.185185185185185e-02;
        const double d20 = .413359788359788e-02;
        const double d30 = .649434156378601e-03;
        const double d40 = -.861888290916712e-03;
        const double d50 = -.336798553366358e-03;
        const double d60 = .531307936463992e-03;
        const double d70 = .344367606892378e-03;
        const double rt2pin = .398942280401433e0;
        const double rtpi = 1.77245385090552e0;
        const double third = .333333333333333e0;


        double[] big =  {
                20.0e0,14.0e0,10.0e0
            }
            ;
        double[] d0 =  {
                .833333333333333e-01,-.148148148148148e-01,.115740740740741e-02,
                .352733686067019e-03,-.178755144032922e-03,.391926317852244e-04,
                -.218544851067999e-05,-.185406221071516e-05,.829671134095309e-06,
                -.176659527368261e-06,.670785354340150e-08,.102618097842403e-07,
                -.438203601845335e-08
            }
            ;
        double[] d1 =  {
                -.347222222222222e-02,.264550264550265e-02,-.990226337448560e-03,
                .205761316872428e-03,-.401877572016461e-06,-.180985503344900e-04,
                .764916091608111e-05,-.161209008945634e-05,.464712780280743e-08,
                .137863344691572e-06,-.575254560351770e-07,.119516285997781e-07
            }
            ;
        double[] d2 =  {
                -.268132716049383e-02,.771604938271605e-03,.200938786008230e-05,
                -.107366532263652e-03,.529234488291201e-04,-.127606351886187e-04,
                .342357873409614e-07,.137219573090629e-05,-.629899213838006e-06,
                .142806142060642e-06
            }
            ;
        double[] d3 =  {
                .229472093621399e-03,-.469189494395256e-03,.267720632062839e-03,
                -.756180167188398e-04,-.239650511386730e-06,.110826541153473e-04,
                -.567495282699160e-05,.142309007324359e-05
            }
            ;
        double[] d4 =  {
                .784039221720067e-03,-.299072480303190e-03,-.146384525788434e-05,
                .664149821546512e-04,-.396836504717943e-04,.113757269706784e-04
            }
            ;
        double[] d5 =  {
                -.697281375836586e-04,.277275324495939e-03,-.199325705161888e-03,
                .679778047793721e-04
            }
            ;
        double[] d6 =  {
                -.592166437353694e-03,.270878209671804e-03
            }
            ;
        double[] e00 =  {
                .25e-3,.25e-1,.14e0
            }
            ;
        double[] x00 =  {
                31.0e0,17.0e0,9.7e0
            }
            ;
        const int K1 = 1;
        const int K2 = 0;
        double a2n,
            c,
            c0,
            c1,
            c2,
            c3,
            c4,
            c5,
            c6,
            l,
            r,
            sum,
            t,
            t1,
            tol,
            u,
            w,
            z;
        int m, max, n;
        double[] wk = new double[20];

        //
        //  E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST
        //  NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
        //
        double e = dpmpar(K1);
        if (a < 0.0e0 || x < 0.0e0)
        {
            goto S430;
        }

        switch (a)
        {
            case 0.0e0 when x == 0.0e0:
                goto S430;
        }

        switch (a * x)
        {
            case 0.0e0:
                goto S420;
        }

        int iop = ind + 1;
        if (iop != 1 && iop != 2)
        {
            iop = 3;
        }

        double acc = Math.Max(acc0[iop - 1], e);
        double e0 = e00[iop - 1];
        double x0 = x00[iop - 1];
        switch (a)
        {
            //
            //  SELECT THE APPROPRIATE ALGORITHM
            //
            case >= 1.0e0:
                goto S10;
            case 0.5e0:
                goto S390;
        }

        switch (x)
        {
            case < 1.1e0:
                goto S160;
        }

        t1 = a * Math.Log(x) - x;
        u = a * Math.Exp(t1);
        switch (u)
        {
            case 0.0e0:
                goto S380;
        }

        r = u * (1.0e0 + gam1(a));
        goto S250;
        S10:
        if (a >= big[iop - 1])
        {
            goto S30;
        }

        if (a > x || x >= x0)
        {
            goto S20;
        }

        double twoa = a + a;
        m = (int)Math.Truncate(twoa);
        if (Math.Abs(twoa - m) > double.Epsilon)
        {
            goto S20;
        }

        int i = m / 2;
        if (Math.Abs(a - i) <= double.Epsilon)
        {
            goto S210;
        }

        goto S220;
        S20:
        t1 = a * Math.Log(x) - x;
        r = Math.Exp(t1) / gamma_x(a);
        goto S40;
        S30:
        l = x / a;
        switch (l)
        {
            case 0.0e0:
                goto S370;
        }

        double s = 0.5e0 + (0.5e0 - l);
        z = rlog(l);
        if (z >= 700.0e0 / a)
        {
            goto S410;
        }

        double y = a * z;
        double rta = Math.Sqrt(a);
        if (Math.Abs(s) <= e0 / rta)
        {
            goto S330;
        }

        switch (Math.Abs(s))
        {
            case <= 0.4e0:
                goto S270;
        }

        t = Math.Pow(1.0e0 / a, 2.0);
        t1 = (((0.75e0 * t - 1.0e0) * t + 3.5e0) * t - 105.0e0) / (a * 1260.0e0);
        t1 -= y;
        r = rt2pin * rta * Math.Exp(t1);
        S40:
        switch (r)
        {
            case 0.0e0:
                goto S420;
        }

        if (x <= Math.Max(a, alog10))
        {
            goto S50;
        }

        if (x < x0)
        {
            goto S250;
        }

        goto S100;
        S50:
        //
        //  TAYLOR SERIES FOR P/R
        //
        double apn = a + 1.0e0;
        t = x / apn;
        wk[0] = t;
        for (n = 2; n <= 20; n++)
        {
            apn += 1.0e0;
            t *= x / apn;
            switch (t)
            {
                case <= 1e-3:
                    goto S70;
            }

            wk[n - 1] = t;
        }

        n = 20;
        S70:
        sum = t;
        tol = 0.5e0 * acc;
        S80:
        apn += 1.0e0;
        t *= x / apn;
        sum += t;
        if (t > tol)
        {
            goto S80;
        }

        max = n - 1;
        for (m = 1; m <= max; m++)
        {
            n -= 1;
            sum += wk[n - 1];
        }

        ans = r / a * (1.0e0 + sum);
        qans = 0.5e0 + (0.5e0 - ans);
        return;
        S100:
        //
        //  ASYMPTOTIC EXPANSION
        //
        double amn = a - 1.0e0;
        t = amn / x;
        wk[0] = t;
        for (n = 2; n <= 20; n++)
        {
            amn -= 1.0e0;
            t *= amn / x;
            switch (Math.Abs(t))
            {
                case <= 1e-3:
                    goto S120;
            }

            wk[n - 1] = t;
        }

        n = 20;
        S120:
        sum = t;
        S130:
        if (Math.Abs(t) <= acc)
        {
            goto S140;
        }

        amn -= 1.0e0;
        t *= amn / x;
        sum += t;
        goto S130;
        S140:
        max = n - 1;
        for (m = 1; m <= max; m++)
        {
            n -= 1;
            sum += wk[n - 1];
        }

        qans = r / x * (1.0e0 + sum);
        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S160:
        //
        //  TAYLOR SERIES FOR P(A,X)/X**A
        //
        double an = 3.0e0;
        c = x;
        sum = x / (a + 3.0e0);
        tol = 3.0e0 * acc / (a + 1.0e0);
        S170:
        an += 1.0e0;
        c = -(c * (x / an));
        t = c / (a + an);
        sum += t;
        if (Math.Abs(t) > tol)
        {
            goto S170;
        }

        double j = a * x * ((sum / 6.0e0 - 0.5e0 / (a + 2.0e0)) * x + 1.0e0 / (a + 1.0e0));
        z = a * Math.Log(x);
        double h = gam1(a);
        double g = 1.0e0 + h;
        switch (x)
        {
            case < 0.25e0:
                goto S180;
        }

        if (a < x / 2.59e0)
        {
            goto S200;
        }

        goto S190;
        S180:
        switch (z)
        {
            case > -.13394e0:
                goto S200;
        }

        S190:
        w = Math.Exp(z);
        ans = w * g * (0.5e0 + (0.5e0 - j));
        qans = 0.5e0 + (0.5e0 - ans);
        return;
        S200:
        l = rexp(z);
        w = 0.5e0 + (0.5e0 + l);
        qans = (w * j - l) * g - h;
        switch (qans)
        {
            case < 0.0e0:
                goto S380;
        }

        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S210:
        //
        //  FINITE SUMS FOR Q WHEN A .GE. 1 AND 2*A IS AN INTEGER
        //
        sum = Math.Exp(-x);
        t = sum;
        n = 1;
        c = 0.0e0;
        goto S230;
        S220:
        double rtx = Math.Sqrt(x);
        sum = error_fc(K2, rtx);
        t = Math.Exp(-x) / (rtpi * rtx);
        n = 0;
        c = -0.5e0;
        S230:
        if (n == i)
        {
            goto S240;
        }

        n += 1;
        c += 1.0e0;
        t = x * t / c;
        sum += t;
        goto S230;
        S240:
        qans = sum;
        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S250:
        //
        //  CONTINUED FRACTION EXPANSION
        //
        tol = Math.Max(5.0e0 * e, acc);
        double a2nm1 = a2n = 1.0e0;
        double b2nm1 = x;
        double b2n = x + (1.0e0 - a);
        c = 1.0e0;
        S260:
        a2nm1 = x * a2n + c * a2nm1;
        b2nm1 = x * b2n + c * b2nm1;
        double am0 = a2nm1 / b2nm1;
        c += 1.0e0;
        double cma = c - a;
        a2n = a2nm1 + cma * a2n;
        b2n = b2nm1 + cma * b2n;
        double an0 = a2n / b2n;
        if (Math.Abs(an0 - am0) >= tol * an0)
        {
            goto S260;
        }

        qans = r * an0;
        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S270:
        //
        //  GENERAL TEMME EXPANSION
        //
        if (Math.Abs(s) <= 2.0e0 * e && a * e * e > 3.28e-3)
        {
            goto S430;
        }

        c = Math.Exp(-y);
        double T3 = Math.Sqrt(y);
        w = 0.5e0 * error_fc(K1, T3);
        u = 1.0e0 / a;
        z = l switch
        {
            < 1.0e0 => -z,
            _ => Math.Sqrt(z + z)
        };

        int T4 = iop - 2;
        switch (T4)
        {
            case < 0:
                goto S280;
            case 0:
                goto S290;
            default:
                goto S300;
        }
        S280:
        switch (Math.Abs(s))
        {
            case <= 1e-3:
                goto S340;
        }

        c0 = ((((((((((((d0[12] * z + d0[11]) * z + d0[10]) * z + d0[9]) * z + d0[8]) * z + d0[7]) * z + d0[
            6]) * z + d0[5]) * z + d0[4]) * z + d0[3]) * z + d0[2]) * z + d0[1]) * z + d0[0]) * z - third;
        c1 = (((((((((((d1[11] * z + d1[10]) * z + d1[9]) * z + d1[8]) * z + d1[7]) * z + d1[6]) * z + d1[5]
            ) * z + d1[4]) * z + d1[3]) * z + d1[2]) * z + d1[1]) * z + d1[0]) * z + d10;
        c2 = (((((((((d2[9] * z + d2[8]) * z + d2[7]) * z + d2[6]) * z + d2[5]) * z + d2[4]) * z + d2[3]) * z +
                d2[2]) * z + d2[1]) * z + d2[0]) * z + d20;
        c3 = (((((((d3[7] * z + d3[6]) * z + d3[5]) * z + d3[4]) * z + d3[3]) * z + d3[2]) * z + d3[1]) * z +
              d3[0]) * z + d30;
        c4 = (((((d4[5] * z + d4[4]) * z + d4[3]) * z + d4[2]) * z + d4[1]) * z + d4[0]) * z + d40;
        c5 = (((d5[3] * z + d5[2]) * z + d5[1]) * z + d5[0]) * z + d50;
        c6 = (d6[1] * z + d6[0]) * z + d60;
        t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0;
        goto S310;
        S290:
        c0 = (((((d0[5] * z + d0[4]) * z + d0[3]) * z + d0[2]) * z + d0[1]) * z + d0[0]) * z - third;
        c1 = (((d1[3] * z + d1[2]) * z + d1[1]) * z + d1[0]) * z + d10;
        c2 = d2[0] * z + d20;
        t = (c2 * u + c1) * u + c0;
        goto S310;
        S300:
        t = ((d0[2] * z + d0[1]) * z + d0[0]) * z - third;
        S310:
        switch (l)
        {
            case < 1.0e0:
                goto S320;
        }

        qans = c * (w + rt2pin * t / rta);
        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S320:
        ans = c * (w - rt2pin * t / rta);
        qans = 0.5e0 + (0.5e0 - ans);
        return;
        S330:
        switch (a * e * e)
        {
            //
            //  TEMME EXPANSION FOR L = 1
            //
            case > 3.28e-3:
                goto S430;
        }

        c = 0.5e0 + (0.5e0 - y);
        w = (0.5e0 - Math.Sqrt(y) * (0.5e0 + (0.5e0 - y / 3.0e0)) / rtpi) / c;
        u = 1.0e0 / a;
        z = l switch
        {
            < 1.0e0 => -z,
            _ => Math.Sqrt(z + z)
        };

        int T5 = iop - 2;
        switch (T5)
        {
            case < 0:
                goto S340;
            case 0:
                goto S350;
            default:
                goto S360;
        }
        S340:
        c0 = ((((((d0[6] * z + d0[5]) * z + d0[4]) * z + d0[3]) * z + d0[2]) * z + d0[1]) * z + d0[0]) * z -
             third;
        c1 = (((((d1[5] * z + d1[4]) * z + d1[3]) * z + d1[2]) * z + d1[1]) * z + d1[0]) * z + d10;
        c2 = ((((d2[4] * z + d2[3]) * z + d2[2]) * z + d2[1]) * z + d2[0]) * z + d20;
        c3 = (((d3[3] * z + d3[2]) * z + d3[1]) * z + d3[0]) * z + d30;
        c4 = (d4[1] * z + d4[0]) * z + d40;
        c5 = (d5[1] * z + d5[0]) * z + d50;
        c6 = d6[0] * z + d60;
        t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) * u + c0;
        goto S310;
        S350:
        c0 = (d0[1] * z + d0[0]) * z - third;
        c1 = d1[0] * z + d10;
        t = (d20 * u + c1) * u + c0;
        goto S310;
        S360:
        t = d0[0] * z - third;
        goto S310;
        S370:
        //
        //  SPECIAL CASES
        //
        ans = 0.0e0;
        qans = 1.0e0;
        return;
        S380:
        ans = 1.0e0;
        qans = 0.0e0;
        return;
        S390:
        switch (x)
        {
            case >= 0.25e0:
                goto S400;
        }

        double T6 = Math.Sqrt(x);
        ans = error_f(T6);
        qans = 0.5e0 + (0.5e0 - ans);
        return;
        S400:
        double T7 = Math.Sqrt(x);
        qans = error_fc(K2, T7);
        ans = 0.5e0 + (0.5e0 - qans);
        return;
        S410:
        if (Math.Abs(s) <= 2.0e0 * e)
        {
            goto S430;
        }

        S420:
        if (x <= a)
        {
            goto S370;
        }

        goto S380;
        S430:
        //
        //  ERROR RETURN
        //
        ans = 2.0e0;
    }

    public static void gamma_inc_inv(double a, ref double x, ref double x0, ref double p, ref double q,
            ref int ierr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
        //
        //  Discussion:
        //
        //    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
        //    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.
        //    Schroder iteration is employed.  The routine attempts to compute X
        //    to 10 significant digits if this is possible for the particular computer
        //    arithmetic being used.
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
        //    Alfred H Morris, Jr
        //
        //  Parameters:
        //
        //    Input, double *A, the parameter in the incomplete gamma
        //    ratio.  A must be positive.
        //
        //    Output, double *X, the computed point for which the
        //    incomplete gamma functions have the values P and Q.
        //
        //    Input, double *X0, an optional initial approximation
        //    for the solution X.  If the user does not want to supply an
        //    initial approximation, then X0 should be set to 0, or a negative
        //    value.
        //
        //    Input, double *P, *Q, the values of the incomplete gamma
        //    functions, for which the corresponding argument is desired.
        //
        //    Output, int *IERR, error flag.
        //    0, the solution was obtained. Iteration was not used.
        //    0 < K, The solution was obtained. IERR iterations were performed.
        //    -2, A <= 0
        //    -3, No solution was obtained. The ratio Q/A is too large.
        //    -4, P + Q /= 1
        //    -6, 20 iterations were performed. The most recent value obtained
        //        for X is given.  This cannot occur if X0 <= 0.
        //    -7, Iteration failed. No value is given for X.
        //        This may occur when X is approximately 0.
        //    -8, A value for X has been obtained, but the routine is not certain
        //        of its accuracy.  Iteration cannot be performed in this
        //        case. If X0 <= 0, this can occur only when P or Q is
        //        approximately 0. If X0 is positive then this can occur when A is
        //        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
        //
    {
        const double a0 = 3.31125922108741e0;
        const double a1 = 11.6616720288968e0;
        const double a2 = 4.28342155967104e0;
        const double a3 = .213623493715853e0;
        const double b1 = 6.61053765625462e0;
        const double b2 = 6.40691597760039e0;
        const double b3 = 1.27364489782223e0;
        const double b4 = .036117081018842e0;
        const double c = .577215664901533e0;
        const double ln10 = 2.302585e0;
        const double tol = 1e-5;
        double[] amin =  {
                500.0e0,100.0e0
            }
            ;
        double[] bmin =  {
                1e-28,1e-13
            }
            ;
        double[] dmin =  {
                1e-06,1e-04
            }
            ;
        double[] emin =  {
                2e-03,6e-03
            }
            ;
        double[] eps0 =  {
                1e-10,1e-08
            }
            ;
        const int K1 = 1;
        const int K2 = 2;
        const int K3 = 3;
        const int K8 = 0;
        double am1,
            b = 0,
            d,
            h,
            pn = 0,
            qn = 0,
            r,
            s = 0,
            u,
            w,
            y,
            z;

        //
        //  E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
        //            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
        //            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
        //            LARGEST POSITIVE NUMBER.
        //
        double e = dpmpar(K1);
        double xmin = dpmpar(K2);
        double xmax = dpmpar(K3);
        x = 0.0e0;
        switch (a)
        {
            case <= 0.0e0:
                goto S300;
        }

        double t = p + q - 1e0;
        if (Math.Abs(t) > e)
        {
            goto S320;
        }

        ierr = 0;
        switch (p)
        {
            case 0.0e0:
                return;
        }

        switch (q)
        {
            case 0.0e0:
                goto S270;
        }

        switch (a)
        {
            case 1.0e0:
                goto S280;
        }

        double e2 = 2.0e0 * e;
        double amax = 0.4e-10 / (e * e);
        int iop = e switch
        {
            > 1e-10 => 2,
            _ => 1
        };

        double eps = eps0[iop - 1];
        double xn = x0;
        switch (x0)
        {
            case > 0.0e0:
                goto S160;
        }

        switch (a)
        {
            //
            //        SELECTION OF THE INITIAL APPROXIMATION XN OF X
            //                       WHEN A < 1
            //
            case > 1.0e0:
                goto S80;
        }

        double T4 = a + 1.0e0;
        double g = gamma_x(T4);
        double qg = q * g;
        switch (qg)
        {
            case 0.0e0:
                goto S360;
        }

        b = qg / a;
        if (qg > 0.6e0 * a)
        {
            goto S40;
        }

        if (a >= 0.30e0 || b < 0.35e0)
        {
            goto S10;
        }

        t = Math.Exp(-(b + c));
        u = t * Math.Exp(t);
        xn = t * Math.Exp(u);
        goto S160;
        S10:
        switch (b)
        {
            case >= 0.45e0:
                goto S40;
            case 0.0e0:
                goto S360;
        }

        y = -Math.Log(b);
        s = 0.5e0 + (0.5e0 - a);
        z = Math.Log(y);
        t = y - s * z;
        switch (b)
        {
            case < 0.15e0:
                goto S20;
        }

        xn = y - s * Math.Log(t) - Math.Log(1.0e0 + s / (t + 1.0e0));
        goto S220;
        S20:
        switch (b)
        {
            case <= 0.01e0:
                goto S30;
        }

        u = ((t + 2.0e0 * (3.0e0 - a)) * t + (2.0e0 - a) * (3.0e0 - a)) / ((t + (5.0e0 - a)) * t + 2.0e0);
        xn = y - s * Math.Log(t) - Math.Log(u);
        goto S220;
        S30:
        double c1 = -(s * z);
        double c2 = -(s * (1.0e0 + c1));
        double c3 = s * ((0.5e0 * c1 + (2.0e0 - a)) * c1 + (2.5e0 - 1.5e0 * a));
        double c4 = -(s * (((c1 / 3.0e0 + (2.5e0 - 1.5e0 * a)) * c1 + ((a - 6.0e0) * a + 7.0e0)) * c1 + (
            (11.0e0 * a - 46.0) * a + 47.0e0) / 6.0e0));
        double c5 = -(s * ((((-(c1 / 4.0e0) + (11.0e0 * a - 17.0e0) / 6.0e0) * c1 + ((-(3.0e0 * a) + 13.0e0) *
            a - 13.0e0)) * c1 + 0.5e0 * (((2.0e0 * a - 25.0e0) * a + 72.0e0) * a - 61.0e0)) * c1 + ((
            (25.0e0 * a - 195.0e0) * a + 477.0e0) * a - 379.0e0) / 12.0e0));
        xn = (((c5 / y + c4) / y + c3) / y + c2) / y + c1 + y;
        switch (a)
        {
            case > 1.0e0:
                goto S220;
        }

        if (b > bmin[iop - 1])
        {
            goto S220;
        }

        x = xn;
        return;
        S40:
        switch (b * q)
        {
            case > 1e-8:
                goto S50;
        }

        xn = Math.Exp(-(q / a + c));
        goto S70;
        S50:
        switch (p)
        {
            case <= 0.9e0:
                goto S60;
        }

        double T5 = -q;
        xn = Math.Exp((alnrel(T5) + gamma_ln1(a)) / a);
        goto S70;
        S60:
        xn = Math.Exp(Math.Log(p * g) / a);
        S70:
        switch (xn)
        {
            case 0.0e0:
                goto S310;
        }

        t = 0.5e0 + (0.5e0 - xn / (a + 1.0e0));
        xn /= t;
        goto S160;
        S80:
        switch (q)
        {
            //
            //  SELECTION OF THE INITIAL APPROXIMATION XN OF X WHEN A .GT. 1
            //
            case <= 0.5e0:
                goto S90;
        }

        w = Math.Log(p);
        goto S100;
        S90:
        w = Math.Log(q);
        S100:
        t = Math.Sqrt(-(2.0e0 * w));
        s = q switch
        {
            > 0.5e0 => -s,
            _ => t - (((a3 * t + a2) * t + a1) * t + a0) / ((((b4 * t + b3) * t + b2) * t + b1) * t + 1.0e0)
        };

        double rta = Math.Sqrt(a);
        double s2 = s * s;
        xn = a + s * rta + (s2 - 1.0e0) / 3.0e0 + s * (s2 - 7.0e0) / (36.0e0 * rta) - ((3.0e0 * s2 + 7.0e0) *
            s2 - 16.0e0) / (810.0e0 * a) + s * ((9.0e0 * s2 + 256.0e0) * s2 - 433.0e0) / (38880.0e0 * a *
            rta);
        xn = Math.Max(xn, 0.0e0);
        if (a < amin[iop - 1])
        {
            goto S110;
        }

        x = xn;
        d = 0.5e0 + (0.5e0 - x / a);
        if (Math.Abs(d) <= dmin[iop - 1])
        {
            return;
        }

        S110:
        switch (p)
        {
            case <= 0.5e0:
                goto S130;
        }

        if (xn < 3.0e0 * a)
        {
            goto S220;
        }

        y = -(w + gamma_log(a));
        d = Math.Max(2.0e0, a * (a - 1.0e0));
        if (y < ln10 * d)
        {
            goto S120;
        }

        s = 1.0e0 - a;
        z = Math.Log(y);
        goto S30;
        S120:
        t = a - 1.0e0;
        double T6 = -(t / (xn + 1.0e0));
        xn = y + t * Math.Log(xn) - alnrel(T6);
        double T7 = -(t / (xn + 1.0e0));
        xn = y + t * Math.Log(xn) - alnrel(T7);
        goto S220;
        S130:
        double ap1 = a + 1.0e0;
        if (xn > 0.70e0 * ap1)
        {
            goto S170;
        }

        w += gamma_log(ap1);
        if (xn > 0.15e0 * ap1)
        {
            goto S140;
        }

        double ap2 = a + 2.0e0;
        double ap3 = a + 3.0e0;
        x = Math.Exp((w + x) / a);
        x = Math.Exp((w + x - Math.Log(1.0e0 + x / ap1 * (1.0e0 + x / ap2))) / a);
        x = Math.Exp((w + x - Math.Log(1.0e0 + x / ap1 * (1.0e0 + x / ap2))) / a);
        x = Math.Exp((w + x - Math.Log(1.0e0 + x / ap1 * (1.0e0 + x / ap2 * (1.0e0 + x / ap3)))) / a);
        xn = x;
        if (xn > 1e-2 * ap1)
        {
            goto S140;
        }

        if (xn <= emin[iop - 1] * ap1)
        {
            return;
        }

        goto S170;
        S140:
        double apn = ap1;
        t = xn / apn;
        double sum = 1.0e0 + t;
        S150:
        apn += 1.0e0;
        t *= xn / apn;
        sum += t;
        switch (t)
        {
            case > 1e-4:
                goto S150;
        }

        t = w - Math.Log(sum);
        xn = Math.Exp((xn + t) / a);
        xn *= 1.0e0 - (a * Math.Log(xn) - xn - t) / (a - xn);
        goto S170;
        S160:
        switch (p)
        {
            //
            //  SCHRODER ITERATION USING P
            //
            case > 0.5e0:
                goto S220;
        }

        S170:
        if (p <= 1e10 * xmin)
        {
            goto S350;
        }

        am1 = a - 0.5e0 - 0.5e0;
        S180:
        if (a <= amax)
        {
            goto S190;
        }

        d = 0.5e0 + (0.5e0 - xn / a);
        if (Math.Abs(d) <= e2)
        {
            goto S350;
        }

        S190:
        switch (ierr)
        {
            case >= 20:
                goto S330;
        }

        ierr += 1;
        gamma_inc(a, xn, ref pn, ref qn, K8);
        if (pn == 0.0e0 || qn == 0.0e0)
        {
            goto S350;
        }

        r = rcomp(a, xn);
        switch (r)
        {
            case 0.0e0:
                goto S350;
        }

        t = (pn - p) / r;
        w = 0.5e0 * (am1 - xn);
        switch (Math.Abs(t))
        {
            case <= 0.1e0 when Math.Abs(w * t) <= 0.1e0:
                goto S200;
        }

        x = xn * (1.0e0 - t);
        switch (x)
        {
            case <= 0.0e0:
                goto S340;
        }

        d = Math.Abs(t);
        goto S210;
        S200:
        h = t * (1.0e0 + w * t);
        x = xn * (1.0e0 - h);
        switch (x)
        {
            case <= 0.0e0:
                goto S340;
        }

        switch (Math.Abs(w))
        {
            case >= 1.0e0 when Math.Abs(w) * t * t <= eps:
                return;
        }

        d = Math.Abs(h);
        S210:
        xn = x;
        if (d > tol)
        {
            goto S180;
        }

        if (d <= eps)
        {
            return;
        }

        if (Math.Abs(p - pn) <= tol * p)
        {
            return;
        }

        goto S180;
        S220:
        //
        //  SCHRODER ITERATION USING Q
        //
        if (q <= 1e10 * xmin)
        {
            goto S350;
        }

        am1 = a - 0.5e0 - 0.5e0;
        S230:
        if (a <= amax)
        {
            goto S240;
        }

        d = 0.5e0 + (0.5e0 - xn / a);
        if (Math.Abs(d) <= e2)
        {
            goto S350;
        }

        S240:
        switch (ierr)
        {
            case >= 20:
                goto S330;
        }

        ierr += 1;
        gamma_inc(a, xn, ref pn, ref qn, K8);
        if (pn == 0.0e0 || qn == 0.0e0)
        {
            goto S350;
        }

        r = rcomp(a, xn);
        switch (r)
        {
            case 0.0e0:
                goto S350;
        }

        t = (q - qn) / r;
        w = 0.5e0 * (am1 - xn);
        switch (Math.Abs(t))
        {
            case <= 0.1e0 when Math.Abs(w * t) <= 0.1e0:
                goto S250;
        }

        x = xn * (1.0e0 - t);
        switch (x)
        {
            case <= 0.0e0:
                goto S340;
        }

        d = Math.Abs(t);
        goto S260;
        S250:
        h = t * (1.0e0 + w * t);
        x = xn * (1.0e0 - h);
        switch (x)
        {
            case <= 0.0e0:
                goto S340;
        }

        switch (Math.Abs(w))
        {
            case >= 1.0e0 when Math.Abs(w) * t * t <= eps:
                return;
        }

        d = Math.Abs(h);
        S260:
        xn = x;
        if (d > tol)
        {
            goto S230;
        }

        if (d <= eps)
        {
            return;
        }

        if (Math.Abs(q - qn) <= tol * q)
        {
            return;
        }

        goto S230;
        S270:
        //
        //  SPECIAL CASES
        //
        x = xmax;
        return;
        S280:
        switch (q)
        {
            case < 0.9e0:
                goto S290;
        }

        double T9 = -p;
        x = -alnrel(T9);
        return;
        S290:
        x = -Math.Log(q);
        return;
        S300:
        //
        //  ERROR RETURN
        //
        ierr = -2;
        return;
        S310:
        ierr = -3;
        return;
        S320:
        ierr = -4;
        return;
        S330:
        ierr = -6;
        return;
        S340:
        ierr = -7;
        return;
        S350:
        x = xn;
        ierr = -8;
        return;
        S360:
        x = xmax;
        ierr = -8;
    }

    public static void gamma_inc_values(ref int n_data, ref double a, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
        //
        //  Discussion:
        //
        //    The (normalized) incomplete Gamma function P(A,X) is defined as:
        //
        //      PN(A,X) = 1/GAMMA(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
        //
        //    With this definition, for all A and X,
        //
        //      0 <= PN(A,X) <= 1
        //
        //    and
        //
        //      PN(A,INFINITY) = 1.0
        //
        //    Mathematica can compute this value as
        //
        //      1 - GammaRegularized[A,X]
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
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double *A, the parameter of the function.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 20;

        double[] a_vec =  {
                0.1E+00, 0.1E+00, 0.1E+00, 0.5E+00,
                0.5E+00, 0.5E+00, 1.0E+00, 1.0E+00,
                1.0E+00, 1.1E+00, 1.1E+00, 1.1E+00,
                2.0E+00, 2.0E+00, 2.0E+00, 6.0E+00,
                6.0E+00, 11.0E+00, 26.0E+00, 41.0E+00
            }
            ;
        double[] fx_vec =  {
                0.7420263E+00, 0.9119753E+00, 0.9898955E+00, 0.2931279E+00,
                0.7656418E+00, 0.9921661E+00, 0.0951626E+00, 0.6321206E+00,
                0.9932621E+00, 0.0757471E+00, 0.6076457E+00, 0.9933425E+00,
                0.0091054E+00, 0.4130643E+00, 0.9931450E+00, 0.0387318E+00,
                0.9825937E+00, 0.9404267E+00, 0.4863866E+00, 0.7359709E+00
            }
            ;
        double[] x_vec =  {
                3.1622777E-02, 3.1622777E-01, 1.5811388E+00, 7.0710678E-02,
                7.0710678E-01, 3.5355339E+00, 0.1000000E+00, 1.0000000E+00,
                5.0000000E+00, 1.0488088E-01, 1.0488088E+00, 5.2440442E+00,
                1.4142136E-01, 1.4142136E+00, 7.0710678E+00, 2.4494897E+00,
                1.2247449E+01, 1.6583124E+01, 2.5495098E+01, 4.4821870E+01
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            a = 0.0E+00;
            x = 0.0E+00;
            fx = 0.0E+00;
        }
        else
        {
            a = a_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double gamma_ln1(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
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
        //  Input:
        //
        //    double *A, defines the argument of the function.
        //
        //  Output:
        //
        //    double GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
        //
    {
        const double p0 = .577215664901533e+00;
        const double p1 = .844203922187225e+00;
        const double p2 = -.168860593646662e+00;
        const double p3 = -.780427615533591e+00;
        const double p4 = -.402055799310489e+00;
        const double p5 = -.673562214325671e-01;
        const double p6 = -.271935708322958e-02;
        const double q1 = .288743195473681e+01;
        const double q2 = .312755088914843e+01;
        const double q3 = .156875193295039e+01;
        const double q4 = .361951990101499e+00;
        const double q5 = .325038868253937e-01;
        const double q6 = .667465618796164e-03;
        const double r0 = .422784335098467e+00;
        const double r1 = .848044614534529e+00;
        const double r2 = .565221050691933e+00;
        const double r3 = .156513060486551e+00;
        const double r4 = .170502484022650e-01;
        const double r5 = .497958207639485e-03;
        const double s1 = .124313399877507e+01;
        const double s2 = .548042109832463e+00;
        const double s3 = .101552187439830e+00;
        const double s4 = .713309612391000e-02;
        const double s5 = .116165475989616e-03;
        double gamln1, w;

        switch (a)
        {
            case >= 0.6e0:
                goto S10;
        }

        w = ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a + p0) / ((((((q6 * a + q5) * a +
                                                                                        q4) * a + q3) * a + q2) * a + q1) * a + 1.0e0);
        gamln1 = -(a * w);
        return gamln1;
        S10:
        double x = a - 0.5e0 - 0.5e0;
        w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) /
            (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x
             + 1.0e0);
        gamln1 = x * w;
        return gamln1;
    }

    public static double gamma_log(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
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
        //    Alfred H Morris, Jr
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
        //    Input, double *A, the argument of the function.
        //    A should be positive.
        //
        //    Output, double GAMMA_LOG, the value of ln ( Gamma ( A ) ).
        //
    {
        const double c0 = .833333333333333e-01;
        const double c1 = -.277777777760991e-02;
        const double c2 = .793650666825390e-03;
        const double c3 = -.595202931351870e-03;
        const double c4 = .837308034031215e-03;
        const double c5 = -.165322962780713e-02;
        const double d = .418938533204673e0;
        double gamln, t, w;
        int i;

        switch (a)
        {
            case > 0.8e0:
                goto S10;
        }

        gamln = gamma_ln1(a) - Math.Log(a);
        return gamln;
        S10:
        switch (a)
        {
            case > 2.25e0:
                goto S20;
        }

        t = a - 0.5e0 - 0.5e0;
        gamln = gamma_ln1(t);
        return gamln;
        S20:
        switch (a)
        {
            case >= 10.0e0:
                goto S40;
        }

        int n = (int) (a - 1.25e0);
        t = a;
        w = 1.0e0;
        for (i = 1; i <= n; i++)
        {
            t -= 1.0e0;
            w = t * w;
        }

        double T1 = t - 1.0e0;
        gamln = gamma_ln1(T1) + Math.Log(w);
        return gamln;
        S40:
        t = Math.Pow(1.0e0 / a, 2.0);
        w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a;
        gamln = d + w + (a - 0.5e0) * (Math.Log(a) - 1.0e0);
        return gamln;
    }

    public static void gamma_rat1(double a, double x, double r, ref double p, ref double q,
            double eps)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
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
        //    Input, double *A, *X, the parameters of the functions.
        //    It is assumed that A <= 1.
        //
        //    Input, double *R, the value exp(-X) * X**A / Gamma(A).
        //
        //    Output, double *P, *Q, the values of P(A,X) and Q(A,X).
        //
        //    Input, double *EPS, the tolerance.
        //
    {
        double a2n;
        const int K2 = 0;
        double am0 = 0, an = 0, an0 = 0, b2n = 0, b2nm1 = 0, c = 0,
            cma = 0, g = 0, h = 0, j = 0, l = 0, sum = 0, t = 0,
            tol = 0, w = 0, z = 0, T1 = 0, T3 = 0;

        switch (a * x)
        {
            case 0.0e0:
                goto S120;
        }

        switch (a)
        {
            case 0.5e0:
                goto S100;
        }

        switch (x)
        {
            case < 1.1e0:
                goto S10;
        }

        goto S60;
        S10:
        //
        //             TAYLOR SERIES FOR P(A,X)/X**A
        //
        an = 3.0e0;
        c = x;
        sum = x / (a + 3.0e0);
        tol = 0.1e0 * eps / (a + 1.0e0);
        S20:
        an += 1.0e0;
        c = -(c * (x / an));
        t = c / (a + an);
        sum += t;
        if (Math.Abs(t) > tol)
        {
            goto S20;
        }

        j = a * x * ((sum / 6.0e0 - 0.5e0 / (a + 2.0e0)) * x + 1.0e0 / (a + 1.0e0));
        z = a * Math.Log(x);
        h = gam1(a);
        g = 1.0e0 + h;
        switch (x)
        {
            case < 0.25e0:
                goto S30;
        }

        if (a < x / 2.59e0)
        {
            goto S50;
        }

        goto S40;
        S30:
        switch (z)
        {
            case > -.13394e0:
                goto S50;
        }

        S40:
        w = Math.Exp(z);
        p = w * g * (0.5e0 + (0.5e0 - j));
        q = 0.5e0 + (0.5e0 - p);
        return;
        S50:
        l = rexp(z);
        w = 0.5e0 + (0.5e0 + l);
        q = (w * j - l) * g - h;
        switch (q)
        {
            case < 0.0e0:
                goto S90;
        }

        p = 0.5e0 + (0.5e0 - q);
        return;
        S60:
        //
        //  CONTINUED FRACTION EXPANSION
        //
        double a2nm1 = a2n = 1.0e0;
        b2nm1 = x;
        b2n = x + (1.0e0 - a);
        c = 1.0e0;
        S70:
        a2nm1 = x * a2n + c * a2nm1;
        b2nm1 = x * b2n + c * b2nm1;
        am0 = a2nm1 / b2nm1;
        c += 1.0e0;
        cma = c - a;
        a2n = a2nm1 + cma * a2n;
        b2n = b2nm1 + cma * b2n;
        an0 = a2n / b2n;
        if (Math.Abs(an0 - am0) >= eps * an0)
        {
            goto S70;
        }

        q = r * an0;
        p = 0.5e0 + (0.5e0 - q);
        return;
        S80:
        //
        //  SPECIAL CASES
        //
        p = 0.0e0;
        q = 1.0e0;
        return;
        S90:
        p = 1.0e0;
        q = 0.0e0;
        return;
        S100:
        switch (x)
        {
            case >= 0.25e0:
                goto S110;
        }

        T1 = Math.Sqrt(x);
        p = error_f(T1);
        q = 0.5e0 + (0.5e0 - p);
        return;
        S110:
        T3 = Math.Sqrt(x);
        q = error_fc(K2, T3);
        p = 0.5e0 + (0.5e0 - q);
        return;
        S120:
        if (x <= a)
        {
            goto S80;
        }

        goto S90;
    }

    public static void gamma_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_VALUES returns some values of the Gamma function.
        //
        //  Definition:
        //
        //    GAMMA(Z) = Integral ( 0 <= T < Infinity) T^(Z-1) EXP(-T) dT
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
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz and Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    US Department of Commerce, 1964.
        //
        //  Parameters:
        //
        //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double *X, the argument of the function.
        //
        //    Output, double *FX, the value of the function.
        //
    {
        const int N_MAX = 18;

        double[] fx_vec =  {
                4.590845E+00, 2.218160E+00, 1.489192E+00, 1.164230E+00,
                1.0000000000E+00, 0.9513507699E+00, 0.9181687424E+00, 0.8974706963E+00,
                0.8872638175E+00, 0.8862269255E+00, 0.8935153493E+00, 0.9086387329E+00,
                0.9313837710E+00, 0.9617658319E+00, 1.0000000000E+00, 3.6288000E+05,
                1.2164510E+17, 8.8417620E+30
            }
            ;
        double[] x_vec =  {
                0.2E+00, 0.4E+00, 0.6E+00, 0.8E+00,
                1.0E+00, 1.1E+00, 1.2E+00, 1.3E+00,
                1.4E+00, 1.5E+00, 1.6E+00, 1.7E+00,
                1.8E+00, 1.9E+00, 2.0E+00, 10.0E+00,
                20.0E+00, 30.0E+00
            }
            ;

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            x = 0.0E+00;
            fx = 0.0E+00;
        }
        else
        {
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double gamma_x(double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_X evaluates the gamma function.
        //
        //  Discussion:
        //
        //    This routine was renamed from "GAMMA" to avoid a conflict with the
        //    C/C++ math library routine.
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
        //    Alfred H Morris, Jr
        //
        //  Input:
        //
        //    double *A, the argument of the Gamma function.
        //
        //  Output:
        //
        //    double GAMMA_X, the value of the Gamma function.
        //
    {
        const double d = .41893853320467274178e0;
        const double r1 = .820756370353826e-03;
        const double r2 = -.595156336428591e-03;
        const double r3 = .793650663183693e-03;
        const double r4 = -.277777777770481e-02;
        const double r5 = .833333333333333e-01;
        double[] p =  {
                .539637273585445e-03,.261939260042690e-02,.204493667594920e-01,
                .730981088720487e-01,.279648642639792e+00,.553413866010467e+00,1.0e0
            }
            ;
        double[] q =  {
                -.832979206704073e-03,.470059485860584e-02,.225211131035340e-01,
                -.170458969313360e+00,-.567902761974940e-01,.113062953091122e+01,1.0e0
            }
            ;
        const int K2 = 3;
        const int K3 = 0;
        double s = 0, t = 0;
        int i, j;

        double Xgamm = 0.0e0;
        double x = a;
        switch (Math.Abs(a))
        {
            case >= 15.0e0:
                goto S110;
        }

        //
        //  EVALUATION OF GAMMA(A) FOR ABS(A) < 15
        //
        t = 1.0e0;
        int m = (int)Math.Truncate(a) - 1;
        //
        //  LET T BE THE PRODUCT OF A-J WHEN A .GE. 2
        //
        int T1 = m;
        switch (T1)
        {
            case < 0:
                goto S40;
            case 0:
                goto S30;
            default:
                goto S10;
        }
        S10:
        for (j = 1; j <= m; j++)
        {
            x -= 1.0e0;
            t = x * t;
        }

        S30:
        x -= 1.0e0;
        goto S80;
        S40:
        //
        //  LET T BE THE PRODUCT OF A+J WHEN A < 1
        //
        t = a;
        switch (a)
        {
            case > 0.0e0:
                goto S70;
        }

        m = -m - 1;
        switch (m)
        {
            case 0:
                goto S60;
        }

        for (j = 1; j <= m; j++)
        {
            x += 1.0e0;
            t = x * t;
        }

        S60:
        x += 0.5e0 + 0.5e0;
        t = x * t;
        switch (t)
        {
            case 0.0e0:
                return Xgamm;
        }

        S70:
        switch (Math.Abs(t))
        {
            //
            //  THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
            //  CODE MAY BE OMITTED IF DESIRED.
            //
            case >= 1e-30:
                goto S80;
        }

        if (Math.Abs(t) * dpmpar(K2) <= 1.0001e0)
        {
            return Xgamm;
        }

        Xgamm = 1.0e0 / t;
        return Xgamm;
        S80:
        //
        //  COMPUTE GAMMA(1 + X) FOR  0 <= X < 1
        //
        double top = p[0];
        double bot = q[0];
        for (i = 1; i < 7; i++)
        {
            top = p[i] + x * top;
            bot = q[i] + x * bot;
        }

        Xgamm = top / bot;
        switch (a)
        {
            //
            //  TERMINATION
            //
            case < 1.0e0:
                goto S100;
        }

        Xgamm *= t;
        return Xgamm;
        S100:
        Xgamm /= t;
        return Xgamm;
        S110:
        switch (Math.Abs(a))
        {
            //
            //  EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15
            //
            case >= 1e3:
                return Xgamm;
        }

        switch (a)
        {
            case > 0.0e0:
                goto S120;
        }

        x = -a;
        int n = (int) x;
        t = t switch
        {
            > 0.9e0 => 1.0e0 - t,
            _ => x - n
        };

        s = (n % 2) switch
        {
            0 => -s,
            _ => Math.Sin(Math.PI * t) / Math.PI
        };

        switch (s)
        {
            case 0.0e0:
                return Xgamm;
        }

        S120:
        //
        //  COMPUTE THE MODIFIED ASYMPTOTIC SUM
        //
        t = 1.0e0 / (x * x);
        double g = ((((r1 * t + r2) * t + r3) * t + r4) * t + r5) / x;
        //
        //  ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X)
        //  BUT LESS ACCURACY WILL NORMALLY BE OBTAINED.
        //
        double lnx = Math.Log(x);
        //
        //  FINAL ASSEMBLY
        //
        double z = x;
        g = d + g + (z - 0.5e0) * (lnx - 1e0);
        double w = g;
        t = g - w;
        if (w > 0.99999e0 * exparg(K3))
        {
            return Xgamm;
        }

        Xgamm = a switch
        {
            < 0.0e0 => 1.0e0 / (Xgamm * s) / x,
            _ => Math.Exp(w) * (1.0e0 + t)
        };

        return Xgamm;
    }

}