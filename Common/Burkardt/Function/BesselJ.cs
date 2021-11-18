using System;

namespace Burkardt.Function;

public static class BesselJ
{
    public static double besselj ( double order, double x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BESSELJ evaluates the Bessel J function at an arbitrary real order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ORDER, the order.
        //    0.0 <= ORDER.
        //
        //    Input, double, X, the evaluation point. 
        //
        //    Output, double BESSELJ, the value.
        //
    {
        double alpha;
        double[] b;
        int n;
        int nb;
        int ncalc = 0;
        double value = 0;

        n = ( int ) order;
        nb = n + 1;
        alpha = order - n;
        b = new double[nb];

        rjbesl ( x, alpha, nb, ref b, ref ncalc );

        value = b[n];

        return value;
    }
    public static void bessel_jx_values(ref int n_data, ref double nu, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BESSEL_JX_VALUES returns some values of the Jx Bessel function.
        //
        //  Discussion:
        //
        //    This set of data considers the less common case in which the
        //    index of the Bessel function Jn is actually not an integer.
        //    We may suggest this case by occasionally replacing the symbol
        //    "Jn" by "Jx".
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      BesselJ[n,x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &NU, the order of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 28;

        double[] fx_vec =
            {
                0.3544507442114011E+00,
                0.6713967071418031E+00,
                0.5130161365618278E+00,
                0.3020049060623657E+00,
                0.06500818287737578E+00,
                -0.3421679847981618E+00,
                -0.1372637357550505E+00,
                0.1628807638550299E+00,
                0.2402978391234270E+00,
                0.4912937786871623E+00,
                -0.1696513061447408E+00,
                0.1979824927558931E+00,
                -0.1094768729883180E+00,
                0.04949681022847794E+00,
                0.2239245314689158E+00,
                0.2403772011113174E+00,
                0.1966584835818184E+00,
                0.02303721950962553E+00,
                0.3314145508558904E+00,
                0.5461734240402840E+00,
                -0.2616584152094124E+00,
                0.1296035513791289E+00,
                -0.1117432171933552E+00,
                0.03142623570527935E+00,
                0.1717922192746527E+00,
                0.3126634069544786E+00,
                0.1340289119304364E+00,
                0.06235967135106445E+00
            }
            ;
        double[] nu_vec =
            {
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                1.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                2.50E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                1.25E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00,
                2.75E+00
            }
            ;

        double[] x_vec =
            {
                0.2E+00,
                1.0E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00,
                1.0E+00,
                2.0E+00,
                5.0E+00,
                10.0E+00,
                50.0E+00
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
            nu = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            nu = nu_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }

    }

    public static void rjbesl(double x, double alpha, int nb, ref double[] b, ref int ncalc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RJBESL evaluates a sequence of Bessel J functions.
        //
        //  Discussion:
        //
        //    This routine calculates Bessel functions J sub(N+ALPHA) (X)
        //    for non-negative argument X, and non-negative order N+ALPHA.
        //
        //    This program is based on a program written by David Sookne
        //    that computes values of the Bessel functions J or I of real
        //    argument and integer order.  Modifications include the restriction
        //    of the computation to the J Bessel function of non-negative real
        //    argument, the extension of the computation to arbitrary positive
        //    order, and the elimination of most underflow.
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    Original FORTRAN77 version by William Cody.
        //    C++ version by John Burkardt.
        //
        //  Reference: 
        //
        //    F Olver, David Sookne,
        //    A Note on Backward Recurrence Algorithms," 
        //    Math. Comp.,
        //    Volume 26, 1972, pages 941-947.
        //
        //    David Sookne,
        //    Bessel Functions of Real Argument and Integer Order,
        //    NBS Journal of Res. B,
        //    Volume 77B, 1973, pages 125-132.
        //
        //  Parameters:
        //
        //    Input, double X, non-negative real argument for which
        //    J's are to be calculated.
        //
        //    Input, double ALPHA, fractional part of order for which
        //    J's or exponentially scaled J'r (J*exp(X)) are
        //    to be calculated.  0 <= ALPHA < 1.0.
        //
        //    Input, int NB, number of functions to be calculated, 
        //    NB > 0.  The first function calculated is of order ALPHA, and the
        //    last is of order (NB - 1 + ALPHA).
        //
        //    Output, double B(NB).  If RJBESL
        //    terminates normally (NCALC=NB), the vector B contains the
        //    functions J/ALPHA/(X) through J/NB-1+ALPHA/(X), or the
        //    corresponding exponentially scaled functions.
        //
        //    Output, int &NCALC, indicates possible errors.
        //    Before using the vector B, the user should check that
        //    NCALC=NB, i.e., all orders have been calculated to
        //    the desired accuracy.  See Error Returns below.
        //
        //  Internal Parameters:
        //
        //    IT = Number of bits in the mantissa of a working precision variable
        //
        //    NSIG   = Decimal significance desired.  Should be set to
        //    INT(LOG10(2)*it+1).  Setting NSIG lower will result
        //    in decreased accuracy while setting NSIG higher will
        //    increase CPU time without increasing accuracy.  The
        //    truncation error is limited to a relative error of
        //    T=.5*10**(-NSIG).
        //
        //    Then the following machine-dependent constants must be declared
        //    in DATA statements.  IEEE values are provided as a default.
        //
        //    ENTEN  = 10.0 ** K, where K is the largest integer such that
        //    ENTEN is machine-representable in working precision.
        //
        //    ENSIG  = 10.0 ** NSIG
        //
        //    RTNSIG = 10.0 ** (-K) for the smallest integer K such that K >= NSIG/4
        //
        //    ENMTEN = Smallest ABS(X) such that X/4 does not underflow
        //
        //    XLARGE = Upper limit on the magnitude of X.  If ABS(X)=N,
        //    then at least N iterations of the backward recursion
        //    will be executed.  The value of 10.0 ** 4 is used on
        //    every machine.
        //
        //  Error returns:
        //
        //    In case of an error,  NCALC != NB, and not all J's are
        //    calculated to the desired accuracy.
        //
        //    NCALC < 0:  An argument is out of range. For example,
        //    NBES <= 0, ALPHA < 0 or > 1, or X is too large.
        //    In this case, B(1) is set to zero, the remainder of the
        //    B-vector is not calculated, and NCALC is set to
        //    MIN(NB,0)-1 so that NCALC != NB.
        //
        //    NB > NCALC > 0: Not all requested function values could
        //    be calculated accurately.  This usually occurs because NB is
        //    much larger than ABS(X).  In this case, B(N) is calculated
        //    to the desired accuracy for N <= NCALC, but precision
        //    is lost for NCALC < N <= NB.  If B(N) does not vanish
        //    for N > NCALC (because it is too small to be represented),
        //    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
        //    significant figures of B(N) can be trusted.
        //
    {
        double alpem;
        double alp2em;
        double capp;
        double capq;
        const double eighth = 0.125E+00;
        double em;
        double en;
        const double enmten = 8.90E-308;
        const double ensig = 1.0E+16;
        const double enten = 1.0E+308;
        double[] fact = {
                1.0E+00,
                1.0E+00,
                2.0E+00,
                6.0E+00,
                24.0E+00,
                1.2E+02,
                7.2E+02,
                5.04E+03,
                4.032E+04,
                3.6288E+05,
                3.6288E+06,
                3.99168E+07,
                4.790016E+08,
                6.2270208E+09,
                8.71782912E+10,
                1.307674368E+12,
                2.0922789888E+13,
                3.55687428096E+14,
                6.402373705728E+15,
                1.21645100408832E+17,
                2.43290200817664E+18,
                5.109094217170944E+19,
                1.1240007277776076E+21,
                2.585201673888497664E+22,
                6.2044840173323943936E+23
            }
            ;
        const double four = 4.0E+00;
        double gnu;
        const double half = 0.5E+00;
        double halfx;
        int i;
        int j;
        bool jump;
        int k;
        int l;
        int m;
        int magx;
        int n;
        int nbmx;
        int nend;
        int nstart;
        const double one = 1.0E+00;
        const double one30 = 130.0E+00;
        double p;
        const double pi2 = 0.636619772367581343075535E+00;
        double plast;
        double pold;
        double psave;
        double psavel;
        const double rtnsig = 1.0E-04;
        double s;
        double sum;
        double t;
        double t1;
        double tempa;
        double tempb;
        double tempc;
        double test;
        const double three = 3.0E+00;
        const double three5 = 35.0E+00;
        double tover;
        const double two = 2.0E+00;
        const double twofiv = 25.0E+00;
        const double twopi1 = 6.28125E+00;
        const double twopi2 = 1.935307179586476925286767E-03;
        double xc;
        double xin;
        double xk;
        const double xlarge = 1.0E+04;
        double xm;
        double vcos;
        double vsin;
        double z;
        const double zero = 0.0E+00;

        jump = false;
        //
        //  Check for out of range arguments.
        //
        magx = (int) x;

        switch (nb)
        {
            case > 0 when zero <= x && x <= xlarge && zero <= alpha && alpha < one:
            {
                //
                //  Initialize result array to zero.
                //
                ncalc = nb;
                for (i = 1; i <= nb; i++)
                {
                    b[i - 1] = zero;
                }

                switch (x)
                {
                    //
                    //  Branch to use 2-term ascending series for small X and asymptotic
                    //  form for large X when NB is not too large.
                    //
                    case < rtnsig:
                    {
                        //
                        //  Two-term ascending series for small X.
                        //
                        tempa = one;
                        alpem = one + alpha;

                        halfx = x switch
                        {
                            > enmten => half * x,
                            _ => zero
                        };

                        if (Math.Abs(alpha - zero) > double.Epsilon)
                        {
                            tempa = Math.Pow(halfx, alpha) / (alpha * Helpers.Gamma(alpha));
                        }

                        tempb = (x + one) switch
                        {
                            > one => -halfx * halfx,
                            _ => zero
                        };

                        b[0] = tempa + tempa * tempb / alpem;

                        if (Math.Abs(x - zero) > double.Epsilon && Math.Abs(b[0] - zero) <= double.Epsilon)
                        {
                            ncalc = 0;
                        }

                        if (nb != 1)
                        {
                            switch (x)
                            {
                                case <= zero:
                                {
                                    for (n = 2; n <= nb; n++)
                                    {
                                        b[n - 1] = zero;
                                    }

                                    break;
                                }
                                //
                                default:
                                {
                                    tempc = halfx;

                                    if (Math.Abs(tempb - zero) > double.Epsilon)
                                    {
                                        tover = enmten / tempb;
                                    }
                                    else
                                    {
                                        tover = (enmten + enmten) / x;
                                    }

                                    for (n = 2; n <= nb; n++)
                                    {
                                        tempa /= alpem;
                                        alpem += one;

                                        tempa *= tempc;
                                        if (tempa <= tover * alpem)
                                        {
                                            tempa = zero;
                                        }

                                        b[n - 1] = tempa + tempa * tempb / alpem;

                                        ncalc = b[n - 1] switch
                                        {
                                            zero when n < ncalc => n - 1,
                                            _ => ncalc
                                        };
                                    }

                                    break;
                                }
                            }
                        }

                        break;
                    }
                    //
                    //  Asymptotic series for 21.0 < X.
                    //
                    case > twofiv when nb <= magx + 1:
                    {
                        xc = Math.Sqrt(pi2 / x);
                        xin = Math.Pow(eighth / x, 2);

                        m = x switch
                        {
                            < three5 => 11,
                            < one30 => 8,
                            _ => 4
                        };

                        xm = four * m;
                        //
                        //  Argument reduction for SIN and COS routines.
                        //
                        t = Math.Truncate(x / (twopi1 + twopi2) + half);
                        z = x - t * twopi1 - t * twopi2 - (alpha + half) / pi2;
                        vsin = Math.Sin(z);
                        vcos = Math.Cos(z);
                        gnu = alpha + alpha;

                        for (i = 1; i <= 2; i++)
                        {
                            s = (xm - one - gnu) * (xm - one + gnu) * xin * half;
                            t = (gnu - (xm - three)) * (gnu + (xm - three));
                            capp = s * t / fact[2 * m];
                            t1 = (gnu - (xm + one)) * (gnu + (xm + one));
                            capq = s * t1 / fact[2 * m + 1];
                            xk = xm;
                            k = m + m;
                            t1 = t;

                            for (j = 2; j <= m; j++)
                            {
                                xk -= four;
                                s = (xk - one - gnu) * (xk - one + gnu);
                                t = (gnu - (xk - three)) * (gnu + (xk - three));
                                capp = (capp + one / fact[k - 2]) * s * t * xin;
                                capq = (capq + one / fact[k - 1]) * s * t1 * xin;
                                k -= 2;
                                t1 = t;
                            }

                            capp += one;
                            capq = (capq + one) * (gnu * gnu - one) * (eighth / x);
                            b[i - 1] = xc * (capp * vcos - capq * vsin);

                            switch (nb)
                            {
                                case 1:
                                    return;
                            }

                            t = vsin;
                            vsin = -vcos;
                            vcos = t;
                            gnu += two;
                        }

                        switch (nb)
                        {
                            //
                            //  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1
                            //
                            case > 2:
                            {
                                gnu = alpha + alpha + two;
                                for (j = 3; j <= nb; j++)
                                {
                                    b[j - 1] = gnu * b[j - 2] / x - b[j - 3];
                                    gnu += two;
                                }

                                break;
                            }
                        }

                        break;
                    }
                    //
                    default:
                    {
                        nbmx = nb - magx;
                        n = magx + 1;
                        en = n + n + (alpha + alpha);
                        plast = one;
                        p = en / x;
                        //
                        //  Calculate general significance test.
                        //
                        test = ensig + ensig;
                        switch (nbmx)
                        {
                            //
                            //  Calculate P*S until N = NB-1.  Check for possible overflow.
                            //
                            case >= 3:
                            {
                                tover = enten / ensig;
                                nstart = magx + 2;
                                nend = nb - 1;
                                en = nstart + nstart - two + (alpha + alpha);

                                for (k = nstart; k <= nend; k++)
                                {
                                    n = k;
                                    en += two;
                                    pold = plast;
                                    plast = p;
                                    p = en * plast / x - pold;
                                    //
                                    //  To avoid overflow, divide P*S by TOVER.  Calculate P*S until 1 < ABS(P).
                                    //
                                    if (tover < p)
                                    {
                                        tover = enten;
                                        p /= tover;
                                        plast /= tover;
                                        psave = p;
                                        psavel = plast;
                                        nstart = n + 1;

                                        for (;;)
                                        {
                                            n += 1;
                                            en += two;
                                            pold = plast;
                                            plast = p;
                                            p = en * plast / x - pold;
                                            if (one < p)
                                            {
                                                break;
                                            }
                                        }

                                        tempb = en / x;
                                        //
                                        //  Calculate backward test and find NCALC, the highest N such that
                                        //  the test is passed.
                                        //
                                        test = pold * plast * (half - half / (tempb * tempb));
                                        test /= ensig;
                                        p = plast * tover;
                                        n -= 1;
                                        en -= two;
                                        if (nb < n)
                                        {
                                            nend = nb;
                                        }
                                        else
                                        {
                                            nend = n;
                                        }

                                        for (l = nstart; l <= nend; l++)
                                        {
                                            pold = psavel;
                                            psavel = psave;
                                            psave = en * psavel / x - pold;
                                            if (test < psave * psavel)
                                            {
                                                ncalc = l - 1;
                                                jump = true;
                                                break;
                                            }
                                        }

                                        if (jump)
                                        {
                                            break;
                                        }

                                        ncalc = nend;
                                        jump = true;
                                        break;
                                    }
                                }

                                switch (jump)
                                {
                                    case false:
                                    {
                                        n = nend;
                                        en = n + n + (alpha + alpha);
                                        //
                                        //  Calculate special significance test for 2 < NBMX.
                                        //
                                        if (test < Math.Sqrt(plast * ensig) * Math.Sqrt(p + p))
                                        {
                                            test = Math.Sqrt(plast * ensig) * Math.Sqrt(p + p);
                                        }

                                        break;
                                    }
                                }

                                break;
                            }
                        }

                        switch (jump)
                        {
                            //
                            //  Calculate P*S until significance test passes.
                            //
                            case false:
                            {
                                for (;;)
                                {
                                    n += 1;
                                    en += two;
                                    pold = plast;
                                    plast = p;
                                    p = en * plast / x - pold;

                                    if (test <= p)
                                    {
                                        break;
                                    }
                                }

                                break;
                            }
                        }

                        //
                        //  Initialize the backward recursion and the normalization sum.
                        //
                        n += 1;
                        en += two;
                        tempb = zero;
                        tempa = one / p;
                        m = 2 * n - 4 * (n / 2);
                        sum = zero;
                        em = (double)n / 2;
                        alpem = em - one + alpha;
                        alp2em = em + em + alpha;
                        if (m != 0)
                        {
                            sum = tempa * alpem * alp2em / em;
                        }

                        nend = n - nb;
                        switch (nend)
                        {
                            //
                            //  Recur backward via difference equation, calculating (but not
                            //  storing) B, until N = NB.
                            //
                            case > 0:
                            {
                                for (l = 1; l <= nend; l++)
                                {
                                    n -= 1;
                                    en -= two;
                                    tempc = tempb;
                                    tempb = tempa;
                                    tempa = en * tempb / x - tempc;
                                    m = 2 - m;

                                    if (m != 0)
                                    {
                                        em -= one;
                                        alp2em = em + em + alpha;
                                        if (n == 1)
                                        {
                                            break;
                                        }

                                        alpem = alpem switch
                                        {
                                            zero => one,
                                            _ => em - one + alpha
                                        };

                                        sum = (sum + tempa * alp2em) * alpem / em;
                                    }
                                }

                                break;
                            }
                        }

                        //
                        //  Store B[NB-1].
                        //
                        b[n - 1] = tempa;

                        switch (nend)
                        {
                            case >= 0 when nb <= 1:
                            {
                                alp2em = (alpha + one) switch
                                {
                                    one => one,
                                    _ => alpha
                                };

                                sum += b[0] * alp2em;

                                if (Math.Abs(alpha + one - one) > double.Epsilon)
                                {
                                    sum = sum * Helpers.Gamma(alpha) * Math.Pow(x * half, -alpha);
                                }

                                tempa = enmten;

                                switch (sum)
                                {
                                    case > one:
                                        tempa *= sum;
                                        break;
                                }

                                for (n = 1; n <= nb; n++)
                                {
                                    if (Math.Abs(b[n - 1]) < tempa)
                                    {
                                        b[n - 1] = zero;
                                    }

                                    b[n - 1] /= sum;
                                }

                                return;
                            }
                            //
                            //  Calculate and store B[NB-2].
                            //
                            case >= 0:
                            {
                                n -= 1;
                                en -= two;
                                b[n - 1] = en * tempa / x - tempb;

                                switch (n)
                                {
                                    case 1:
                                    {
                                        em -= one;
                                        alp2em = alp2em switch
                                        {
                                            zero => one,
                                            _ => em + em + alpha
                                        };

                                        sum += b[0] * alp2em;
                                        //
                                        //  Normalize.  Divide all B by sum.
                                        //
                                        if (Math.Abs(alpha + one - one) > double.Epsilon)
                                        {
                                            sum = sum * Helpers.Gamma(alpha) * Math.Pow(x * half, -alpha);
                                        }

                                        tempa = enmten;

                                        switch (sum)
                                        {
                                            case > one:
                                                tempa *= sum;
                                                break;
                                        }

                                        for (n = 1; n <= nb; n++)
                                        {
                                            if (Math.Abs(b[n - 1]) < tempa)
                                            {
                                                b[n - 1] = zero;
                                            }

                                            b[n - 1] /= sum;
                                        }

                                        return;
                                    }
                                }

                                m = 2 - m;

                                if (m != 0)
                                {
                                    em -= one;
                                    alp2em = em + em + alpha;
                                    alpem = alpem switch
                                    {
                                        zero => one,
                                        _ => em - one + alpha
                                    };

                                    sum = (sum + b[n - 1] * alp2em) * alpem / em;
                                }

                                break;
                            }
                        }

                        nend = n - 2;
                        //
                        //  Calculate via difference equation and store B, until N = 2.
                        //
                        if (nend != 0)
                        {
                            for (l = 1; l <= nend; l++)
                            {
                                n -= 1;
                                en -= two;
                                b[n - 1] = en * b[n] / x - b[n + 1];
                                m = 2 - m;

                                if (m != 0)
                                {
                                    em -= one;
                                    alp2em = em + em + alpha;
                                    alpem = alpem switch
                                    {
                                        zero => one,
                                        _ => em - one + alpha
                                    };

                                    sum = (sum + b[n - 1] * alp2em) * alpem / em;
                                }
                            }
                        }

                        //
                        //  Calculate B[0].
                        //
                        b[0] = two * (alpha + one) * b[1] / x - b[2];

                        em -= one;
                        alp2em = alp2em switch
                        {
                            zero => one,
                            _ => em + em + alpha
                        };

                        sum += b[0] * alp2em;
                        //
                        //  Normalize.  Divide all B by sum.
                        //
                        if (Math.Abs(alpha + one - one) > double.Epsilon)
                        {
                            sum = sum * Helpers.Gamma(alpha) * Math.Pow(x * half, -alpha);
                        }

                        tempa = enmten;

                        switch (sum)
                        {
                            case > one:
                                tempa *= sum;
                                break;
                        }

                        for (n = 1; n <= nb; n++)
                        {
                            if (Math.Abs(b[n - 1]) < tempa)
                            {
                                b[n - 1] = zero;
                            }

                            b[n - 1] /= sum;
                        }

                        break;
                    }
                }

                break;
            }
            //
            default:
            {
                b[0] = zero;
                ncalc = nb switch
                {
                    < 0 => nb - 1,
                    _ => -1
                };

                break;
            }
        }
    }

}