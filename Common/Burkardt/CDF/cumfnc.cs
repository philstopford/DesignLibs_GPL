using System;

namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void cumfnc(double f, double dfn, double dfd, double pnonc,
            ref double cum, ref double ccum)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUMFNC evaluates the cumulative noncentral F distribution.
        //
        //  Discussion:
        //
        //    This routine computes the noncentral F distribution with DFN and DFD
        //    degrees of freedom and noncentrality parameter PNONC.
        //
        //    The series is calculated backward and forward from J = LAMBDA/2
        //    (this is the term with the largest Poisson weight) until
        //    the convergence criterion is met.
        //
        //    The sum continues until a succeeding term is less than EPS
        //    times the sum (or the sum is less than 1.0e-20).  EPS is
        //    set to 1.0e-4 in a data statement which can be changed.
        //
        //    The original version of this routine allowed the input values
        //    of DFN and DFD to be negative (nonsensical) or zero (which
        //    caused numerical overflow.)  I have forced both these values
        //    to be at least 1.
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
        //    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
        //
        //  Parameters:
        //
        //    Input, double *F, the upper limit of integration.
        //
        //    Input, double *DFN, *DFD, the number of degrees of freedom
        //    in the numerator and denominator.  Both DFN and DFD must be positive,
        //    and normally would be integers.  This routine requires that they
        //    be no less than 1.
        //
        //    Input, double *PNONC, the noncentrality parameter.
        //
        //    Output, double *CUM, *CCUM, the noncentral F CDF and
        //    complementary CDF.
        //
    {
        double half = 0.5e0;
        double done = 1.0e0;

        double adn = 0;
        double aup = 0;
        double b = 0;
        double betdn = 0;
        double betup = 0;
        double centwt = 0;
        double dnterm = 0;
        double dsum = 0;
        double dummy = 0;
        double eps = 1.0e-4;
        double prod = 0;
        double xx = 0;
        double yy = 0;

        double sum = 0,
            upterm = 0,
            xmult = 0,
            xnonc = 0;
        int i = 0, icent = 0, ierr = 0;
        double T1 = 0, T2 = 0, T3 = 0, T4 = 0, T5 = 0, T6 = 0;

        bool qsmall(double x)
        {
            return sum < 1.0e-20 || x < eps * sum;
        }

        switch (f)
        {
            case <= 0.0e0:
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
        }

        switch ((pnonc < 1.0e-10))
        {
            case false:
                goto S20;
        }

        //
        //  Handle case in which the non-centrality parameter is
        //  (essentially) zero.
        //
        cumf(f, dfn, dfd, ref cum, ref ccum);
        return;
        S20:
        xnonc = pnonc / 2.0e0;
        icent = icent switch
        {
            0 => 1,
            //
            //  Calculate the central term of the poisson weighting factor.
            //
            _ => (int) xnonc
        };

        //
        //  Compute central weight term
        //
        T1 = icent + 1;
        centwt = Math.Exp(-xnonc + icent * Math.Log(xnonc) - gamma_log(T1));
        //
        //  Compute central incomplete beta term
        //  Assure that minimum of arg to beta and 1 - arg is computed
        //  accurately.
        //
        prod = dfn * f;
        dsum = dfd + prod;
        yy = dfd / dsum;
        if (yy > half)
        {
            xx = prod / dsum;
            yy = done - xx;
        }
        else
        {
            xx = done - yy;
        }

        T2 = dfn * half + icent;
        T3 = dfd * half;
        beta_inc(T2, T3, xx, yy, ref betdn, ref dummy, ref ierr);
        adn = dfn / 2.0e0 + icent;
        aup = adn;
        b = dfd / 2.0e0;
        betup = betdn;
        sum = centwt * betdn;
        //
        //  Now sum terms backward from icent until convergence or all done
        //
        xmult = centwt;
        i = icent;
        T4 = adn + b;
        T5 = adn + 1.0e0;
        dnterm = Math.Exp(gamma_log(T4) - gamma_log(T5)
                                        - gamma_log(b) + adn * Math.Log(xx) + b * Math.Log(yy));
        S30:
        if (qsmall(xmult * betdn) || i <= 0)
        {
            goto S40;
        }

        xmult *= (i / xnonc);
        i -= 1;
        adn -= 1.0;
        dnterm = (adn + 1.0) / ((adn + b) * xx) * dnterm;
        betdn += dnterm;
        sum += xmult * betdn;
        goto S30;
        S40:
        i = icent + 1;
        //
        //  Now sum forwards until convergence
        //
        xmult = centwt;
        switch (aup - 1.0 + b)
        {
            case 0:
                upterm = Math.Exp(-gamma_log(aup)
                                  - gamma_log(b) + (aup - 1.0) * Math.Log(xx) +
                                  b * Math.Log(yy));
                break;
            default:
                T6 = aup - 1.0 + b;
                upterm = Math.Exp(gamma_log(T6) - gamma_log(aup)
                                                - gamma_log(b) + (aup - 1.0) * Math.Log(xx) + b *
                    Math.Log(yy));
                break;
        }

        goto S60;
        S50:
        if (qsmall(xmult * betup))
        {
            goto S70;
        }

        S60:
        xmult *= (xnonc / i);
        i += 1;
        aup += 1.0;
        upterm = (aup + b - 2.0e0) * xx / (aup - 1.0) * upterm;
        betup -= upterm;
        sum += xmult * betup;
        goto S50;
        S70:
        cum = sum;
        ccum = 0.5e0 + (0.5e0 - cum);
    }
}