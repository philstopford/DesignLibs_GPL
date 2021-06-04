using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cumbet(double x, double y, double a, double b, ref double cum,
                ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMBET evaluates the cumulative incomplete beta distribution.
            //
            //  Discussion:
            //
            //    This routine calculates the CDF to X of the incomplete beta distribution
            //    with parameters A and B.  This is the integral from 0 to x
            //    of (1/B(a,b))f(t)) where f(t) = t^(a-1) * (1-t)^(b-1)
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
            //    A R Didonato and Alfred Morris,
            //    Algorithm 708:
            //    Significant Digit Computation of the Incomplete Beta Function Ratios.
            //    ACM Transactions on Mathematical Software,
            //    Volume 18, Number 3, September 1992, pages 360-373.
            //
            //  Input:
            //
            //    double *X, the upper limit of integration.
            //
            //    double *Y, the value of 1-X.
            //
            //    double *A, *B, the parameters of the distribution.
            //
            //  Output:
            //
            //    double *CUM, *CCUM, the values of the cumulative
            //    density function and complementary cumulative density function.
            //
        {
            int ierr = 0;

            if (x <= 0.0)
            {
                cum = 0.0;
                ccum = 1.0;
            }
            else if (y <= 0.0)
            {
                cum = 1.0;
                ccum = 0.0;
            }
            else
            {
                beta_inc(a, b, x, y, ref cum, ref ccum, ref ierr);
            }

        }

        public static void cumbin(double s, double xn, double pr, double ompr,
                ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMBIN evaluates the cumulative binomial distribution.
            //
            //  Discussion:
            //
            //    This routine returns the probability of 0 to S successes in XN binomial
            //    trials, each of which has a probability of success, PR.
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
            //    1966, Formula 26.5.24.
            //
            //  Parameters:
            //
            //    Input, double *S, the upper limit of summation.
            //
            //    Input, double *XN, the number of trials.
            //
            //    Input, double *PR, the probability of success in one trial.
            //
            //    Input, double *OMPR, equals ( 1 - PR ).
            //
            //    Output, double *CUM, the cumulative binomial distribution.
            //
            //    Output, double *CCUM, the complement of the cumulative
            //    binomial distribution.
            //
        {
            double T1;
            double T2;

            if (s < xn)
            {
                T1 = s + 1.0;
                T2 = xn - s;
                cumbet(pr, ompr, T1, T2, ref ccum, ref cum);
            }
            else
            {
                cum = 1.0;
                ccum = 0.0;
            }

        }

        public static void cumchi(double x, double df, ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMCHI evaluates the cumulative chi-square distribution.
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
            //    Input, double *X, the upper limit of integration.
            //
            //    Input, double *DF, the degrees of freedom of the
            //    chi-square distribution.
            //
            //    Output, double *CUM, the cumulative chi-square distribution.
            //
            //    Output, double *CCUM, the complement of the cumulative
            //    chi-square distribution.
            //
        {
            double a;
            double xx;

            a = df * 0.5;
            xx = x * 0.5;
            cumgam(xx, a, ref cum, ref ccum);

        }

        public static void cumchn(double x, double df, double pnonc, ref double cum,
                ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMCHN evaluates the cumulative noncentral chi-square distribution.
            //
            //  Discussion:
            //
            //    Calculates the cumulative noncentral chi-square
            //    distribution, i.e., the probability that a random variable
            //    which follows the noncentral chi-square distribution, with
            //    noncentrality parameter PNONC and continuous degrees of
            //    freedom DF, is less than or equal to X.
            //
            //    Thanks to Tongzhou Wang, who corrected the line
            //      "sumadj = sum + adj;"
            //    to
            //      "sumadj = sumadj + adj;"
            //    13 February 2021
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
            //    1966, Formula 26.4.25.
            //
            //  Input:
            //
            //    double *X, the upper limit of integration.
            //
            //    double *DF, the number of degrees of freedom.
            //
            //    double *PNONC, the noncentrality parameter of
            //    the noncentral chi-square distribution.
            //
            //  Output:
            //
            //    double *CUM, *CCUM, the CDF and complementary
            //    CDF of the noncentral chi-square distribution.
            //
            //  Local:
            //
            //    double EPS, the convergence criterion.  The sum
            //    stops when a term is less than EPS*SUM.
            //
            //    int NTIRED, the maximum number of terms to be evaluated
            //    in each sum.
            //
            //    bool QCONV, is TRUE if convergence was achieved, that is,
            //    the program did not stop on NTIRED criterion.
            //
        {
            double adj;
            double centaj;
            double centwt;
            double chid2;
            double dfd2;
            double eps = 1.0e-5;
            int i;
            int icent;
            int iterb;
            int iterf;
            double lcntaj;
            double lcntwt;
            double lfact;
            int ntired = 1000;
            double pcent = 0;
            double pterm;
            double sum;
            double sumadj;
            double term;
            double wt;
            double xnonc;
            double T1;
            double T2;
            double T3;

            if (x <= 0.0e0)
            {
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
            }

            double dg(double i)
            {
                return (df + 2.0e0 * (double) (i));
            }

            bool qsmall(double xx)
            {
                return (sum < 1.0e-20 || (xx) < eps*sum);
            }

            bool qtired(double i)
            {
                return ((i) > ntired);
            }
            
            //
            //  When non-centrality parameter is (essentially) zero,
            //  use cumulative chi-square distribution
            //
            if (pnonc <= 1.0e-10)
            {
                cumchi(x, df, ref cum, ref ccum);
                return;
            }

            xnonc = pnonc / 2.0e0;
            //
            //  The following code calculates the weight, chi-square, and
            //  adjustment term for the central term in the infinite series.
            //  The central term is the one in which the poisson weight is
            //  greatest.  The adjustment term is the amount that must
            //  be subtracted from the chi-square to move up two degrees
            //  of freedom.
            //
            icent = (int)(xnonc);
            if (icent == 0)
            {
                icent = 1;
            }

            chid2 = x / 2.0e0;
            //
            //  Calculate central weight term
            //
            T1 = (double) (icent + 1);
            lfact = gamma_log(T1);
            lcntwt = -xnonc + (double) icent * Math.Log(xnonc) - lfact;
            centwt = Math.Exp(lcntwt);
            //
            //  Calculate central chi-square
            //
            T2 = dg(icent);
            cumchi(x, T2, ref pcent, ref ccum);
            //
            //  Calculate central adjustment term
            //
            dfd2 = dg(icent) / 2.0e0;
            T3 = 1.0e0 + dfd2;
            lfact = gamma_log(T3);
            lcntaj = dfd2 * Math.Log(chid2) - chid2 - lfact;
            centaj = Math.Exp(lcntaj);
            sum = centwt * pcent;
            //
            //  Sum backwards from the central term towards zero.
            //  Quit whenever either
            //  (1) the zero term is reached, or
            //  (2) the term gets small relative to the sum, or
            //  (3) More than NTIRED terms are totaled.
            //
            iterb = 0;
            sumadj = 0.0e0;
            adj = centaj;
            wt = centwt;
            i = icent;
            goto S40;
            S30:
            if (qtired(iterb) || qsmall(term) || i == 0) goto S50;
            S40:
            dfd2 = dg(i) / 2.0e0;
            //
            //  Adjust chi-square for two fewer degrees of freedom.
            //  The adjusted value ends up in PTERM.
            //
            adj = adj * dfd2 / chid2;
            sumadj = sumadj + adj;
            pterm = pcent + sumadj;
            //
            //  Adjust poisson weight for J decreased by one
            //
            wt = wt * ((double) i / xnonc);
            term = wt * pterm;
            sum = sum + term;
            i = i - 1;
            iterb = iterb + 1;
            goto S30;
            S50:
            iterf = 0;
            //
            //  Now sum forward from the central term towards infinity.
            //  Quit when either
            //  (1) the term gets small relative to the sum, or
            //  (2) More than NTIRED terms are totaled.
            //
            sumadj = adj = centaj;
            wt = centwt;
            i = icent;
            goto S70;
            S60:
            if (qtired(iterf) || qsmall(term)) goto S80;
            S70:
            //
            //  Update weights for next higher J
            //
            wt = wt * (xnonc / (double) (i + 1));
            //
            //  Calculate PTERM and add term to sum
            //
            pterm = pcent - sumadj;
            term = wt * pterm;
            sum = sum + term;
            //
            //  Update adjustment term for DF for next iteration
            //
            i = i + 1;
            dfd2 = dg(i) / 2.0e0;
            adj = adj * chid2 / dfd2;
            //
            //  Previously, this line read "sumadj = sum + adj;"
            //  Correction submitted by Tongzhou Wang, 13 February 2021.
            //
            sumadj = sumadj + adj;
            iterf = iterf + 1;
            goto S60;
            S80:
            cum = sum;
            ccum = 0.5e0 + (0.5e0 - cum);
            return;
        }

        public static void cumf(double f, double dfn, double dfd, ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMF evaluates the cumulative F distribution.
            //
            //  Discussion:
            //
            //    CUMF computes the integral from 0 to F of the F density with DFN
            //    numerator and DFD denominator degrees of freedom.
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
            //    1966, Formula 26.5.28.
            //
            //  Parameters:
            //
            //    Input, double *F, the upper limit of integration.
            //
            //    Input, double *DFN, *DFD, the number of degrees of
            //    freedom for the numerator and denominator.
            //
            //    Output, double *CUM, *CCUM, the value of the F CDF and
            //    the complementary F CDF.
            //
        {
            double half = 0.5e0;
            double done = 1.0e0;

            double dsum;
            int ierr = 0;
            double prod;
            double T1;
            double T2;
            double xx;
            double yy;

            if (f <= 0.0e0)
            {
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
            }

            prod = dfn * f;
            //
            //  XX is such that the incomplete beta with parameters
            //  DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
            //  YY is 1 - XX
            //  Calculate the smaller of XX and YY accurately
            //
            dsum = dfd + prod;
            xx = dfd / dsum;

            if (xx > half)
            {
                yy = prod / dsum;
                xx = done - yy;
            }
            else
            {
                yy = done - xx;
            }

            T1 = dfd * half;
            T2 = dfn * half;
            beta_inc(T1, T2, xx, yy, ref ccum, ref cum, ref ierr);
        }

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

            double adn;
            double aup;
            double b;
            double betdn = 0;
            double betup;
            double centwt;
            double dnterm;
            double dsum;
            double dummy = 0;
            double eps = 1.0e-4;
            double prod;
            double xx;
            double yy;

            double sum,
                upterm,
                xmult,
                xnonc;
            int i, icent, ierr = 0;
            double T1, T2, T3, T4, T5, T6;

            if (f <= 0.0e0)
            {
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
            }

            bool qsmall(double x)
            {
                return (sum < 1.0e-20 || (x) < eps*sum);
            }

            if (!(pnonc < 1.0e-10)) goto S20;
            //
            //  Handle case in which the non-centrality parameter is
            //  (essentially) zero.
            //
            cumf(f, dfn, dfd, ref cum, ref ccum);
            return;
            S20:
            xnonc = pnonc / 2.0e0;
            //
            //  Calculate the central term of the poisson weighting factor.
            //
            icent = (int) xnonc;
            if (icent == 0) icent = 1;
            //
            //  Compute central weight term
            //
            T1 = (double) (icent + 1);
            centwt = Math.Exp(-xnonc + (double) icent * Math.Log(xnonc) - gamma_log(T1));
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
            else xx = done - yy;

            T2 = dfn * half + (double) icent;
            T3 = dfd * half;
            beta_inc(T2, T3, xx, yy, ref betdn, ref dummy, ref ierr);
            adn = dfn / 2.0e0 + (double) icent;
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
            if (qsmall(xmult * betdn) || i <= 0) goto S40;
            xmult = xmult * ((double) i / xnonc);
            i = i - 1;
            adn = adn - 1.0;
            dnterm = (adn + 1.0) / ((adn + b) * xx) * dnterm;
            betdn = betdn + dnterm;
            sum = sum + (xmult * betdn);
            goto S30;
            S40:
            i = icent + 1;
            //
            //  Now sum forwards until convergence
            //
            xmult = centwt;
            if (aup - 1.0 + b == 0)
                upterm = Math.Exp(-gamma_log(aup)
                             - gamma_log(b) + (aup - 1.0) * Math.Log(xx) +
                             b * Math.Log(yy));
            else
            {
                T6 = aup - 1.0 + b;
                upterm = Math.Exp(gamma_log(T6) - gamma_log(aup)
                                            - gamma_log(b) + (aup - 1.0) * Math.Log(xx) + b *
                    Math.Log(yy));
            }

            goto S60;
            S50:
            if (qsmall(xmult * betup)) goto S70;
            S60:
            xmult = xmult * (xnonc / (double) i);
            i = i + 1;
            aup = aup + 1.0;
            upterm = (aup + b - 2.0e0) * xx / (aup - 1.0) * upterm;
            betup = betup - upterm;
            sum = sum + (xmult * betup);
            goto S50;
            S70:
            cum = sum;
            ccum = 0.5e0 + (0.5e0 - cum);
            return;
        }

        public static void cumgam(double x, double a, ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMGAM evaluates the cumulative incomplete gamma distribution.
            //
            //  Discussion:
            //
            //    This routine computes the cumulative distribution function of the
            //    incomplete gamma distribution, i.e., the integral from 0 to X of
            //
            //      (1/GAM(A))*EXP(-T)*T^(A-1) DT
            //
            //    where GAM(A) is the complete gamma function of A, i.e.,
            //
            //      GAM(A) = integral from 0 to infinity of EXP(-T)*T^(A-1) DT
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
            //    Input, double *X, the upper limit of integration.
            //
            //    Input, double *A, the shape parameter of the incomplete
            //    Gamma distribution.
            //
            //    Output, double *CUM, *CCUM, the incomplete Gamma CDF and
            //    complementary CDF.
            //
        {
            int K1 = 0;

            if (x <= 0.0e0)
            {
                cum = 0.0e0;
                ccum = 1.0e0;
            }
            else
            {
                gamma_inc(a, x, ref cum, ref ccum, K1);
            }
        }

        public static void cumnbn(double s, double xn, double pr, double ompr,
                ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMNBN evaluates the cumulative negative binomial distribution.
            //
            //  Discussion:
            //
            //    This routine returns the probability that there will be F or
            //    fewer failures before there are S successes, with each binomial
            //    trial having a probability of success PR.
            //
            //    Prob(# failures = F | S successes, PR)  =
            //                        ( S + F - 1 )
            //                        (            ) * PR^S * (1-PR)^F
            //                        (      F     )
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
            //    1966, Formula 26.5.26.
            //
            //  Parameters:
            //
            //    Input, double *F, the number of failures.
            //
            //    Input, double *S, the number of successes.
            //
            //    Input, double *PR, *OMPR, the probability of success on
            //    each binomial trial, and the value of (1-PR).
            //
            //    Output, double *CUM, *CCUM, the negative binomial CDF,
            //    and the complementary CDF.
            //
        {
            double T1;

            T1 = s + 1e0;
            cumbet(pr, ompr, xn, T1, ref cum, ref ccum);
        }

        public static void cumnor(double arg, ref double result, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMNOR computes the cumulative normal distribution.
            //
            //  Discussion:
            //
            //    This function evaluates the normal distribution function:
            //
            //                              / x
            //                     1       |       -tt/2
            //          P(x) = ----------- |      e       dt
            //                 sqrt(2 pi)  |
            //                             /-oo
            //
            //    This transportable program uses rational functions that
            //    theoretically approximate the normal distribution function to
            //    at least 18 significant decimal digits.  The accuracy achieved
            //    depends on the arithmetic system, the compiler, the intrinsic
            //    functions, and proper selection of the machine-dependent
            //    constants.
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
            //    William Cody
            //
            //  Reference:
            //
            //    William Cody,
            //    Rational Chebyshev approximations for the error function,
            //    Mathematics of Computation,
            //    1969, pages 631-637.
            //
            //    William Cody,
            //    Algorithm 715:
            //    SPECFUN - A Portable FORTRAN Package of Special Function Routines
            //      and Test Drivers,
            //    ACM Transactions on Mathematical Software,
            //    Volume 19, 1993, pages 22-32.
            //
            //  Parameters:
            //
            //    Input, double *ARG, the upper limit of integration.
            //
            //    Output, double *CUM, *CCUM, the Normal density CDF and
            //    complementary CDF.
            //
            //  Local:
            //
            //    double EPS, the argument below which anorm(x)
            //    may be represented by 0.5D+00 and above which  xx  will not underflow.
            //    A conservative value is the largest machine number X
            //    such that   1.0D+00 + X = 1.0D+00   to machine precision.
            //
        {
            double[] a =  {
                2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
                1.8154981253343561249e04,6.5682337918207449113e-2
            }
            ;
            double[] b =  {
                4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
                4.5507789335026729956e04
            }
            ;
            double[] c =  {
                3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
                5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
                1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
            }
            ;
            double[] d =  {
                2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
                6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
                3.8912003286093271411e04,1.9685429676859990727e04
            }
            ;
            double half = 0.5e0;
            double[] p =  {
                2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
                1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
            }
            ;
            double one = 1.0e0;
            double[] q =  {
                1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
                3.78239633202758244e-3,7.29751555083966205e-5
            }
            ;
            double sixten = 1.60e0;
            double sqrpi = 3.9894228040143267794e-1;
            double thrsh = 0.66291e0;
            double root32 = 5.656854248e0;
            double zero = 0.0e0;
            int K1 = 1;
            int K2 = 2;
            int i;
            double del, eps, temp, x, xden, xnum, y, xsq, min;
            //
            //  Machine dependent constants
            //
            eps = dpmpar(K1) * 0.5e0;
            min = dpmpar(K2);
            x = arg;
            y = Math.Abs(x);
            if (y <= thrsh)
            {
                //
                //  Evaluate  anorm  for  |X| <= 0.66291
                //
                xsq = zero;
                if (y > eps) xsq = x * x;
                xnum = a[4] * xsq;
                xden = xsq;
                for (i = 0; i < 3; i++)
                {
                    xnum = (xnum + a[i]) * xsq;
                    xden = (xden + b[i]) * xsq;
                }

                result = x * (xnum + a[3]) / (xden + b[3]);
                temp = result;
                result = half + temp;
                ccum = half - temp;
            }
            //
            //  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
            //
            else if (y <= root32)
            {
                xnum = c[8] * y;
                xden = y;
                for (i = 0; i < 7; i++)
                {
                    xnum = (xnum + c[i]) * y;
                    xden = (xden + d[i]) * y;
                }

                result = (xnum + c[7]) / (xden + d[7]);
                xsq = (double)(int)(y * sixten) / sixten;
                del = (y - xsq) * (y + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (x > zero)
                {
                    temp = result;
                    result = ccum;
                    ccum = temp;
                }
            }
            //
            //  Evaluate  anorm  for |X| > sqrt(32)
            //
            else
            {
                result = zero;
                xsq = one / (x * x);
                xnum = p[5] * xsq;
                xden = xsq;
                for (i = 0; i < 4; i++)
                {
                    xnum = (xnum + p[i]) * xsq;
                    xden = (xden + q[i]) * xsq;
                }

                result = xsq * (xnum + p[4]) / (xden + q[4]);
                result = (sqrpi - result) / y;
                xsq = (double)(int)(x * sixten) / sixten;
                del = (x - xsq) * (x + xsq);
                result = Math.Exp(-(xsq * xsq * half)) * Math.Exp(-(del * half)) * result;
                ccum = one - result;
                if (x > zero)
                {
                    temp = result;
                    result = ccum;
                    ccum = temp;
                }
            }

            if (result < min) result = 0.0e0;
            //
            //  Fix up for negative argument, erf, etc.
            //
            if (ccum < min) ccum = 0.0e0;
        }

        public static void cumpoi(double s, double xlam, ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMPOI evaluates the cumulative Poisson distribution.
            //
            //  Discussion:
            //
            //    CUMPOI returns the probability of S or fewer events in a Poisson
            //    distribution with mean XLAM.
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
            //    Handbook of Mathematical Functions,
            //    Formula 26.4.21.
            //
            //  Parameters:
            //
            //    Input, double *S, the upper limit of cumulation of the
            //    Poisson density function.
            //
            //    Input, double *XLAM, the mean of the Poisson distribution.
            //
            //    Output, double *CUM, *CCUM, the Poisson density CDF and
            //    complementary CDF.
            //
        {
            double chi;
            double df;

            df = 2.0e0 * (s + 1.0e0);
            chi = 2.0e0 * xlam;
            cumchi(chi, df, ref ccum, ref cum);

            return;
        }

        public static void cumt(double t, double df, ref double cum, ref double ccum)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUMT evaluates the cumulative T distribution.
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
            //    Handbook of Mathematical Functions,
            //    Formula 26.5.27.
            //
            //  Parameters:
            //
            //    Input, double *T, the upper limit of integration.
            //
            //    Input, double *DF, the number of degrees of freedom of
            //    the T distribution.
            //
            //    Output, double *CUM, *CCUM, the T distribution CDF and
            //    complementary CDF.
            //
        {
            double a = 0;
            double dfptt;
            double K2 = 0.5e0;
            double oma = 0;
            double T1;
            double tt;
            double xx;
            double yy;

            tt = (t) * (t);
            dfptt = (df) + tt;
            xx = df / dfptt;
            yy = tt / dfptt;
            T1 = 0.5e0 * (df);
            cumbet(xx, yy, T1, K2, ref a, ref oma);

            if (t <= 0.0e0)
            {
                cum = 0.5e0 * a;
                ccum = oma + (cum);
            }
            else
            {
                ccum = 0.5e0 * a;
                cum = oma + (ccum);
            }

        }

    }
}