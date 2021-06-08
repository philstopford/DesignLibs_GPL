using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
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

            double adj = 0;
            double centaj = 0;
            double centwt = 0;
            double chid2 = 0;
            double dfd2 = 0;
            double eps = 1.0e-5;
            int i = 0;
            int icent = 0;
            int iterb = 0;
            int iterf = 0;
            double lcntaj = 0;
            double lcntwt = 0;
            double lfact = 0;
            int ntired = 1000;
            double pcent = 0;
            double pterm = 0;
            double sum;
            double sumadj = 0;
            double term = 0;
            double wt = 0;
            double xnonc = 0;
            double T1 = 0;
            double T2 = 0;
            double T3 = 0;

            double dg(int i_)
            {
                return (df + 2.0e0 * (double) (i_));
            }

            bool qsmall(double xx)
            {
                return (sum < 1.0e-20 || (xx) < eps * sum);
            }

            bool qtired(int i_)
            {
                return ((i_) > ntired);
            }

            if (x <= 0.0e0)
            {
                cum = 0.0e0;
                ccum = 1.0e0;
                return;
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
        }
    }
}