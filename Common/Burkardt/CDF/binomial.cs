using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfbin(int which, ref double p, ref double q, ref double s, ref double xn,
                ref double pr, ref double ompr, ref int status, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFBIN evaluates the CDF of the Binomial distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the binomial distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    P is the probablility of S or fewer successes in XN binomial trials,
            //    each trial having an individual probability of success of PR.
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
            //    Input, int *WHICH, indicates which of argument values is to
            //    be calculated from the others.
            //    1: Calculate P and Q from S, XN, PR and OMPR;
            //    2: Calculate S from P, Q, XN, PR and OMPR;
            //    3: Calculate XN from P, Q, S, PR and OMPR;
            //    4: Calculate PR and OMPR from P, Q, S and XN.
            //
            //    Input/output, double *P, the cumulation, from 0 to S,
            //    of the binomial distribution.  If P is an input value, it should
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *S, the number of successes observed.
            //    Whether this is an input or output value, it should lie in the
            //    range [0,XN].
            //
            //    Input/output, double *XN, the number of binomial trials.
            //    If this is an input value it should lie in the range: (0, +infinity).
            //    If it is an output value it will be searched for in the
            //    range [1.0D-300, 1.0D+300].
            //
            //    Input/output, double *PR, the probability of success in each
            //    binomial trial.  Whether this is an input or output value, it should
            //    lie in the range: [0,1].
            //
            //    Input/output, double *OMPR, equal to 1-PR.  Whether this is an
            //    input or output value, it should lie in the range [0,1].  Also, it should
            //    be the case that PR + OMPR = 1.
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +4, if PR + OMPR /= 1.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double atol = (1.0e-50);
            double tol = (1.0e-8);
            double zero = (1.0e-300);
            double inf = 1.0e300;
            double one = 1.0e0;

            double ccum = 0;
            double cum = 0;
            double fx = 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K3 = 0.5e0;
            double K4 = 5.0e0;
            double K11 = 1.0e0;
            double pq;
            double prompr = 0;
            bool qhi = false;
            bool qleft = false;
            bool qporq = false;
            double T5;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T12;
            double T13;
            double xhi = 0;
            double xlo = 0;

            status = 0;
            bound = 0.0;

            E0000Data e0000Data = new E0000Data();
            E0001Data e0001Data = new E0001Data();
            
            //
            //  Check arguments
            //
            if (!(which < 1 && which > 4)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 4.0e0;
            S20:
            status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //     P
            //
            if (!(p < 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = 1.0e0;
            S50:
            status = -2;
            return;
            S70:
            S60:
            if (which == 1) goto S110;
            //
            //     Q
            //
            if (!(q < 0.0e0 || q > 1.0e0)) goto S100;
            if (!(q < 0.0e0)) goto S80;
            bound = 0.0e0;
            goto S90;
            S80:
            bound = 1.0e0;
            S90:
            status = -3;
            return;
            S110:
            S100:
            if (which == 3) goto S130;
            //
            //     XN
            //
            if (!(xn <= 0.0e0)) goto S120;
            bound = 0.0e0;
            status = -5;
            return;
            S130:
            S120:
            if (which == 2) goto S170;
            //
            //     S
            //
            if (!(s < 0.0e0 || (which != 3 && s > xn))) goto S160;
            if (!(s < 0.0e0)) goto S140;
            bound = 0.0e0;
            goto S150;
            S140:
            bound = xn;
            S150:
            status = -4;
            return;
            S170:
            S160:
            if (which == 4) goto S210;
            //
            //     PR
            //
            if (!(pr < 0.0e0 || pr > 1.0e0)) goto S200;
            if (!(pr < 0.0e0)) goto S180;
            bound = 0.0e0;
            goto S190;
            S180:
            bound = 1.0e0;
            S190:
            status = -6;
            return;
            S210:
            S200:
            if (which == 4) goto S250;
            //
            //     OMPR
            //
            if (!(ompr < 0.0e0 || ompr > 1.0e0)) goto S240;
            if (!(ompr < 0.0e0)) goto S220;
            bound = 0.0e0;
            goto S230;
            S220:
            bound = 1.0e0;
            S230:
            status = -7;
            return;
            S250:
            S240:
            if (which == 1) goto S290;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S280;
            if (!(pq < 0.0e0)) goto S260;
            bound = 0.0e0;
            goto S270;
            S260:
            bound = 1.0e0;
            S270:
            status = 3;
            return;
            S290:
            S280:
            if (which == 4) goto S330;
            //
            //     PR + OMPR
            //
            prompr = pr + ompr;
            if (!(Math.Abs(prompr - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S320;
            if (!(prompr < 0.0e0)) goto S300;
            bound = 0.0e0;
            goto S310;
            S300:
            bound = 1.0e0;
            S310:
            status = 4;
            return;
            S330:
            S320:
            if (!(which == 1)) qporq = p <= q;
            //
            //     Select the minimum of P or Q
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P
                //
                cumbin(s, xn, pr, ompr, ref p, ref q);
                status = 0;
            }
            else if (2 == which)
            {
                //
                //  Calculating S
                //
                s = 5.0e0;
                T5 = atol;
                T6 = tol;
                e0000Data.zsmall = K2;
                e0000Data.zbig = xn;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T5;
                e0000Data.zrelto = T6;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S340:
                if (!(status == 1)) goto S370;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S350;
                fx = cum - p;
                goto S360;
                S350:
                fx = ccum - q;
                S360:
                e0000Data.status = status;
                dinvr(ref e0000Data);
                goto S340;
                S370:
                if (!(status == -1)) goto S400;
                if (!qleft) goto S380;
                status = 1;
                bound = 0.0e0;
                goto S390;
                S380:
                status = 2;
                bound = xn;
                S400:
                S390: ;
            }
            else if (3 == which)
            {
                //
                //  Calculating XN
                //
                xn = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                e0000Data.zsmall = T7;
                e0000Data.zbig = T8;
                e0000Data.zabsst = K3;
                e0000Data.zrelst = K3;
                e0000Data.zstpmu = K4;
                e0000Data.zabsto = T9;
                e0000Data.zrelto = T10;
                dstinv(ref e0000Data);
                status = 0;
                e0000Data.status = status;
                dinvr(ref e0000Data);
                S410:
                if (!(status == 1)) goto S440;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                if (!qporq) goto S420;
                fx = cum - p;
                goto S430;
                S420:
                fx = ccum - q;
                S430:
                e0000Data.status = status;
                dinvr(ref e0000Data);
                goto S410;
                S440:
                if (!(status == -1)) goto S470;
                if (!qleft) goto S450;
                status = 1;
                bound = zero;
                goto S460;
                S450:
                status = 2;
                bound = inf;
                S470:
                S460: ;
            }
            else if (4 == which)
            {
                //
                //  Calculating PR and OMPR
                //
                T12 = atol;
                T13 = tol;
                e0001Data.zxlo = K2;
                e0001Data.zxhi = K11;
                e0001Data.zabstl = T12;
                e0001Data.zreltl = T13;
                dstzr(ref e0001Data);
                if (!qporq) goto S500;
                status = 0;
                dzror(ref e0000Data, ref e0001Data);
                ompr = one - pr;
                S480:
                if (!(status == 1)) goto S490;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                fx = cum - p;
                dzror(ref e0000Data, ref e0001Data);
                ompr = one - pr;
                goto S480;
                S490:
                goto S530;
                S500:
                status = 0;
                dzror(ref e0000Data, ref e0001Data);
                pr = one - ompr;
                S510:
                if (!(status == 1)) goto S520;
                cumbin(s, xn, pr, ompr, ref cum, ref ccum);
                fx = ccum - q;
                dzror(ref e0000Data, ref e0001Data);
                pr = one - ompr;
                goto S510;
                S530:
                S520:
                if (!(status == -1)) goto S560;
                if (!qleft) goto S540;
                status = 1;
                bound = 0.0e0;
                goto S550;
                S540:
                status = 2;
                bound = 1.0e0;
                S550: ;
            }

            S560:
            return;
        }

        public static void negative_binomial_cdf_values(ref int n_data, ref int f, ref int s, ref double p,
                ref double cdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    negative_binomial_cdf_values() returns values of the negative binomial CDF.
            //
            //  Discussion:
            //
            //    Assume that a coin has a probability P of coming up heads on
            //    any one trial.  Suppose that we plan to flip the coin until we
            //    achieve a total of S heads.  If we let F represent the number of
            //    tails that occur in this process, then the value of F satisfies
            //    a negative binomial PDF:
            //
            //      PDF(F,S,P) = Choose ( F from F+S-1 ) * P^S * (1-P)^F
            //
            //    The negative binomial CDF is the probability that there are F or
            //    fewer failures upon the attainment of the S-th success.  Thus,
            //
            //      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
            //
            //    Thanks to Sebastiono Vigna for reporting an error in which the
            //    expression "if ( n_data < 0 )" was used, 05 May 2021.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    F C Powell,
            //    Statistical Tables for Sociology, Biology and Physical Sciences,
            //    Cambridge University Press, 1982.
            //
            //  Input:
            //
            //    int *N_DATA.  The user sets N_DATA to 0 before the first call.
            //
            //  Output:
            //
            //    int *N_DATA.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    int *F, the maximum number of failures.
            //
            //    int *S, the number of successes.
            //
            //    double *P, the probability of a success on one trial.
            //
            //    double *CDF, the probability of at most F failures before the
            //    S-th success.
            //
        {
            int N_MAX = 27;

            double[] cdf_vec =
                {
                    0.6367, 0.3633, 0.1445,
                    0.5000, 0.2266, 0.0625,
                    0.3438, 0.1094, 0.0156,
                    0.1792, 0.0410, 0.0041,
                    0.0705, 0.0109, 0.0007,
                    0.9862, 0.9150, 0.7472,
                    0.8499, 0.5497, 0.2662,
                    0.6513, 0.2639, 0.0702,
                    1.0000, 0.0199, 0.0001
                }
                ;
            int[] f_vec =
                {
                    4, 3, 2,
                    3, 2, 1,
                    2, 1, 0,
                    2, 1, 0,
                    2, 1, 0,
                    11, 10, 9,
                    17, 16, 15,
                    9, 8, 7,
                    2, 1, 0
                }
                ;
            double[] p_vec =
                {
                    0.50, 0.50, 0.50,
                    0.50, 0.50, 0.50,
                    0.50, 0.50, 0.50,
                    0.40, 0.40, 0.40,
                    0.30, 0.30, 0.30,
                    0.30, 0.30, 0.30,
                    0.10, 0.10, 0.10,
                    0.10, 0.10, 0.10,
                    0.01, 0.01, 0.01
                }
                ;
            int[] s_vec =
                {
                    4, 5, 6,
                    4, 5, 6,
                    4, 5, 6,
                    4, 5, 6,
                    4, 5, 6,
                    1, 2, 3,
                    1, 2, 3,
                    1, 2, 3,
                    0, 1, 2
                }
                ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                f = 0;
                s = 0;
                p = 0.0E+00;
                cdf = 0.0E+00;
            }
            else
            {
                f = f_vec[n_data - 1];
                s = s_vec[n_data - 1];
                p = p_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
            }
        }

        public static void binomial_cdf_values(ref int n_data, ref int a, ref double b, ref int x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A,B) is the probability of at most X successes in A trials,
            //    given that the probability of success on a single trial is B.
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int *A, double *B, the parameters of the function.
            //
            //    Output, int *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 17;

            int[] a_vec =  {
                2, 2, 2, 2,
                2, 4, 4, 4,
                4, 10, 10, 10,
                10, 10, 10, 10,
                10
            }
            ;
            double[] b_vec =  {
                0.05E+00, 0.05E+00, 0.05E+00, 0.50E+00,
                0.50E+00, 0.25E+00, 0.25E+00, 0.25E+00,
                0.25E+00, 0.05E+00, 0.10E+00, 0.15E+00,
                0.20E+00, 0.25E+00, 0.30E+00, 0.40E+00,
                0.50E+00
            }
            ;
            double[] fx_vec =  {
                0.9025E+00, 0.9975E+00, 1.0000E+00, 0.2500E+00,
                0.7500E+00, 0.3164E+00, 0.7383E+00, 0.9492E+00,
                0.9961E+00, 0.9999E+00, 0.9984E+00, 0.9901E+00,
                0.9672E+00, 0.9219E+00, 0.8497E+00, 0.6331E+00,
                0.3770E+00
            }
            ;
            int[] x_vec =  {
                0, 1, 2, 0,
                1, 0, 1, 2,
                3, 4, 4, 4,
                4, 4, 4, 4,
                4
            }
            ;

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0.0E+00;
                x = 0;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}