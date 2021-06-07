using System;

namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void cdfchi(int which, ref double p, ref double q, ref double x, ref double df,
                ref int status_, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFCHI evaluates the CDF of the chi square distribution.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the chi square distribution
            //    given the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The CDF of the chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF [ ChiSquareDistribution [ DF ], X ]
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
            //    1966, Formula 26.4.19.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X and DF;
            //    2: Calculate X from P, Q and DF;
            //    3: Calculate DF from P, Q and X.
            //
            //    Input/output, double *P, the integral from 0 to X of
            //    the chi-square distribution.  If this is an input value, it should
            //    lie in the range [0,1].
            //
            //    Input/output, double *Q, equal to 1-P.  If Q is an input
            //    value, it should lie in the range [0,1].  If Q is an output value,
            //    it will lie in the range [0,1].
            //
            //    Input/output, double *X, the upper limit of integration
            //    of the chi-square distribution.  If this is an input
            //    value, it should lie in the range: [0, +infinity).  If it is an output
            //    value, it will be searched for in the range: [0,1.0D+300].
            //
            //    Input/output, double *DF, the degrees of freedom of the
            //    chi-square distribution.  If this is an input value, it should lie
            //    in the range: (0, +infinity).  If it is an output value, it will be
            //    searched for in the range: [ 1.0D-300, 1.0D+300].
            //
            //    Output, int *STATUS, reports the status of the computation.
            //     0, if the calculation completed correctly;
            //    -I, if the input parameter number I is out of range;
            //    +1, if the answer appears to be lower than lowest search bound;
            //    +2, if the answer appears to be higher than greatest search bound;
            //    +3, if P + Q /= 1;
            //    +10, an error was returned from CUMGAM.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double inf = 1.0e300;

            double ccum = 0;
            double cum = 0;
            int K1 = 1;
            double K2 = 0.0e0;
            double K4 = 0.5e0;
            double K5 = 5.0e0;
            double porq = 0;
            double pq;
            bool qporq = false;
            double T3;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11;

            bound = 0.0;

            E0000E0001 eData = new E0000E0001();
            eData.status = 0;

            
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 3)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 3.0e0;
            S20:
            eData.status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //  P
            //
            if (!(p < 0.0e0 || p > 1.0e0)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = 1.0e0;
            S50:
            eData.status = -2;
            return;
            S70:
            S60:
            if (which == 1) goto S110;
            //
            //     Q
            //
            if (!(q <= 0.0e0 || q > 1.0e0)) goto S100;
            if (!(q <= 0.0e0)) goto S80;
            bound = 0.0e0;
            goto S90;
            S80:
            bound = 1.0e0;
            S90:
            eData.status = -3;
            return;
            S110:
            S100:
            if (which == 2) goto S130;
            //
            //     X
            //
            if (!(x < 0.0e0)) goto S120;
            bound = 0.0e0;
            eData.status = -4;
            return;
            S130:
            S120:
            if (which == 3) goto S150;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S140;
            bound = 0.0e0;
            eData.status = -5;
            return;
            S150:
            S140:
            if (which == 1) goto S190;
            //
            //     P + Q
            //
            pq = p + q;
            if (!(Math.Abs(pq - 0.5e0 - 0.5e0) > 3.0e0 * dpmpar(K1))) goto S180;
            if (!(pq < 0.0e0)) goto S160;
            bound = 0.0e0;
            goto S170;
            S160:
            bound = 1.0e0;
            S170:
            eData.status = 3;
            return;
            S190:
            S180:
            if (which == 1) goto S220;
            //
            //  Select the minimum of P or Q
            //
            qporq = p <= q;
            if (!qporq) goto S200;
            porq = p;
            goto S210;
            S200:
            porq = q;
            S220:
            S210:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //  Calculating P and Q
                //
                eData.status = 0;
                cumchi(x, df, ref p, ref q);
                if (porq > 1.5e0)
                {
                    eData.status = 10;
                    return;
                }
            }
            else if (2 == which)
            {
                //
                //  Calculating X
                //
                x = 5.0e0;
                T3 = inf;
                T6 = atol;
                T7 = tol;
                eData.dstinv(K2, T3, K4, K4, K5, T6, T7);
                eData.status = 0;
                eData.dinvr();
                S230:
                if (!(eData.status == 1)) goto S270;
                cumchi(x, df, ref cum, ref ccum);
                if (!qporq) goto S240;
                eData.fx = cum - p;
                goto S250;
                S240:
                eData.fx = ccum - q;
                S250:
                if (!(eData.fx + porq > 1.5e0)) goto S260;
                eData.status = 10;
                return;
                S260:
                eData.dinvr();
                goto S230;
                S270:
                if (!(eData.status == -1)) goto S300;
                if (!eData.qleft) goto S280;
                eData.status = 1;
                bound = 0.0e0;
                goto S290;
                S280:
                eData.status = 2;
                bound = inf;
                S300:
                S290: ;
            }
            else if (3 == which)
            {
                //
                //  Calculating DF
                //
                df = 5.0e0;
                T8 = zero;
                T9 = inf;
                T10 = atol;
                T11 = tol;
                eData.dstinv(T8, T9, K4, K4, K5, T10, T11);
                eData.status = 0;
                eData.dinvr();
                S310:
                if (!(eData.status == 1)) goto S350;
                cumchi(x, df, ref cum, ref ccum);
                if (!qporq) goto S320;
                eData.fx = cum - p;
                goto S330;
                S320:
                eData.fx = ccum - q;
                S330:
                if (!(eData.fx + porq > 1.5e0)) goto S340;
                eData.status = 10;
                return;
                S340:
                eData.dinvr();
                goto S310;
                S350:
                if (!(eData.status == -1)) goto S380;
                if (!eData.qleft) goto S360;
                eData.status = 1;
                bound = zero;
                goto S370;
                S360:
                eData.status = 2;
                bound = inf;
                S370: ;
            }

            S380:
            return;
        }

        public static void cdfchn(int which, ref double p, ref double q, double x, ref double df,
                ref double pnonc, ref int status_, ref double bound)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CDFCHN evaluates the CDF of the Noncentral Chi-Square.
            //
            //  Discussion:
            //
            //    This routine calculates any one parameter of the noncentral chi-square
            //    distribution given values for the others.
            //
            //    The value P of the cumulative distribution function is calculated
            //    directly.
            //
            //    Computation of the other parameters involves a seach for a value that
            //    produces the desired value of P.  The search relies on the
            //    monotonicity of P with respect to the other parameters.
            //
            //    The computation time required for this routine is proportional
            //    to the noncentrality parameter (PNONC).  Very large values of
            //    this parameter can consume immense computer resources.  This is
            //    why the search range is bounded by 10,000.
            //
            //    The CDF of the noncentral chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
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
            //    1966, Formula 26.5.25.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input, int *WHICH, indicates which argument is to be calculated
            //    from the others.
            //    1: Calculate P and Q from X, DF and PNONC;
            //    2: Calculate X from P, DF and PNONC;
            //    3: Calculate DF from P, X and PNONC;
            //    4: Calculate PNONC from P, X and DF.
            //
            //    Input/output, double *P, the integral from 0 to X of
            //    the noncentral chi-square distribution.  If this is an input
            //    value, it should lie in the range: [0, 1.0-1.0D-16).
            //
            //    Input/output, double *Q, is generally not used by this
            //    subroutine and is only included for similarity with other routines.
            //    However, if P is to be computed, then a value will also be computed
            //    for Q.
            //
            //    Input, double *X, the upper limit of integration of the
            //    noncentral chi-square distribution.  If this is an input value, it
            //    should lie in the range: [0, +infinity).  If it is an output value,
            //    it will be sought in the range: [0,1.0D+300].
            //
            //    Input/output, double *DF, the number of degrees of freedom
            //    of the noncentral chi-square distribution.  If this is an input value,
            //    it should lie in the range: (0, +infinity).  If it is an output value,
            //    it will be searched for in the range: [ 1.0D-300, 1.0D+300].
            //
            //    Input/output, double *PNONC, the noncentrality parameter of
            //    the noncentral chi-square distribution.  If this is an input value, it
            //    should lie in the range: [0, +infinity).  If it is an output value,
            //    it will be searched for in the range: [0,1.0D+4]
            //
            //    Output, int *STATUS, reports on the calculation.
            //    0, if calculation completed correctly;
            //    -I, if input parameter number I is out of range;
            //    1, if the answer appears to be lower than the lowest search bound;
            //    2, if the answer appears to be higher than the greatest search bound.
            //
            //    Output, double *BOUND, is only defined if STATUS is nonzero.
            //    If STATUS is negative, then this is the value exceeded by parameter I.
            //    if STATUS is 1 or 2, this is the search bound that was exceeded.
            //
        {
            double tent4 = 1.0e4;
            double tol = (1.0e-8);
            double atol = (1.0e-50);
            double zero = (1.0e-300);
            double one = (1.0e0 - 1.0e-16);
            double inf = 1.0e300;

            double cum = 0;
            double ccum = 0;
            double K1 = 0.0e0;
            double K3 = 0.5e0;
            double K4 = 5.0e0;
            double T2;
            double T5;
            double T6;
            double T7;
            double T8;
            double T9;
            double T10;
            double T11;
            double T12;
            double T13;

            bound = 0.0;

            E0000E0001 eData = new E0000E0001();
            eData.status = 0;
            
            //
            //  Check arguments
            //
            if (!(which < 1 || which > 4)) goto S30;
            if (!(which < 1)) goto S10;
            bound = 1.0e0;
            goto S20;
            S10:
            bound = 4.0e0;
            S20:
            eData.status = -1;
            return;
            S30:
            if (which == 1) goto S70;
            //
            //  P
            //
            if (!(p < 0.0e0 || p > one)) goto S60;
            if (!(p < 0.0e0)) goto S40;
            bound = 0.0e0;
            goto S50;
            S40:
            bound = one;
            S50:
            eData.status = -2;
            return;
            S70:
            S60:
            if (which == 2) goto S90;
            //
            //     X
            //
            if (!(x < 0.0e0)) goto S80;
            bound = 0.0e0;
            eData.status = -4;
            return;
            S90:
            S80:
            if (which == 3) goto S110;
            //
            //     DF
            //
            if (!(df <= 0.0e0)) goto S100;
            bound = 0.0e0;
            eData.status = -5;
            return;
            S110:
            S100:
            if (which == 4) goto S130;
            //
            //  PNONC
            //
            if (!(pnonc < 0.0e0)) goto S120;
            bound = 0.0e0;
            eData.status = -6;
            return;
            S130:
            S120:
            //
            //     Calculate ANSWERS
            //
            if (1 == which)
            {
                //
                //     Calculating P and Q
                //
                cumchn(x, df, pnonc, ref p, ref q);
                eData.status = 0;
            }
            else if (2 == which)
            {
                //
                //     Calculating X
                //
                x = 5.0e0;
                T2 = inf;
                T5 = atol;
                T6 = tol;
                eData.dstinv(K1, T2, K3, K3, K4, T5, T6);
                eData.status = 0;
                eData.dinvr();
                S140:
                if (!(eData.status == 1)) goto S150;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                eData.fx = cum - p;
                eData.dinvr();
                goto S140;
                S150:
                if (!(eData.status == -1)) goto S180;
                if (!eData.qleft) goto S160;
                eData.status = 1;
                bound = 0.0e0;
                goto S170;
                S160:
                eData.status = 2;
                bound = inf;
                S180:
                S170: ;
            }
            else if (3 == which)
            {
                //
                //     Calculating DF
                //
                df = 5.0e0;
                T7 = zero;
                T8 = inf;
                T9 = atol;
                T10 = tol;
                eData.dstinv(T7, T8, K3, K3, K4, T9, T10);
                eData.status = 0;
                eData.dinvr();
                S190:
                if (!(eData.status == 1)) goto S200;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                eData.fx = cum - p;
                eData.dinvr();
                goto S190;
                S200:
                if (!(eData.status == -1)) goto S230;
                if (!eData.qleft) goto S210;
                eData.status = 1;
                bound = zero;
                goto S220;
                S210:
                eData.status = 2;
                bound = inf;
                S230:
                S220: ;
            }
            else if (4 == which)
            {
                //
                //     Calculating PNONC
                //
                pnonc = 5.0e0;
                T11 = tent4;
                T12 = atol;
                T13 = tol;
                eData.dstinv(K1, T11, K3, K3, K4, T12, T13);
                eData.status = 0;
                eData.dinvr();
                S240:
                if (!(eData.status == 1)) goto S250;
                cumchn(x, df, pnonc, ref cum, ref ccum);
                eData.fx = cum - p;
                eData.dinvr();
                goto S240;
                S250:
                if (!(eData.status == -1)) goto S280;
                if (!eData.qleft) goto S260;
                eData.status = 1;
                bound = zero;
                goto S270;
                S260:
                eData.status = 2;
                bound = tent4;
                S270: ;
            }

            S280:
            return;
        }

        public static void chi_noncentral_cdf_values(ref int n_data, ref double x, ref double lambda,
                ref int df, ref double cdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
            //
            //  Discussion:
            //
            //    The CDF of the noncentral chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
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
            //    Output, double *LAMBDA, the noncentrality parameter.
            //
            //    Output, int *DF, the number of degrees of freedom.
            //
            //    Output, double *CDF, the noncentral chi CDF.
            //
        {
            int N_MAX = 27;

            double[] cdf_vec =  {
                0.839944E+00, 0.695906E+00, 0.535088E+00,
                0.764784E+00, 0.620644E+00, 0.469167E+00,
                0.307088E+00, 0.220382E+00, 0.150025E+00,
                0.307116E-02, 0.176398E-02, 0.981679E-03,
                0.165175E-01, 0.202342E-03, 0.498448E-06,
                0.151325E-01, 0.209041E-02, 0.246502E-03,
                0.263684E-01, 0.185798E-01, 0.130574E-01,
                0.583804E-01, 0.424978E-01, 0.308214E-01,
                0.105788E+00, 0.794084E-01, 0.593201E-01
            }
            ;
            int[] df_vec =  {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                60, 80, 100,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10
            }
            ;
            double[] lambda_vec =  {
                0.5E+00, 0.5E+00, 0.5E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                5.0E+00, 5.0E+00, 5.0E+00,
                20.0E+00, 20.0E+00, 20.0E+00,
                30.0E+00, 30.0E+00, 30.0E+00,
                5.0E+00, 5.0E+00, 5.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00
            }
            ;
            double[] x_vec =  {
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                60.000E+00, 60.000E+00, 60.000E+00,
                0.050E+00, 0.050E+00, 0.050E+00,
                4.000E+00, 4.000E+00, 4.000E+00,
                5.000E+00, 5.000E+00, 5.000E+00,
                6.000E+00, 6.000E+00, 6.000E+00
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
                x = 0.0E+00;
                lambda = 0.0E+00;
                df = 0;
                cdf = 0.0E+00;
            }
            else
            {
                x = x_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                df = df_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
            }
        }

        public static void chi_square_cdf_values(ref int n_data, ref int a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
            //
            //  Discussion:
            //
            //    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
            //    commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[ChiSquareDistribution[DF], X ]
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int *A, the parameter of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 21;

            int[] a_vec =  {
                1, 2, 1, 2,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 3, 3, 3,
                3, 3, 10, 10,
                10
            }
            ;
            double[] fx_vec =  {
                0.0796557E+00, 0.00498752E+00, 0.112463E+00, 0.00995017E+00,
                0.472911E+00, 0.181269E+00, 0.0597575E+00, 0.0175231E+00,
                0.682689E+00, 0.393469E+00, 0.198748E+00, 0.090204E+00,
                0.0374342E+00, 0.427593E+00, 0.608375E+00, 0.738536E+00,
                0.828203E+00, 0.88839E+00, 0.000172116E+00, 0.00365985E+00,
                0.0185759E+00
            }
            ;
            double[] x_vec =  {
                0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
                0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
                5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
                3.00E+00
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

    }
}