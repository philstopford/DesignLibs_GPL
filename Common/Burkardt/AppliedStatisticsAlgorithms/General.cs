using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double alnfac(int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ALNFAC computes the logarithm of the factorial of N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the factorial.
            //
            //    Output, double ALNFAC, the logarithm of the factorial of N.
            //
        {
            return Helpers.LogGamma((double) (n + 1));
        }

        public static double alnorm(double x, bool upper)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ALNORM computes the cumulative density of the standard normal distribution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by David Hill.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    David Hill,
            //    Algorithm AS 66:
            //    The Normal Integral,
            //    Applied Statistics,
            //    Volume 22, Number 3, 1973, pages 424-427.
            //
            //  Parameters:
            //
            //    Input, double X, is one endpoint of the semi-infinite interval
            //    over which the integration takes place.
            //
            //    Input, bool UPPER, determines whether the upper or lower
            //    interval is to be integrated:
            //    .TRUE.  => integrate from X to + Infinity;
            //    .FALSE. => integrate from - Infinity to X.
            //
            //    Output, double ALNORM, the integral of the standard normal
            //    distribution over the desired interval.
            //
        {
            double a1 = 5.75885480458;
            double a2 = 2.62433121679;
            double a3 = 5.92885724438;
            double b1 = -29.8213557807;
            double b2 = 48.6959930692;
            double c1 = -0.000000038052;
            double c2 = 0.000398064794;
            double c3 = -0.151679116635;
            double c4 = 4.8385912808;
            double c5 = 0.742380924027;
            double c6 = 3.99019417011;
            double con = 1.28;
            double d1 = 1.00000615302;
            double d2 = 1.98615381364;
            double d3 = 5.29330324926;
            double d4 = -15.1508972451;
            double d5 = 30.789933034;
            double ltone = 7.0;
            double p = 0.398942280444;
            double q = 0.39990348504;
            double r = 0.398942280385;
            bool up;
            double utzero = 18.66;
            double value;
            double y;
            double z;

            up = upper;
            z = x;

            if (z < 0.0)
            {
                up = !up;
                z = -z;
            }

            if (ltone < z && ((!up) || utzero < z))
            {
                if (up)
                {
                    value = 0.0;
                }
                else
                {
                    value = 1.0;
                }

                return value;
            }

            y = 0.5 * z * z;

            if (z <= con)
            {
                value = 0.5 - z * (p - q * y
                    / (y + a1 + b1
                        / (y + a2 + b2
                            / (y + a3))));
            }
            else
            {
                value = r * Math.Exp(-y)
                        / (z + c1 + d1
                            / (z + c2 + d2
                                / (z + c3 + d3
                                    / (z + c4 + d4
                                        / (z + c5 + d5
                                            / (z + c6))))));
            }

            if (!up)
            {
                value = 1.0 - value;
            }

            return value;
        }

        public static double prncst(double st, int idf, double d, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PRNCST computes the lower tail of noncentral T distribution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by BE Cooper.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    BE Cooper,
            //    Algorithm AS 5:
            //    The Integral of the Non-Central T-Distribution,
            //    Applied Statistics,
            //    Volume 17, Number 2, 1968, page 193.
            //
            //  Parameters:
            //
            //    Input, double ST, the argument.
            //
            //    Input, int IDF, the number of degrees of freedom.
            //
            //    Input, double D, the noncentrality parameter.
            //
            //    Output, int *IFAULT, error flag.
            //    0, no error occurred.
            //    nonzero, an error occurred.
            //
            //    Output, double PRNCST, the value of the lower tail of
            //    the noncentral T distribution.
            //
            //  Local Parameters:
            //
            //    Local, double G1, 1.0 / sqrt(2.0 * pi)
            //
            //    Local, double G2, 1.0 / (2.0 * pi)
            //
            //    Local, double G3, sqrt(2.0 * pi)
            //
        {
            double a;
            double ak;
            double b;
            double da;
            double drb;
            double emin = 12.5;
            double f;
            double fk;
            double fkm1;
            double fmkm1;
            double fmkm2;
            double g1 = 0.3989422804;
            double g2 = 0.1591549431;
            double g3 = 2.5066282746;
            int ioe;
            int k;
            double rb;
            double sum;
            double value;

            f = (double) (idf);
            //
            //  For very large IDF, use the normal approximation.
            //
            if (100 < idf)
            {
                ifault = 1;

                a = Math.Sqrt(0.5 * f) * Math.Exp(Helpers.LogGamma(0.5 * (f - 1.0))
                                                  - Helpers.LogGamma(0.5 * f)) * d;

                value = alnorm((st - a) / Math.Sqrt(f * (1.0 + d * d)
                    / (f - 2.0) - a * a), false);
                return value;
            }

            ifault = 0;
            ioe = (idf % 2);
            a = st / Math.Sqrt(f);
            b = f / (f + st * st);
            rb = Math.Sqrt(b);
            da = d * a;
            drb = d * rb;

            if (idf == 1)
            {
                value = alnorm(drb, true) + 2.0 * tfn(drb, a);
                return value;
            }

            sum = 0.0;

            if (Math.Abs(drb) < emin)
            {
                fmkm2 = a * rb * Math.Exp(-0.5 * drb * drb)
                        * alnorm(a * drb, false) * g1;
            }
            else
            {
                fmkm2 = 0.0;
            }

            fmkm1 = b * da * fmkm2;

            if (Math.Abs(d) < emin)
            {
                fmkm1 = fmkm1 + b * a * g2 * Math.Exp(-0.5 * d * d);
            }

            if (ioe == 0)
            {
                sum = fmkm2;
            }
            else
            {
                sum = fmkm1;
            }

            ak = 1.0;
            fk = 2.0;

            for (k = 2; k <= idf - 2; k = k + 2)
            {
                fkm1 = fk - 1.0;
                fmkm2 = b * (da * ak * fmkm1 + fmkm2) * fkm1 / fk;
                ak = 1.0 / (ak * fkm1);
                fmkm1 = b * (da * ak * fmkm2 + fmkm1) * fk / (fk + 1.0);

                if (ioe == 0)
                {
                    sum = sum + fmkm2;
                }
                else
                {
                    sum = sum + fmkm1;
                }

                ak = 1.0 / (ak * fk);
                fk = fk + 2.0;
            }

            if (ioe == 0)
            {
                value = alnorm(d, true) + sum * g3;
            }
            else
            {
                value = alnorm(drb, true) + 2.0 * (sum + tfn(drb, a));
            }

            return value;
        }


        public static double tfn(double x, double fx)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TFN calculates the T-function of Owen.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by JC Young, Christoph Minder.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    MA Porter, DJ Winstanley,
            //    Remark AS R30:
            //    A Remark on Algorithm AS76:
            //    An Integral Useful in Calculating Noncentral T and Bivariate
            //    Normal Probabilities,
            //    Applied Statistics,
            //    Volume 28, Number 1, 1979, page 113.
            //
            //    JC Young, Christoph Minder,
            //    Algorithm AS 76: 
            //    An Algorithm Useful in Calculating Non-Central T and 
            //    Bivariate Normal Distributions,
            //    Applied Statistics,
            //    Volume 23, Number 3, 1974, pages 455-457.
            //
            //  Parameters:
            //
            //    Input, double X, FX, the parameters of the function.
            //
            //    Output, double TFN, the value of the T-function.
            //
        {
            int NG = 5;

            double fxs;
            int i;
            double[] r =
            {
                0.1477621,
                0.1346334,
                0.1095432,
                0.0747257,
                0.0333357
            };
            double r1;
            double r2;
            double rt;
            double tp = 0.159155;
            double tv1 = 1.0E-35;
            double tv2 = 15.0;
            double tv3 = 15.0;
            double tv4 = 1.0E-05;
            double[] u =
            {
                0.0744372,
                0.2166977,
                0.3397048,
                0.4325317,
                0.4869533
            };
            double value;
            double x1;
            double x2;
            double xs;
            //
            //  Test for X near zero.
            //
            if (Math.Abs(x) < tv1)
            {
                value = tp * Math.Atan(fx);
                return value;
            }

            //
            //  Test for large values of abs(X).
            //
            if (tv2 < Math.Abs(x))
            {
                value = 0.0;
                return value;
            }

            //
            //  Test for FX near zero.
            //
            if (Math.Abs(fx) < tv1)
            {
                value = 0.0;
                return value;
            }

            //
            //  Test whether abs ( FX ) is so large that it must be truncated.
            //
            xs = -0.5 * x * x;
            x2 = fx;
            fxs = fx * fx;
            //
            //  Computation of truncation point by Newton iteration.
            //
            if (tv3 <= Math.Log(1.0 + fxs) - xs * fxs)
            {
                x1 = 0.5 * fx;
                fxs = 0.25 * fxs;

                for (;;)
                {
                    rt = fxs + 1.0;

                    x2 = x1 + (xs * fxs + tv3 - Math.Log(rt))
                        / (2.0 * x1 * (1.0 / rt - xs));

                    fxs = x2 * x2;

                    if (Math.Abs(x2 - x1) < tv4)
                    {
                        break;
                    }

                    x1 = x2;
                }
            }

            //
            //  Gaussian quadrature.
            //
            rt = 0.0;
            for (i = 0; i < NG; i++)
            {
                r1 = 1.0 + fxs * Math.Pow(0.5 + u[i], 2);
                r2 = 1.0 + fxs * Math.Pow(0.5 - u[i], 2);

                rt = rt + r[i] * (Math.Exp(xs * r1) / r1 + Math.Exp(xs * r2) / r2);
            }

            value = rt * x2 * tp;

            return value;
        }

        public static void rnorm(ref int seed, ref double u1, ref double u2)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RNORM returns two independent standard random normal deviates.
            //
            //  Discussion:
            //
            //    This routine sets U1 and U2 to two independent standardized 
            //    random normal deviates.   This is a version of the 
            //    method given in Knuth.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 April 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by William Smith, Ronald Hocking.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Donald Knuth,
            //    The Art of Computer Programming,
            //    Volume 2, Seminumerical Algorithms,
            //    Third Edition,
            //    Addison Wesley, 1997,
            //    ISBN: 0201896842,
            //    LC: QA76.6.K64.
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double &U1, &U2, two standard random normal deviates.
            //
        {
            for (;;)
            {
                double x = UniformRNG.r8_uniform_01(ref seed);
                double y = UniformRNG.r8_uniform_01(ref seed);
                x = 2.0 * x - 1.0;
                y = 2.0 * y - 1.0;
                double s = x * x + y * y;

                if (s <= 1.0)
                {
                    s = Math.Sqrt(-2.0 * Math.Log(s) / s);
                    u1 = x * s;
                    u2 = y * s;
                    break;
                }
            }
        }

        public static double ppnd(double p, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PPND produces the normal deviate value corresponding to lower tail area = P.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by J Beasley, S Springer.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    J Beasley, S Springer,
            //    Algorithm AS 111:
            //    The Percentage Points of the Normal Distribution,
            //    Applied Statistics,
            //    Volume 26, Number 1, 1977, pages 118-121.
            //
            //  Parameters:
            //
            //    Input, double P, the value of the cumulative probability
            //    densitity function.  0 < P < 1.
            //
            //    Output, integer *IFAULT, error flag.
            //    0, no error.
            //    1, P <= 0 or P >= 1.  PPND is returned as 0.
            //
            //    Output, double PPND, the normal deviate value with the property that
            //    the probability of a standard normal deviate being less than or
            //    equal to PPND is P.
            //
        {
            double a0 = 2.50662823884;
            double a1 = -18.61500062529;
            double a2 = 41.39119773534;
            double a3 = -25.44106049637;
            double b1 = -8.47351093090;
            double b2 = 23.08336743743;
            double b3 = -21.06224101826;
            double b4 = 3.13082909833;
            double c0 = -2.78718931138;
            double c1 = -2.29796479134;
            double c2 = 4.85014127135;
            double c3 = 2.32121276858;
            double d1 = 3.54388924762;
            double d2 = 1.63706781897;
            double r;
            double split = 0.42;
            double value;

            ifault = 0;
            //
            //  0.08 < P < 0.92
            //
            if (Math.Abs(p - 0.5) <= split)
            {
                r = (p - 0.5) * (p - 0.5);

                value = (p - 0.5) * (((
                                          a3 * r
                                          + a2) * r
                                      + a1) * r
                                     + a0) / ((((
                                                    b4 * r
                                                    + b3) * r
                                                + b2) * r
                                               + b1) * r
                                              + 1.0);
            }
            //
            //  P < 0.08 or P > 0.92,
            //  R = min ( P, 1-P )
            //
            else if (0.0 < p && p < 1.0)
            {
                if (0.5 < p)
                {
                    r = Math.Sqrt(-Math.Log(1.0 - p));
                }
                else
                {
                    r = Math.Sqrt(-Math.Log(p));
                }

                value = (((
                              c3 * r
                              + c2) * r
                          + c1) * r
                         + c0) / ((
                                      d2 * r
                                      + d1) * r
                                  + 1.0);

                if (p < 0.5)
                {
                    value = -value;
                }
            }
            //
            //  P <= 0.0 or 1.0 <= P
            //
            else
            {
                ifault = 1;
                value = 0.0;
            }

            return value;
        }

        public static void normal_01_cdf_values(ref int n_data, ref double x, ref double fx)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NormalDistribution [ 0, 1 ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 August 2004
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
            int N_MAX = 17;

            double[] fx_vec =
                {
                    0.5000000000000000E+00,
                    0.5398278372770290E+00,
                    0.5792597094391030E+00,
                    0.6179114221889526E+00,
                    0.6554217416103242E+00,
                    0.6914624612740131E+00,
                    0.7257468822499270E+00,
                    0.7580363477769270E+00,
                    0.7881446014166033E+00,
                    0.8159398746532405E+00,
                    0.8413447460685429E+00,
                    0.9331927987311419E+00,
                    0.9772498680518208E+00,
                    0.9937903346742239E+00,
                    0.9986501019683699E+00,
                    0.9997673709209645E+00,
                    0.9999683287581669E+00
                }
                ;

            double[] x_vec =
                {
                    0.0000000000000000E+00,
                    0.1000000000000000E+00,
                    0.2000000000000000E+00,
                    0.3000000000000000E+00,
                    0.4000000000000000E+00,
                    0.5000000000000000E+00,
                    0.6000000000000000E+00,
                    0.7000000000000000E+00,
                    0.8000000000000000E+00,
                    0.9000000000000000E+00,
                    0.1000000000000000E+01,
                    0.1500000000000000E+01,
                    0.2000000000000000E+01,
                    0.2500000000000000E+01,
                    0.3000000000000000E+01,
                    0.3500000000000000E+01,
                    0.4000000000000000E+01
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }

        }

        public static double betain(double x, double p, double q, double beta, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETAIN computes the incomplete Beta function ratio.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    KL Majumder, GP Bhattacharjee,
            //    Algorithm AS 63:
            //    The incomplete Beta Integral,
            //    Applied Statistics,
            //    Volume 22, Number 3, 1973, pages 409-411.
            //
            //  Parameters:
            //
            //    Input, double X, the argument, between 0 and 1.
            //
            //    Input, double P, Q, the parameters, which
            //    must be positive.
            //
            //    Input, double BETA, the logarithm of the complete
            //    beta function.
            //
            //    Output, int &IFAULT, error flag.
            //    0, no error.
            //    nonzero, an error occurred.
            //
            //    Output, double BETAIN, the value of the incomplete
            //    Beta function ratio.
            //
        {
            double acu = 0.1E-14;
            bool indx;
            double pp;
            double qq;
            double xx;

            double value = x;
            ifault = 0;
            //
            //  Check the input arguments.
            //
            if (p <= 0.0 || q <= 0.0)
            {
                ifault = 1;
                return value;
            }

            if (x < 0.0 || 1.0 < x)
            {
                ifault = 2;
                return value;
            }

            //
            //  Special cases.
            //
            if (x == 0.0 || x == 1.0)
            {
                return value;
            }

            //
            //  Change tail if necessary and determine S.
            //
            double psq = p + q;
            double cx = 1.0 - x;

            if (p < psq * x)
            {
                xx = cx;
                cx = x;
                pp = q;
                qq = p;
                indx = true;
            }
            else
            {
                xx = x;
                pp = p;
                qq = q;
                indx = false;
            }

            double term = 1.0;
            double ai = 1.0;
            value = 1.0;
            int ns = (int) (qq + cx * psq);
            //
            //  Use the Soper reduction formula.
            //
            double rx = xx / cx;
            double temp = qq - ai;
            if (ns == 0)
            {
                rx = xx;
            }

            for (;;)
            {
                term = term * temp * rx / (pp + ai);
                value = value + term;
                ;
                temp = Math.Abs(term);

                if (temp <= acu && temp <= acu * value)
                {
                    value = value * Math.Exp(pp * Math.Log(xx)
                        + (qq - 1.0) * Math.Log(cx) - beta) / pp;

                    if (indx)
                    {
                        value = 1.0 - value;
                    }

                    break;
                }

                ai = ai + 1.0;
                ns = ns - 1;

                if (0 <= ns)
                {
                    temp = qq - ai;
                    if (ns == 0)
                    {
                        rx = xx;
                    }
                }
                else
                {
                    temp = psq;
                    psq = psq + 1.0;
                }
            }

            return value;
        }

        public static double gammad(double x, double p, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMAD computes the Incomplete Gamma Integral
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by B Shea.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    B Shea,
            //    Algorithm AS 239:
            //    Chi-squared and Incomplete Gamma Integral,
            //    Applied Statistics,
            //    Volume 37, Number 3, 1988, pages 466-473.
            //
            //  Parameters:
            //
            //    Input, double X, P, the parameters of the incomplete 
            //    gamma ratio.  0 <= X, and 0 < P.
            //
            //    Output, int IFAULT, error flag.
            //    0, no error.
            //    1, X < 0 or P <= 0.
            //
            //    Output, double GAMMAD, the value of the incomplete 
            //    Gamma integral.
            //
        {
            double a;
            double arg;
            double c;
            double elimit = -88.0;
            double oflo = 1.0E+37;
            double plimit = 1000.0;
            double pn1;
            double tol = 1.0E-14;
            double xbig = 1.0E+08;

            double value = 0.0;
            //
            //  Check the input.
            //
            if (x < 0.0)
            {
                ifault = 1;
                return value;
            }

            if (p <= 0.0)
            {
                ifault = 1;
                return value;
            }

            ifault = 0;

            if (x == 0.0)
            {
                value = 0.0;
                return value;
            }

            //
            //  If P is large, use a normal approximation.
            //
            if (plimit < p)
            {
                pn1 = 3.0 * Math.Sqrt(p) * (Math.Pow(x / p, 1.0 / 3.0)
                    + 1.0 / (9.0 * p) - 1.0);

                bool upper = false;
                value = alnorm(pn1, upper);
                return value;
            }

            //
            //  If X is large set value = 1.
            //
            if (xbig < x)
            {
                value = 1.0;
                return value;
            }

            //
            //  Use Pearson's series expansion.
            //  (Note that P is not large enough to force overflow in ALOGAM).
            //  No need to test IFAULT on exit since P > 0.
            //
            if (x <= 1.0 || x < p)
            {
                arg = p * Math.Log(x) - x - Helpers.LogGamma(p + 1.0);
                c = 1.0;
                value = 1.0;
                a = p;

                for (;;)
                {
                    a = a + 1.0;
                    c = c * x / a;
                    value = value + c;

                    if (c <= tol)
                    {
                        break;
                    }
                }

                arg = arg + Math.Log(value);

                if (elimit <= arg)
                {
                    value = Math.Exp(arg);
                }
                else
                {
                    value = 0.0;
                }
            }
            //
            //  Use a continued fraction expansion.
            //
            else
            {
                arg = p * Math.Log(x) - x - Helpers.LogGamma(p);
                a = 1.0 - p;
                double b = a + x + 1.0;
                c = 0.0;
                pn1 = 1.0;
                double pn2 = x;
                double pn3 = x + 1.0;
                double pn4 = x * b;
                value = pn3 / pn4;

                for (;;)
                {
                    a = a + 1.0;
                    b = b + 2.0;
                    c = c + 1.0;
                    double an = a * c;
                    double pn5 = b * pn3 - an * pn1;
                    double pn6 = b * pn4 - an * pn2;

                    if (pn6 != 0.0)
                    {
                        double rn = pn5 / pn6;

                        if (Math.Abs(value - rn) <= Math.Min(tol, tol * rn))
                        {
                            break;
                        }

                        value = rn;
                    }

                    pn1 = pn3;
                    pn2 = pn4;
                    pn3 = pn5;
                    pn4 = pn6;
                    //
                    //  Re-scale terms in continued fraction if terms are large.
                    //
                    if (oflo <= Math.Abs(pn5))
                    {
                        pn1 = pn1 / oflo;
                        pn2 = pn2 / oflo;
                        pn3 = pn3 / oflo;
                        pn4 = pn4 / oflo;
                    }
                }

                arg = arg + Math.Log(value);

                if (elimit <= arg)
                {
                    value = 1.0 - Math.Exp(arg);
                }
                else
                {
                    value = 1.0;
                }
            }

            return value;
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
            //      PN(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T**(A-1) * exp(-T) dT.
            //
            //    With this definition, for all A and X,
            //
            //      0 <= PN(A,X) <= 1
            //
            //    and
            //
            //      PN(A,INFINITY) = 1.0
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      1 - GammaRegularized[A,X]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 November 2004
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
            int N_MAX = 20;

            double[] a_vec =
                {
                    0.10E+00,
                    0.10E+00,
                    0.10E+00,
                    0.50E+00,
                    0.50E+00,
                    0.50E+00,
                    0.10E+01,
                    0.10E+01,
                    0.10E+01,
                    0.11E+01,
                    0.11E+01,
                    0.11E+01,
                    0.20E+01,
                    0.20E+01,
                    0.20E+01,
                    0.60E+01,
                    0.60E+01,
                    0.11E+02,
                    0.26E+02,
                    0.41E+02
                }
                ;

            double[] fx_vec =
                {
                    0.7382350532339351E+00,
                    0.9083579897300343E+00,
                    0.9886559833621947E+00,
                    0.3014646416966613E+00,
                    0.7793286380801532E+00,
                    0.9918490284064973E+00,
                    0.9516258196404043E-01,
                    0.6321205588285577E+00,
                    0.9932620530009145E+00,
                    0.7205974576054322E-01,
                    0.5891809618706485E+00,
                    0.9915368159845525E+00,
                    0.1018582711118352E-01,
                    0.4421745996289254E+00,
                    0.9927049442755639E+00,
                    0.4202103819530612E-01,
                    0.9796589705830716E+00,
                    0.9226039842296429E+00,
                    0.4470785799755852E+00,
                    0.7444549220718699E+00
                }
                ;

            double[] x_vec =
                {
                    0.30E-01,
                    0.30E+00,
                    0.15E+01,
                    0.75E-01,
                    0.75E+00,
                    0.35E+01,
                    0.10E+00,
                    0.10E+01,
                    0.50E+01,
                    0.10E+00,
                    0.10E+01,
                    0.50E+01,
                    0.15E+00,
                    0.15E+01,
                    0.70E+01,
                    0.25E+01,
                    0.12E+02,
                    0.16E+02,
                    0.25E+02,
                    0.45E+02
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
                a = 0.0;
                x = 0.0;
                fx = 0.0;
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