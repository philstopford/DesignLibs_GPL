using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static double alngam(double xvalue, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ALNGAM computes the logarithm of the gamma function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Allan Macleod.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Allan Macleod,
            //    Algorithm AS 245,
            //    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
            //    Applied Statistics,
            //    Volume 38, Number 2, 1989, pages 397-402.
            //
            //  Parameters:
            //
            //    Input, double XVALUE, the argument of the Gamma function.
            //
            //    Output, int &IFAULT, error flag.
            //    0, no error occurred.
            //    1, XVALUE is less than or equal to 0.
            //    2, XVALUE is too big.
            //
            //    Output, double ALNGAM, the logarithm of the gamma function of X.
            //
        {
            double alr2pi = 0.918938533204673;
            double[] r1 =
                {
                    -2.66685511495,
                    -24.4387534237,
                    -21.9698958928,
                    11.1667541262,
                    3.13060547623,
                    0.607771387771,
                    11.9400905721,
                    31.4690115749,
                    15.2346874070
                }
                ;
            double[] r2 =
                {
                    -78.3359299449,
                    -142.046296688,
                    137.519416416,
                    78.6994924154,
                    4.16438922228,
                    47.0668766060,
                    313.399215894,
                    263.505074721,
                    43.3400022514
                }
                ;
            double[] r3 =
                {
                    -2.12159572323E+05,
                    2.30661510616E+05,
                    2.74647644705E+04,
                    -4.02621119975E+04,
                    -2.29660729780E+03,
                    -1.16328495004E+05,
                    -1.46025937511E+05,
                    -2.42357409629E+04,
                    -5.70691009324E+02
                }
                ;
            double[] r4 =
                {
                    0.279195317918525,
                    0.4917317610505968,
                    0.0692910599291889,
                    3.350343815022304,
                    6.012459259764103
                }
                ;
            double value;
            double x;
            double x1;
            double x2;
            double xlge = 510000.0;
            double xlgst = 1.0E+30;
            double y;

            x = xvalue;
            value = 0.0;
            //
            //  Check the input.
            //
            if (xlgst <= x)
            {
                ifault = 2;
                return value;
            }

            if (x <= 0.0)
            {
                ifault = 1;
                return value;
            }

            ifault = 0;
            //
            //  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
            //
            if (x < 1.5)
            {
                if (x < 0.5)
                {
                    value = -Math.Log(x);
                    y = x + 1.0;
                    //
                    //  Test whether X < machine epsilon.
                    //
                    if (y == 1.0)
                    {
                        return value;
                    }
                }
                else
                {
                    value = 0.0;
                    y = x;
                    x = (x - 0.5) - 0.5;
                }

                value = value + x * ((((
                                           r1[4] * y
                                           + r1[3]) * y
                                       + r1[2]) * y
                                      + r1[1]) * y
                                     + r1[0]) / ((((
                                                       y
                                                       + r1[8]) * y
                                                   + r1[7]) * y
                                                  + r1[6]) * y
                                                 + r1[5]);

                return value;
            }

            //
            //  Calculation for 1.5 <= X < 4.0.
            //
            if (x < 4.0)
            {
                y = (x - 1.0) - 1.0;

                value = y * ((((
                                   r2[4] * x
                                   + r2[3]) * x
                               + r2[2]) * x
                              + r2[1]) * x
                             + r2[0]) / ((((
                                               x
                                               + r2[8]) * x
                                           + r2[7]) * x
                                          + r2[6]) * x
                                         + r2[5]);
            }
            //
            //  Calculation for 4.0 <= X < 12.0.
            //
            else if (x < 12.0)
            {
                value = ((((
                               r3[4] * x
                               + r3[3]) * x
                           + r3[2]) * x
                          + r3[1]) * x
                         + r3[0]) / ((((
                                           x
                                           + r3[8]) * x
                                       + r3[7]) * x
                                      + r3[6]) * x
                                     + r3[5]);
            }
            //
            //  Calculation for 12.0 <= X.
            //
            else
            {
                y = Math.Log(x);
                value = x * (y - 1.0) - 0.5 * y + alr2pi;

                if (x <= xlge)
                {
                    x1 = 1.0 / x;
                    x2 = x1 * x1;

                    value = value + x1 * ((
                            r4[2] *
                            x2 + r4[1]) *
                        x2 + r4[0]) / ((
                            x2 + r4[4]) *
                        x2 + r4[3]);
                }
            }

            return value;
        }

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
            //    C++ version by John Burkardt.
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


        public static double gamain(double x, double p, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMAIN computes the incomplete gamma ratio.
            //
            //  Discussion:
            //
            //    A series expansion is used if P > X or X <= 1.  Otherwise, a
            //    continued fraction approximation is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 June 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by G Bhattacharjee.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    G Bhattacharjee,
            //    Algorithm AS 32:
            //    The Incomplete Gamma Integral,
            //    Applied Statistics,
            //    Volume 19, Number 3, 1970, pages 285-287.
            //
            //  Parameters:
            //
            //    Input, double X, P, the parameters of the incomplete 
            //    gamma ratio.  0 <= X, and 0 < P.
            //
            //    Output, int *IFAULT, error flag.
            //    0, no errors.
            //    1, P <= 0.
            //    2, X < 0.
            //    3, underflow.
            //    4, error return from the Log Gamma routine.
            //
            //    Output, double GAMAIN, the value of the incomplete gamma ratio.
            //
        {
            double acu = 1.0E-08;
            double gin;
            double oflo = 1.0E+37;
            double[] pn = new double[6];
            double rn;
            double term;
            double uflo = 1.0E-37;
            double value;

            ifault = 0;
            //
            //  Check the input.
            //
            if (p <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            if (x < 0.0)
            {
                ifault = 2;
                value = 0.0;
                return value;
            }

            if (x == 0.0)
            {
                ifault = 0;
                value = 0.0;
                return value;
            }

            double g = Helpers.LogGamma(p);

            double arg = p * Math.Log(x) - x - g;

            if (arg < Math.Log(uflo))
            {
                ifault = 3;
                value = 0.0;
                return value;
            }

            ifault = 0;
            double factor = Math.Exp(arg);
            //
            //  Calculation by series expansion.
            //
            if (x <= 1.0 || x < p)
            {
                gin = 1.0;
                term = 1.0;
                rn = p;

                for (;;)
                {
                    rn = rn + 1.0;
                    term = term * x / rn;
                    gin = gin + term;

                    if (term <= acu)
                    {
                        break;
                    }
                }

                value = gin * factor / p;
                return value;
            }

            //
            //  Calculation by continued fraction.
            //
            double a = 1.0 - p;
            double b = a + x + 1.0;
            term = 0.0;

            pn[0] = 1.0;
            pn[1] = x;
            pn[2] = x + 1.0;
            pn[3] = x * b;

            gin = pn[2] / pn[3];

            for (;;)
            {
                a = a + 1.0;
                b = b + 2.0;
                term = term + 1.0;
                double an = a * term;
                for (int i = 0; i <= 1; i++)
                {
                    pn[i + 4] = b * pn[i + 2] - an * pn[i];
                }

                if (pn[5] != 0.0)
                {
                    rn = pn[4] / pn[5];
                    double dif = Math.Abs(gin - rn);
                    //
                    //  Absolute error tolerance satisfied?
                    //
                    if (dif <= acu)
                    {
                        //
                        //  Relative error tolerance satisfied?
                        //
                        if (dif <= acu * rn)
                        {
                            value = 1.0 - factor * gin;
                            break;
                        }
                    }

                    gin = rn;
                }

                for (int i = 0; i < 4; i++)
                {
                    pn[i] = pn[i + 2];
                }

                if (oflo <= Math.Abs(pn[4]))
                {
                    for (int i = 0; i < 4; i++)
                    {
                        pn[i] = pn[i] / oflo;
                    }
                }
            }

            return value;
        }

        public static double gammds(double x, double p, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMDS computes the incomplete Gamma integral.
            //
            //  Discussion:
            //
            //    The parameters must be positive.  An infinite series is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Chi Leung Lau.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Chi Leung Lau,
            //    Algorithm AS 147:
            //    A Simple Series for the Incomplete Gamma Integral,
            //    Applied Statistics,
            //    Volume 29, Number 1, 1980, pages 113-114.
            //
            //  Parameters:
            //
            //    Input, double X, P, the arguments of the incomplete
            //    Gamma integral.  X and P must be greater than 0.
            //
            //    Output, int *IFAULT, error flag.
            //    0, no errors.
            //    1, X <= 0 or P <= 0.
            //    2, underflow during the computation.
            //
            //    Output, double GAMMDS, the value of the incomplete
            //    Gamma integral.
            //
        {
            double e = 1.0E-09;
            double uflo = 1.0E-37;
            double value;
            //
            //  Check the input.
            //
            if (x <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            if (p <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            //
            //  LGAMMA is the natural logarithm of the gamma function.
            //
            double arg = p * Math.Log(x) - Helpers.LogGamma(p + 1.0) - x;

            if (arg < Math.Log(uflo))
            {
                value = 0.0;
                ifault = 2;
                return value;
            }

            double f = Math.Exp(arg);

            if (f == 0.0)
            {
                value = 0.0;
                ifault = 2;
                return value;
            }

            ifault = 0;
            //
            //  Series begins.
            //
            double c = 1.0;
            value = 1.0;
            double a = p;

            for (;;)
            {
                a = a + 1.0;
                c = c * x / a;
                value = value + c;

                if (c <= e * value)
                {
                    break;
                }
            }

            value = value * f;

            return value;
        }

        public static double lngamma(double z, ref int ier)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LNGAMMA computes Log(Gamma(X)) using a Lanczos approximation.
            //
            //  Discussion:
            //
            //    This algorithm is not part of the Applied Statistics algorithms.
            //    It is slower but gives 14 or more significant decimal digits
            //    accuracy, except around X = 1 and X = 2.   The Lanczos series from
            //    which this algorithm is derived is interesting in that it is a
            //    convergent series approximation for the gamma function, whereas
            //    the familiar series due to De Moivre (and usually wrongly called
            //    the Stirling approximation) is only an asymptotic approximation, as
            //    is the true and preferable approximation due to Stirling.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Alan Miller.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Cornelius Lanczos,
            //    A precision approximation of the gamma function,
            //    SIAM Journal on Numerical Analysis, B,
            //    Volume 1, 1964, pages 86-96.
            //
            //  Parameters:
            //
            //    Input, double Z, the argument of the Gamma function.
            //
            //    Output, int &IER, error flag.
            //    0, no error occurred.
            //    1, Z is less than or equal to 0.
            //
            //    Output, double LNGAMMA, the logarithm of the gamma function of Z.
            //
        {
            double[] a =
                {
                    0.9999999999995183,
                    676.5203681218835,
                    -1259.139216722289,
                    771.3234287757674,
                    -176.6150291498386,
                    12.50734324009056,
                    -0.1385710331296526,
                    0.9934937113930748E-05,
                    0.1659470187408462E-06
                }
                ;
            int j;
            double lnsqrt2pi = 0.9189385332046727;
            double tmp;
            double value;

            if (z <= 0.0)
            {
                ier = 1;
                value = 0.0;
                return value;
            }

            ier = 0;

            value = 0.0;
            tmp = z + 7.0;
            for (j = 8; 1 <= j; j--)
            {
                value = value + a[j] / tmp;
                tmp = tmp - 1.0;
            }

            value = value + a[0];
            value = Math.Log(value) + lnsqrt2pi - (z + 6.5)
                    + (z - 0.5) * Math.Log(z + 6.5);

            return value;
        }

        public static void normp(double z, ref double p, ref double q, ref double pdf)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMP computes the cumulative density of the standard normal distribution.
            //
            //  Discussion:
            //
            //    This is algorithm 5666 from Hart, et al.
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
            //    Original FORTRAN77 version by Alan Miller.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
            //    Charles Mesztenyi, John Rice, Henry Thacher, 
            //    Christoph Witzgall,
            //    Computer Approximations,
            //    Wiley, 1968,
            //    LC: QA297.C64.
            //
            //  Parameters:
            //
            //    Input, double Z, divides the real line into two 
            //    semi-infinite intervals, over each of which the standard normal 
            //    distribution is to be integrated.
            //
            //    Output, double *P, *Q, the integrals of the standard normal
            //    distribution over the intervals ( - Infinity, Z] and 
            //    [Z, + Infinity ), respectively.
            //
            //    Output, double *PDF, the value of the standard normal distribution
            //    at Z.
            //
        {
            double cutoff = 7.071;
            double p0 = 220.2068679123761;
            double p1 = 221.2135961699311;
            double p2 = 112.0792914978709;
            double p3 = 33.91286607838300;
            double p4 = 6.373962203531650;
            double p5 = 0.7003830644436881;
            double p6 = 0.03526249659989109;
            double q0 = 440.4137358247522;
            double q1 = 793.8265125199484;
            double q2 = 637.3336333788311;
            double q3 = 296.5642487796737;
            double q4 = 86.78073220294608;
            double q5 = 16.06417757920695;
            double q6 = 1.755667163182642;
            double q7 = 0.08838834764831844;
            double root2pi = 2.506628274631001;

            double zabs = Math.Abs(z);
            //
            //  37 < |Z|.
            //
            if (37.0 < zabs)
            {
                pdf = 0.0;
                p = 0.0;
            }
            //
            //  |Z| <= 37.
            //
            else
            {
                double expntl = Math.Exp(-0.5 * zabs * zabs);
                pdf = expntl / root2pi;
                //
                //  |Z| < CUTOFF = 10 / sqrt(2).
                //
                if (zabs < cutoff)
                {
                    p = expntl * ((((((
                                          p6 * zabs
                                          + p5) * zabs
                                      + p4) * zabs
                                     + p3) * zabs
                                    + p2) * zabs
                                   + p1) * zabs
                                  + p0) / (((((((
                                                    q7 * zabs
                                                    + q6) * zabs
                                                + q5) * zabs
                                               + q4) * zabs
                                              + q3) * zabs
                                             + q2) * zabs
                                            + q1) * zabs
                                           + q0);
                }
                //
                //  CUTOFF <= |Z|.
                //
                else
                {
                    p = pdf / (
                        zabs + 1.0 / (
                            zabs + 2.0 / (
                                zabs + 3.0 / (
                                    zabs + 4.0 / (
                                        zabs + 0.65)))));
                }
            }

            if (z < 0.0)
            {
                q = 1.0 - p;
            }
            else
            {
                q = p;
                p = 1.0 - q;
            }
        }

        public static void nprob(double z, ref double p, ref double q, ref double pdf)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPROB computes the cumulative density of the standard normal distribution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2013
            //
            //  Author:
            //
            //    Original FORTRAN77 version by AG Adams.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    AG Adams,
            //    Algorithm 39:
            //    Areas Under the Normal Curve,
            //    Computer Journal,
            //    Volume 12, Number 2, May 1969, pages 197-198.
            //
            //  Parameters:
            //
            //    Input, double Z, divides the real line into 
            //    two semi-infinite intervals, over each of which the standard normal 
            //    distribution is to be integrated.
            //
            //    Output, double *P, *Q, the integrals of the standard normal
            //    distribution over the intervals ( - Infinity, Z] and 
            //    [Z, + Infinity ), respectively.
            //
            //    Output, double *PDF, the value of the standard normal
            //    distribution at Z.
            //
        {
            double a0 = 0.5;
            double a1 = 0.398942280444;
            double a2 = 0.399903438504;
            double a3 = 5.75885480458;
            double a4 = 29.8213557808;
            double a5 = 2.62433121679;
            double a6 = 48.6959930692;
            double a7 = 5.92885724438;
            double b0 = 0.398942280385;
            double b1 = 0.000000038052;
            double b2 = 1.00000615302;
            double b3 = 0.000398064794;
            double b4 = 1.98615381364;
            double b5 = 0.151679116635;
            double b6 = 5.29330324926;
            double b7 = 4.8385912808;
            double b8 = 15.1508972451;
            double b9 = 0.742380924027;
            double b10 = 30.789933034;
            double b11 = 3.99019417011;
            double y;
            double zabs;

            zabs = Math.Abs(z);
            //
            //  |Z| between 0 and 1.28
            //
            if (zabs <= 1.28)
            {
                y = a0 * z * z;
                pdf = Math.Exp(-y) * b0;

                q = a0 - zabs * (a1 - a2 * y
                    / (y + a3 - a4
                        / (y + a5 + a6
                            / (y + a7))));
            }
            //
            //  |Z| between 1.28 and 12.7
            //
            else if (zabs <= 12.7)
            {
                y = a0 * z * z;
                pdf = Math.Exp(-y) * b0;

                q = pdf
                    / (zabs - b1 + b2
                        / (zabs + b3 + b4
                            / (zabs - b5 + b6
                                / (zabs + b7 - b8
                                    / (zabs + b9 + b10
                                        / (zabs + b11))))));
            }
            //
            //  Z far out in tail.
            //
            else
            {
                q = 0.0;
                pdf = 0.0;
            }

            if (z < 0.0)
            {
                p = q;
                q = 1.0 - p;
            }
            else
            {
                p = 1.0 - q;
            }
        }

        public static double ppchi2(double p, double v, double g, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PPCHI2 evaluates the percentage points of the Chi-squared PDF.
            //
            //  Discussion
            //
            //    Incorporates the suggested changes in AS R85 (vol.40(1),
            //    pages 233-5, 1991) which should eliminate the need for the limited
            //    range for P, though these limits have not been removed
            //    from the routine.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 June 2013
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Donald Best, DE Roberts.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Donald Best, DE Roberts,
            //    Algorithm AS 91:
            //    The Percentage Points of the Chi-Squared Distribution,
            //    Applied Statistics,
            //    Volume 24, Number 3, 1975, pages 385-390.
            //
            //  Parameters:
            //
            //    Input, double P,  value of the chi-squared cumulative
            //    probability density function.
            //    0.000002 <= P <= 0.999998.
            //
            //    Input, double V, the parameter of the chi-squared probability
            //    density function.
            //    0 < V.
            //
            //    Input, double G, the value of log ( Gamma ( V / 2 ) ).
            //
            //    Output, int *IFAULT, is nonzero if an error occurred.
            //    0, no error.
            //    1, P is outside the legal range.
            //    2, V is not positive.
            //    3, an error occurred in GAMMAD.
            //    4, the result is probably as accurate as the machine will allow.
            //
            //    Output, double PPCHI2, the value of the chi-squared random
            //    deviate with the property that the probability that a chi-squared random
            //    deviate with parameter V is less than or equal to PPCHI2 is P.
            //
        {
            double a;
            double aa = 0.6931471806;
            double c1 = 0.01;
            double c2 = 0.222222;
            double c3 = 0.32;
            double c4 = 0.4;
            double c5 = 1.24;
            double c6 = 2.2;
            double c7 = 4.67;
            double c8 = 6.66;
            double c9 = 6.73;
            double c10 = 13.32;
            double c11 = 60.0;
            double c12 = 70.0;
            double c13 = 84.0;
            double c14 = 105.0;
            double c15 = 120.0;
            double c16 = 127.0;
            double c17 = 140.0;
            double c18 = 175.0;
            double c19 = 210.0;
            double c20 = 252.0;
            double c21 = 264.0;
            double c22 = 294.0;
            double c23 = 346.0;
            double c24 = 420.0;
            double c25 = 462.0;
            double c26 = 606.0;
            double c27 = 672.0;
            double c28 = 707.0;
            double c29 = 735.0;
            double c30 = 889.0;
            double c31 = 932.0;
            double c32 = 966.0;
            double c33 = 1141.0;
            double c34 = 1182.0;
            double c35 = 1278.0;
            double c36 = 1740.0;
            double c37 = 2520.0;
            double c38 = 5040.0;
            double ch;
            double e = 0.5E-06;
            int if1 = 0;
            int maxit = 20;
            double pmax = 0.999998;
            double pmin = 0.000002;
            double p1;
            double p2;
            double q;
            double t;
            //
            //  Test arguments and initialize.
            //
            double value = -1.0;

            if (p < pmin || pmax < p)
            {
                ifault = 1;
                return value;
            }

            if (v <= 0.0)
            {
                ifault = 2;
                return value;
            }

            ifault = 0;
            double xx = 0.5 * v;
            double c = xx - 1.0;
            //
            //  Starting approximation for small chi-squared
            //
            if (v < -c5 * Math.Log(p))
            {
                ch = Math.Pow(p * xx * Math.Exp(g + xx * aa), 1.0 / xx);

                if (ch < e)
                {
                    value = ch;
                    return value;
                }
            }
            //
            //  Starting approximation for V less than or equal to 0.32
            //
            else if (v <= c3)
            {
                ch = c4;
                a = Math.Log(1.0 - p);

                for (;;)
                {
                    q = ch;
                    p1 = 1.0 + ch * (c7 + ch);
                    p2 = ch * (c9 + ch * (c8 + ch));

                    t = -0.5 + (c7 + 2.0 * ch) / p1 - (c9 + ch * (c10 +
                                                                  3.0 * ch)) / p2;

                    ch = ch - (1.0 - Math.Exp(a + g + 0.5 * ch + c * aa) * p2 / p1) / t;

                    if (Math.Abs(q / ch - 1.0) <= c1)
                    {
                        break;
                    }
                }
            }
            else
            {
                //
                //  Call to algorithm AS 111 - note that P has been tested above.
                //  AS 241 could be used as an alternative.
                //
                double x = ppnd(p, ref ifault);
                //
                //  Starting approximation using Wilson and Hilferty estimate
                //
                p1 = c2 / v;
                ch = v * Math.Pow(x * Math.Sqrt(p1) + 1.0 - p1, 3);
                //
                //  Starting approximation for P tending to 1.
                //
                if (c6 * v + 6.0 < ch)
                {
                    ch = -2.0 * (Math.Log(1.0 - p) - c * Math.Log(0.5 * ch) + g);
                }
            }

            //
            //  Call to algorithm AS 239 and calculation of seven term
            //  Taylor series
            //
            for (int i = 1; i <= maxit; i++)
            {
                q = ch;
                p1 = 0.5 * ch;
                p2 = p - gammad(p1, xx, ref if1);

                if (if1 != 0)
                {
                    ifault = 3;
                    return value;
                }

                t = p2 * Math.Exp(xx * aa + g + p1 - c * Math.Log(ch));
                double b = t / ch;
                a = 0.5 * t - b * c;
                double s1 = (c19 + a * (c17 + a * (c14 + a * (c13 + a * (c12 +
                                                                         c11 * a))))) / c24;
                double s2 = (c24 + a * (c29 + a * (c32 + a * (c33 + c35 * a)))) / c37;
                double s3 = (c19 + a * (c25 + a * (c28 + c31 * a))) / c37;
                double s4 = (c20 + a * (c27 + c34 * a) + c * (c22 + a * (c30 + c36 * a))) / c38;
                double s5 = (c13 + c21 * a + c * (c18 + c26 * a)) / c37;
                double s6 = (c15 + c * (c23 + c16 * c)) / c38;
                ch = ch + t * (1.0 + 0.5 * t * s1 - b * c * (s1 - b *
                    (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));

                if (e < Math.Abs(q / ch - 1.0))
                {
                    value = ch;
                    return value;
                }
            }

            ifault = 4;
            value = ch;

            return value;
        }

        public static double trigamma(double x, ref int ifault)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIGAMMA calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 January 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by BE Schneider.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    BE Schneider,
            //    Algorithm AS 121:
            //    Trigamma Function,
            //    Applied Statistics,
            //    Volume 27, Number 1, pages 97-99, 1978.
            //
            //  Parameters:
            //
            //    Input, double X, the argument of the trigamma function.
            //    0 < X.
            //
            //    Output, int *IFAULT, error flag.
            //    0, no error.
            //    1, X <= 0.
            //
            //    Output, double TRIGAMMA, the value of the trigamma function at X.
            //
        {
            double a = 0.0001;
            double b = 5.0;
            double b2 = 0.1666666667;
            double b4 = -0.03333333333;
            double b6 = 0.02380952381;
            double b8 = -0.03333333333;
            double value;
            double y;
            double z;
            //
            //  Check the input.
            //
            if (x <= 0.0)
            {
                ifault = 1;
                value = 0.0;
                return value;
            }

            ifault = 0;
            z = x;
            //
            //  Use small value approximation if X <= A.
            //
            if (x <= a)
            {
                value = 1.0 / x / x;
                return value;
            }

            //
            //  Increase argument to ( X + I ) >= B.
            //
            value = 0.0;

            while (z < b)
            {
                value = value + 1.0 / z / z;
                z = z + 1.0;
            }

            //
            //  Apply asymptotic formula if argument is B or greater.
            //
            y = 1.0 / z / z;

            value = value + 0.5 *
                y + (1.0
                     + y * (b2
                            + y * (b4
                                   + y * (b6
                                          + y * b8)))) / z;

            return value;
        }


    }
}