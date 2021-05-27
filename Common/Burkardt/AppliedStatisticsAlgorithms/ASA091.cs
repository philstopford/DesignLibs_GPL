using System;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static void chi_square_cdf_values(ref int n_data, ref int a, ref double x, ref double fx)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = ChiSquareDistribution [ df ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 August 2004
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
                0.7965567455405796E-01,
                0.4987520807317687E-02,
                0.1124629160182849E+00,
                0.9950166250831946E-02,
                0.4729107431344619E+00,
                0.1812692469220181E+00,
                0.5975750516063926E-01,
                0.1752309630642177E-01,
                0.6826894921370859E+00,
                0.3934693402873666E+00,
                0.1987480430987992E+00,
                0.9020401043104986E-01,
                0.3743422675270363E-01,
                0.4275932955291202E+00,
                0.6083748237289110E+00,
                0.7385358700508894E+00,
                0.8282028557032669E+00,
                0.8883897749052874E+00,
                0.1721156299558408E-03,
                0.3659846827343712E-02,
                0.1857593622214067E-01
            }
            ;

            double[] x_vec =  {
                0.01E+00,
                0.01E+00,
                0.02E+00,
                0.02E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                5.00E+00,
                6.00E+00,
                1.00E+00,
                2.00E+00,
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




    }
}