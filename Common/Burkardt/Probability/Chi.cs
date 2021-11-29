using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Chi
{
    public static double chi_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_CDF evaluates the Chi CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x <= a)
        {
            cdf = 0.0;
        }
        else
        {
            double y = (x - a) / b;
            double x2 = 0.5 * y * y;
            double p2 = 0.5 * c;

            cdf = typeMethods.r8_gamma_inc(p2, x2);
        }

        return cdf;
    }

    public static double chi_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_CDF_INV inverts the Chi CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double CHI_CDF_INV, the corresponding argument of the CDF.
        //
    {
        const int it_max = 100;
        const double tol = 0.0001;

        double x = 0.0;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = a;
                return x;
            case 1.0:
                x = typeMethods.r8_huge();
                return x;
        }

        double x1 = a;
        double cdf1 = 0.0;

        double x2 = a + 1.0;

        for (;;)
        {
            double cdf2 = chi_cdf(x2, a, b, c);

            if (cdf < cdf2)
            {
                break;
            }

            x2 = a + 2.0 * (x2 - a);
        }

        //
        //  Now use bisection.
        //
        int it = 0;

        for (;;)
        {
            it += 1;

            double x3 = 0.5 * (x1 + x2);
            double cdf3 = chi_cdf(x3, a, b, c);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("CHI_CDF_INV - Fatal error!");
                Console.WriteLine("  Iteration limit exceeded.");
                return x;
            }

            if (cdf3 <= cdf && cdf1 <= cdf || cdf <= cdf3 && cdf <= cdf1)
            {
                x1 = x3;
                cdf1 = cdf3;
            }
            else
            {
                x2 = x3;
            }

        }

        return x;
    }

    public static bool chi_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_CHECK checks the parameters of the Chi CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, bool CHI_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0.");
                return false;
        }

        switch (c)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_CHECK - Warning!");
                Console.WriteLine("  C <= 0.0.");
                return false;
            default:
                return true;
        }
    }

    public static double chi_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_MEAN returns the mean of the Chi PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double MEAN, the mean value.
        //
    {
        double mean = a + Math.Sqrt(2.0) * b * Helpers.Gamma(0.5 * (c + 1.0))
            / Helpers.Gamma(0.5 * c);

        return mean;
    }

    public static double chi_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_PDF evaluates the Chi PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = EXP ( - 0.5 * ( ( X - A ) / B )^2 )
        //      * ( ( X - A ) / B )^( C - 1 ) /
        //      ( 2^( 0.5 * C - 1 ) * B * GAMMA ( 0.5 * C ) )
        //
        //    CHI(A,B,1) is the Half Normal PDF;
        //    CHI(0,B,2) is the Rayleigh PDF;
        //    CHI(0,B,3) is the Maxwell PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x <= a)
        {
            pdf = 0.0;
        }
        else
        {
            double y = (x - a) / b;

            pdf = Math.Exp(-0.5 * y * y) * Math.Pow(y, c - 1.0) /
                  (Math.Pow(2.0, 0.5 * c - 1.0) * b * Helpers.Gamma(0.5 * c));
        }

        return pdf;
    }

    public static double chi_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SAMPLE samples the Chi PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CHI_SAMPLE, a sample of the PDF.
        //
    {
        double x = chi_square_sample(c, ref seed);

        x = a + b * Math.Sqrt(x);

        return x;
    }

    public static double chi_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_VARIANCE returns the variance of the Chi PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double variance = b * b * (c - 2.0 *
            Math.Pow(Helpers.Gamma(0.5 * (c + 1.0)) / Helpers.Gamma(0.5 * c), 2));

        return variance;
    }

    public static double chi_square_cdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_CDF evaluates the Chi squared CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the value of the random deviate.
        //
        //    Input, double A, the parameter of the distribution, usually
        //    the number of degrees of freedom.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        double x2 = 0.5 * x;
        const double a2 = 0.0;
        const double b2 = 1.0;
        double c2 = 0.5 * a;

        double cdf = Gamma.gamma_cdf(x2, a2, b2, c2);

        return cdf;
    }

    public static double chi_square_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_CDF_INV inverts the Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Best and Roberts.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Best, Roberts,
        //    The Percentage Points of the Chi-Squared Distribution,
        //    Algorithm AS 91,
        //    Applied Statistics,
        //    Volume 24, Number ?, pages 385-390, 1975.
        //
        //  Parameters:
        //
        //    Input, double CDF, a value of the chi-squared cumulative
        //    probability density function.
        //    0.000002 <= CDF <= 0.999998.
        //
        //    Input, double A, the parameter of the chi-squared
        //    probability density function.  0 < A.
        //
        //    Output, double CHI_SQUARE_CDF_INV, the value of the chi-squared random deviate
        //    with the property that the probability that a chi-squared random
        //    deviate with parameter A is less than or equal to X is CDF.
        //
    {
        double a2;
        const double aa = 0.6931471806;
        const double c1 = 0.01;
        const double c2 = 0.222222;
        const double c3 = 0.32;
        const double c4 = 0.4;
        const double c5 = 1.24;
        const double c6 = 2.2;
        const double c7 = 4.67;
        const double c8 = 6.66;
        const double c9 = 6.73;
        const double c10 = 13.32;
        const double c11 = 60.0;
        const double c12 = 70.0;
        const double c13 = 84.0;
        const double c14 = 105.0;
        const double c15 = 120.0;
        const double c16 = 127.0;
        const double c17 = 140.0;
        const double c18 = 175.0;
        const double c19 = 210.0;
        const double c20 = 252.0;
        const double c21 = 264.0;
        const double c22 = 294.0;
        const double c23 = 346.0;
        const double c24 = 420.0;
        const double c25 = 462.0;
        const double c26 = 606.0;
        const double c27 = 672.0;
        const double c28 = 707.0;
        const double c29 = 735.0;
        const double c30 = 889.0;
        const double c31 = 932.0;
        const double c32 = 966.0;
        const double c33 = 1141.0;
        const double c34 = 1182.0;
        const double c35 = 1278.0;
        const double c36 = 1740.0;
        const double c37 = 2520.0;
        const double c38 = 5040.0;
        const double cdf_max = 0.999998;
        const double cdf_min = 0.000002;
        double ch;
        const double e = 0.0000005;
        int i;
        const int it_max = 20;
        double p1;
        double p2;
        double q;
        double t;
        double x;
        switch (cdf)
        {
            //
            case < cdf_min:
                x = -1.0;
                Console.WriteLine(" ");
                Console.WriteLine("CHI_SQUARE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < CDF_MIN.");
                return 1;
            case > cdf_max:
                x = -1.0;
                Console.WriteLine(" ");
                Console.WriteLine("CHI_SQUARE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF_MAX < CDF.");
                return 1;
        }

        double xx = 0.5 * a;
        double c = xx - 1.0;
        //
        //  Compute Log ( Gamma ( A/2 ) ).
        //
        double g = Helpers.LogGamma(a / 2.0);
        //
        //  Starting approximation for small chi-squared.
        //
        if (a < -c5 * Math.Log(cdf))
        {
            ch = Math.Pow(cdf * xx * Math.Exp(g + xx * aa), 1.0 / xx);

            if (ch < e)
            {
                x = ch;
                return x;
            }
        }
        //
        //  Starting approximation for A less than or equal to 0.32.
        //
        else if (a <= c3)
        {
            ch = c4;
            a2 = Math.Log(1.0 - cdf);

            for (;;)
            {
                q = ch;
                p1 = 1.0 + ch * (c7 + ch);
                p2 = ch * (c9 + ch * (c8 + ch));

                t = -0.5 + (c7 + 2.0 * ch) / p1 - (c9 + ch * (c10 + 3.0 * ch)) / p2;

                ch -= (1.0 - Math.Exp(a2 + g + 0.5 * ch + c * aa) * p2 / p1) / t;

                if (Math.Abs(q / ch - 1.0) <= c1)
                {
                    break;
                }
            }
        }
        //
        //  Call to algorithm AS 111.
        //  Note that P has been tested above.
        //  AS 241 could be used as an alternative.
        //
        else
        {
            double x2 = Normal.normal_01_cdf_inv(cdf);
            //
            //  Starting approximation using Wilson and Hilferty estimate.
            //
            p1 = c2 / a;
            ch = a * Math.Pow(x2 * Math.Sqrt(p1) + 1.0 - p1, 3);
            //
            //  Starting approximation for P tending to 1.
            //
            if (c6 * a + 6.0 < ch)
            {
                ch = -2.0 * (Math.Log(1.0 - cdf) - c * Math.Log(0.5 * ch) + g);
            }
        }

        //
        //  Call to algorithm AS 239 and calculation of seven term Taylor series.
        //
        for (i = 1; i <= it_max; i++)
        {
            q = ch;
            p1 = 0.5 * ch;
            p2 = cdf - typeMethods.r8_gamma_inc(xx, p1);
            t = p2 * Math.Exp(xx * aa + g + p1 - c * Math.Log(ch));
            double b = t / ch;
            a2 = 0.5 * t - b * c;

            double s1 = (c19 + a2
                * (c17 + a2
                    * (c14 + a2
                        * (c13 + a2
                            * (c12 + a2
                                * c11))))) / c24;

            double s2 = (c24 + a2
                * (c29 + a2
                    * (c32 + a2
                        * (c33 + a2
                            * c35)))) / c37;

            double s3 = (c19 + a2
                * (c25 + a2
                    * (c28 + a2
                        * c31))) / c37;

            double s4 = (c20 + a2
                * (c27 + a2
                    * c34) + c
                * (c22 + a2
                    * (c30 + a2
                        * c36))) / c38;

            double s5 = (c13 + c21 * a2 + c * (c18 + c26 * a2)) / c37;

            double s6 = (c15 + c * (c23 + c16 * c)) / c38;

            ch += t * (1.0 + 0.5 * t * s1 - b * c
                                              * (s1 - b
                                                  * (s2 - b
                                                      * (s3 - b
                                                          * (s4 - b
                                                              * (s5 - b
                                                                  * s6))))));

            if (!(e < Math.Abs(q / ch - 1.0)))
            {
                continue;
            }

            x = ch;
            return x;

        }

        x = ch;
        Console.WriteLine(" ");
        Console.WriteLine("CHI_SQUARE_CDF_INV - Warning!");
        Console.WriteLine("  Convergence not reached.");

        return x;
    }

    public static void chi_square_cdf_values(ref int n_data, ref int a, ref double x, ref double fx )
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &A, the parameter of the function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 21;

        int[] a_vec =
            {
                1, 2, 1, 2,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 3, 3, 3,
                3, 3, 10, 10,
                10
            }
            ;

        double[] fx_vec =
            {
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

        double[] x_vec =
            {
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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static bool chi_square_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_CHECK checks the parameter of the central Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the distribution.
        //    1 <= A.
        //
        //    Output, bool CHI_SQUARE_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case < 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_SQUARE_CHECK - Warning!");
                Console.WriteLine("  A < 1.0.");
                return false;
            default:
                return true;
        }
    }

    public static double chi_square_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_MEAN returns the mean of the central Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the distribution.
        //    1 <= A.
        //
        //    Output, double MEAN, the mean value.
        //
    {
        return a;
    }

    public static double chi_square_pdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_PDF evaluates the central Chi squared PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) =
        //      EXP ( - X / 2 ) * X^((A-2)/2) / ( 2^(A/2) * GAMMA ( A/2 ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X
        //
        //    Input, double A, the parameter of the PDF.
        //    1 <= A.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case < 0.0:
                pdf = 0.0;
                break;
            default:
            {
                double b = a / 2.0;
                pdf = Math.Exp(-0.5 * x) * Math.Pow(x, b - 1.0) / (Math.Pow(2.0, b)
                                                                   * Helpers.Gamma(b));
                break;
            }
        }

        return pdf;
    }

    public static double chi_square_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_SAMPLE samples the central Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1 <= A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CHI_SQUARE_SAMPLE, a sample of the PDF.
        //
    {
        const int it_max = 100;
        double x;

        int n = (int) a;

        if (Math.Abs(n - a) <= double.Epsilon && n <= it_max)
        {
            x = 0.0;
            int i;
            for (i = 1; i <= n; i++)
            {
                double x2 = Normal.normal_01_sample(ref seed);
                x += x2 * x2;
            }
        }
        else
        {
            const double a2 = 0.0;
            const double b2 = 1.0;
            double c2 = a / 2.0;

            x = Gamma.gamma_sample(a2, b2, c2, ref seed);

            x = 2.0 * x;
        }

        return x;
    }

    public static double chi_square_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_VARIANCE returns the variance of the central Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the distribution.
        //    1 <= A.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double variance = 2.0 * a;

        return variance;
    }

    public static void chi_square_noncentral_cdf_values(ref int n_data, ref int df, ref double lambda,
            ref double x, ref double cdf )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralChiSquareDistribution [ df, lambda ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 August 2004
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
        //    Output, int &DF, the number of degrees of freedom.
        //
        //    Output, double &LAMBDA, the noncentrality parameter.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &CDF, the noncentral chi CDF.
        //
    {
        const int N_MAX = 28;

        double[] cdf_vec =
            {
                0.8399444269398261E+00,
                0.6959060300435139E+00,
                0.5350879697078847E+00,
                0.7647841496310313E+00,
                0.6206436532195436E+00,
                0.4691667375373180E+00,
                0.3070884345937569E+00,
                0.2203818092990903E+00,
                0.1500251895581519E+00,
                0.3071163194335791E-02,
                0.1763982670131894E-02,
                0.9816792594625022E-03,
                0.1651753140866208E-01,
                0.2023419573950451E-03,
                0.4984476352854074E-06,
                0.1513252400654827E-01,
                0.2090414910614367E-02,
                0.2465021206048452E-03,
                0.2636835050342939E-01,
                0.1857983220079215E-01,
                0.1305736595486640E-01,
                0.5838039534819351E-01,
                0.4249784402463712E-01,
                0.3082137716021596E-01,
                0.1057878223400849E+00,
                0.7940842984598509E-01,
                0.5932010895599639E-01,
                0.2110395656918684E+00
            }
            ;

        int[] df_vec =
            {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                60, 80, 100,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10,
                8
            }
            ;

        double[] lambda_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                30.0E+00,
                30.0E+00,
                30.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                0.5E+00
            }
            ;

        double[] x_vec =
            {
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                60.000E+00,
                60.000E+00,
                60.000E+00,
                0.050E+00,
                0.050E+00,
                0.050E+00,
                4.000E+00,
                4.000E+00,
                4.000E+00,
                5.000E+00,
                5.000E+00,
                5.000E+00,
                6.000E+00,
                6.000E+00,
                6.000E+00,
                5.000E+00
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
            x = 0.0;
            lambda = 0.0;
            df = 0;
            cdf = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            lambda = lambda_vec[n_data - 1];
            df = df_vec[n_data - 1];
            cdf = cdf_vec[n_data - 1];
        }
    }

    public static bool chi_square_noncentral_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_NONCENTRAL_CHECK checks the parameters of the noncentral Chi Squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the parameter of the PDF.
        //    1.0 <= A.
        //
        //    Input, double B, the noncentrality parameter of the PDF.
        //    0.0 <= B.
        //
        //    Output, bool CHI_SQUARE_NONCENTRAL_CHECK, is true if the parameters
        //    are legal.
        //
    {
        switch (a)
        {
            case < 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_SQUARE_NONCENTRAL_CHECK - Warning!");
                Console.WriteLine("  A < 1.");
                return false;
        }

        switch (b)
        {
            case < 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("CHI_SQUARE_NONCENTRAL_CHECK - Warning!");
                Console.WriteLine("  B < 0.");
                return false;
            default:
                return true;
        }
    }

    public static double chi_square_noncentral_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_NONCENTRAL_MEAN returns the mean of the noncentral Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the parameter of the PDF.
        //    1.0 <= A.
        //
        //    Input, double B, the noncentrality parameter of the PDF.
        //    0.0 <= B.
        //
        //    Output, double MEAN, the mean value.
        //
    {
        double mean = a + b;

        return mean;
    }

    public static double chi_square_noncentral_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_NONCENTRAL_SAMPLE samples the noncentral Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the parameter of the PDF.
        //    1.0 <= A.
        //
        //    Input, double B, the noncentrality parameter of the PDF.
        //    0.0 <= B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CHI_SQUARE_NONCENTRAL_SAMPLE, a sample of the PDF.
        //
    {
        double a1 = a - 1.0;

        double x1 = chi_square_sample(a1, ref seed);

        double a2 = Math.Sqrt(b);
        const double b2 = 1.0;
        double x2 = Normal.normal_sample(a2, b2, ref seed);

        double x = x1 + x2 * x2;

        return x;
    }

    public static double chi_square_noncentral_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHI_SQUARE_NONCENTRAL_VARIANCE: variance of the noncentral Chi squared PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1 <= A.
        //
        //    Input, double B, the noncentrality parameter of the PDF.
        //    0.0 <= B.
        //
        //    Output, double VARIANCE, the variance value.
        //
    {
        double variance = 2.0 * (a + 2.0 * b);

        return variance;
    }

}