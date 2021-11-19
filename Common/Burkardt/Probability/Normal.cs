using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Normal
{
    public static double normal_01_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_CDF evaluates the Normal 01 CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 February 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    A G Adams,
        //    Areas Under the Normal Curve,
        //    Algorithm 39,
        //    Computer j.,
        //    Volume 12, pages 197-198, 1969.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double CDF, the value of the CDF.
        //
    {
        const double a1 = 0.398942280444;
        const double a2 = 0.399903438504;
        const double a3 = 5.75885480458;
        const double a4 = 29.8213557808;
        const double a5 = 2.62433121679;
        const double a6 = 48.6959930692;
        const double a7 = 5.92885724438;
        const double b0 = 0.398942280385;
        const double b1 = 3.8052E-08;
        const double b2 = 1.00000615302;
        const double b3 = 3.98064794E-04;
        const double b4 = 1.98615381364;
        const double b5 = 0.151679116635;
        const double b6 = 5.29330324926;
        const double b7 = 4.8385912808;
        const double b8 = 15.1508972451;
        const double b9 = 0.742380924027;
        const double b10 = 30.789933034;
        const double b11 = 3.99019417011;
        double q;
        double y;
        switch (Math.Abs(x))
        {
            //
            //  |X| <= 1.28.
            //
            case <= 1.28:
                y = 0.5 * x * x;

                q = 0.5 - Math.Abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5
                                                                         + a6 / (y + a7))));
                //
                //  1.28 < |X| <= 12.7
                //
                break;
            case <= 12.7:
                y = 0.5 * x * x;

                q = Math.Exp(-y) * b0 / (Math.Abs(x) - b1
                                         + b2 / (Math.Abs(x) + b3
                                                             + b4 / (Math.Abs(x) - b5
                                                                     + b6 / (Math.Abs(x) + b7
                                                                             - b8 / (Math.Abs(x) + b9
                                                                                 + b10 / (Math.Abs(x) + b11))))));
                //
                //  12.7 < |X|
                //
                break;
            default:
                q = 0.0;
                break;
        }

        double cdf = x switch
        {
            //
            //  Take account of negative X.
            //
            < 0.0 => q,
            _ => 1.0 - q
        };

        return cdf;
    }

    public static double normal_01_cdf_inv(double p)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_CDF_INV inverts the standard normal CDF.
        //
        //  Discussion:
        //
        //    The result is accurate to about 1 part in 10**16.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 December 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Wichura.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Wichura,
        //    The Percentage Points of the Normal Distribution,
        //    Algorithm AS 241,
        //    Applied Statistics,
        //    Volume 37, Number 3, pages 477-484, 1988.
        //
        //  Parameters:
        //
        //    Input, double P, the value of the cumulative probability
        //    densitity function.  0 < P < 1.  If P is outside this range, an
        //    "infinite" value is returned.
        //
        //    Output, double NORMAL_01_CDF_INV, the normal deviate value
        //    with the property that the probability of a standard normal deviate being
        //    less than or equal to this value is P.
        //
    {
        double[] a =  {
                3.3871328727963666080, 1.3314166789178437745e+2,
                1.9715909503065514427e+3, 1.3731693765509461125e+4,
                4.5921953931549871457e+4, 6.7265770927008700853e+4,
                3.3430575583588128105e+4, 2.5090809287301226727e+3
            }
            ;
        double[] b =  {
                1.0, 4.2313330701600911252e+1,
                6.8718700749205790830e+2, 5.3941960214247511077e+3,
                2.1213794301586595867e+4, 3.9307895800092710610e+4,
                2.8729085735721942674e+4, 5.2264952788528545610e+3
            }
            ;
        double[] c =  {
                1.42343711074968357734, 4.63033784615654529590,
                5.76949722146069140550, 3.64784832476320460504,
                1.27045825245236838258, 2.41780725177450611770e-1,
                2.27238449892691845833e-2, 7.74545014278341407640e-4
            }
            ;
        double const1 = 0.180625;
        double const2 = 1.6;
        double[] d =  {
                1.0, 2.05319162663775882187,
                1.67638483018380384940, 6.89767334985100004550e-1,
                1.48103976427480074590e-1, 1.51986665636164571966e-2,
                5.47593808499534494600e-4, 1.05075007164441684324e-9
            }
            ;
        double[] e =  {
                6.65790464350110377720, 5.46378491116411436990,
                1.78482653991729133580, 2.96560571828504891230e-1,
                2.65321895265761230930e-2, 1.24266094738807843860e-3,
                2.71155556874348757815e-5, 2.01033439929228813265e-7
            }
            ;
        double[] f =  {
                1.0, 5.99832206555887937690e-1,
                1.36929880922735805310e-1, 1.48753612908506148525e-2,
                7.86869131145613259100e-4, 1.84631831751005468180e-5,
                1.42151175831644588870e-7, 2.04426310338993978564e-15
            }
            ;
        double r;
        const double r8_huge = 1.0E+30;
        const double split1 = 0.425;
        const double split2 = 5.0;
        double value;

        switch (p)
        {
            case <= 0.0:
                value = -r8_huge;
                return value;
            case >= 1.0:
                value = r8_huge;
                return value;
        }

        double q = p - 0.5;

        if (Math.Abs(q) <= split1)
        {
            r = const1 - q * q;
            value = q * typeMethods.r8poly_value_horner(7, a, r) / typeMethods.r8poly_value_horner(7, b, r);
        }
        else
        {
            r = q switch
            {
                < 0.0 => p,
                _ => 1.0 - p
            };

            switch (r)
            {
                case <= 0.0:
                    value = r8_huge;
                    break;
                default:
                {
                    r = Math.Sqrt(-Math.Log(r));

                    if (r <= split2)
                    {
                        r -= const2;
                        value = typeMethods.r8poly_value_horner(7, c, r) / typeMethods.r8poly_value_horner(7, d, r);
                    }
                    else
                    {
                        r -= split2;
                        value = typeMethods.r8poly_value_horner(7, e, r) / typeMethods.r8poly_value_horner(7, f, r);
                    }

                    break;
                }
            }

            value = q switch
            {
                < 0.0 => -value,
                _ => value
            };
        }

        return value;
    }


    public static double normal_01_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_MEAN returns the mean of the Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double MEAN, the mean of the PDF.
        //
    {
        double mean = 0.0;

        return mean;
    }
        
    public static double normal_01_moment ( int order )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_MOMENT evaluates moments of the Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the moment.
        //    0 <= ORDER.
        //
        //    Output, double NORMAL_01_MOMENT, the value of the moment.
        //
    {
        double value = (order % 2) switch
        {
            0 => typeMethods.r8_factorial2(order - 1),
            _ => 0.0
        };

        return value;
    }

    public static double normal_01_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_PDF evaluates the Normal 01 PDF.
        //
        //  Discussion:
        //
        //    The Normal 01 PDF is also called the "Standard Normal" PDF, or
        //    the Normal PDF with 0 mean and variance 1.
        //
        //    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
            

        double pdf = Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);

        return pdf;
    }

    public static double normal_01_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_SAMPLE samples the standard normal probability distribution.
        //
        //  Discussion:
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //    The Box-Muller method is used, which is efficient, but
        //    generates two values at a time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_01_SAMPLE, a normally distributed random value.
        //
    {
            

        double r1 = UniformRNG.r8_uniform_01(ref seed);
        double r2 = UniformRNG.r8_uniform_01(ref seed);
        double x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * Math.PI * r2);

        return x;
    }

    public static double[] normal_01_samples(int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_SAMPLES returns multiple samples of the standard normal PDF.
        //
        //  Discussion:
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //    The Box-Muller method is used, which is efficient, but
        //    generates two values at a time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of samples.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_01_SAMPLES[N], the samples.
        //
    {
        double[] r1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double[] r2 = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        double[] x = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = Math.Sqrt(-2.0 * Math.Log(r1[i])) * Math.Cos(2.0 * Math.PI * r2[i]);
        }
        return x;
    }

    public static double normal_01_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_VARIANCE returns the variance of the Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        const double variance = 1.0;

        return variance;
    }

    public static double[] normal_01_vector(int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_01_VECTOR samples the standard normal probability distribution.
        //
        //  Discussion:
        //
        //    The standard normal probability distribution function (PDF) has
        //    mean 0 and standard deviation 1.
        //
        //    This routine can generate a vector of values on one call.  It
        //    has the feature that it should provide the same results
        //    in the same order no matter how we break up the task.
        //
        //    Before calling this routine, the user may call RANDOM_SEED
        //    in order to set the seed of the random number generator.
        //
        //    The Box-Muller method is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values desired.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double X(N), a sample of the standard normal PDF.
        //
        //  Local parameters:
        //
        //    Local, real R(N+1), is used to store some uniform random values.
        //    Its dimension is N+1, but really it is only needed to be the
        //    smallest even number greater than or equal to N.
        //
        //    Local, int X_LO, X_HI, records the range of entries of
        //    X that we need to compute.  This starts off as 1:N, but is adjusted
        //    if we have a saved value that can be immediately stored in X(1),
        //    and so on.
        //
    {
        double[] r;

        double[] x = new double[n];
        //
        //  Record the range of X we need to fill in.
        //
        int x_lo = 1;
        int x_hi = n;
        switch (x_hi - x_lo + 1)
        {
            //
            //  If we need just one new value, do that here to avoid null arrays.
            //
            case 1:
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
                break;
            //
            default:
            {
                int i;
                int m;
                switch ((x_hi - x_lo + 1) % 2)
                {
                    case 0:
                    {
                        m = (x_hi - x_lo + 1) / 2;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 2; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        break;
                    }
                    //
                    default:
                    {
                        x_hi -= 1;

                        m = (x_hi - x_lo + 1) / 2 + 1;

                        r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                        for (i = 0; i <= 2 * m - 4; i += 2)
                        {
                            x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                            x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                        }

                        i = 2 * m - 2;

                        x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                        break;
                    }
                }

                break;
            }
        }

        return x;
    }

    public static double normal_ms_mean(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_MEAN returns the mean of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_MS_MEAN, the mean of the PDF.
        //
    {
        double mean = mu;

        return mean;
    }

    public static double normal_ms_moment(int order, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_MOMENT evaluates moments of the Normal PDF.
        //
        //  Discussion:
        //
        //    The formula was posted by John D Cook.
        //
        //    Order  Moment
        //    -----  ------
        //      0    1
        //      1    mu
        //      2    mu^2 +         sigma^2
        //      3    mu^3 +  3 mu   sigma^2
        //      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
        //      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
        //      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
        //      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
        //      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6 
        //           + 105 sigma^8
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the moment.
        //    0 <= ORDER.
        //
        //    Input, double MU, the mean of the distribution.
        //
        //    Input, double SIGMA, the standard deviation of the distribution.
        //
        //    Output, double NORMAL_MS_MOMENT, the value of the central moment.
        //
    {
        int j;

        int j_hi = order / 2;

        double value = 0.0;
        for (j = 0; j <= j_hi; j++)
        {
            value += typeMethods.r8_choose(order, 2 * j)
                     * typeMethods.r8_factorial2(2 * j - 1)
                     * Math.Pow(mu, order - 2 * j) * Math.Pow(sigma, 2 * j);
        }

        return value;
    }

    public static double normal_ms_moment_central(int order, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_MOMENT_CENTRAL evaluates central moments of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the moment.
        //    0 <= ORDER.
        //
        //    Input, double MU, the mean of the distribution.
        //
        //    Input, double SIGMA, the standard deviation of the distribution.
        //
        //    Output, double NORMAL_MS_MOMENT_CENTRAL, the value of the central moment.
        //
    {
        double value = (order % 2) switch
        {
            0 => typeMethods.r8_factorial2(order - 1) * Math.Pow(sigma, order),
            _ => 0.0
        };

        return value;
    }

    public static double normal_ms_moment_central_values(int order, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_MOMENT_CENTRAL_VALUES: moments 0 through 10 of the Normal PDF.
        //
        //  Discussion:
        //
        //    The formula was posted by John D Cook.
        //
        //    Order  Moment
        //    -----  ------
        //      0    1
        //      1    0
        //      2    sigma^2
        //      3    0
        //      4    3 sigma^4
        //      5    0
        //      6    15 sigma^6
        //      7    0
        //      8    105 sigma^8
        //      9    0
        //     10    945 sigma^10
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the moment.
        //    0 <= ORDER <= 10.
        //
        //    Input, double MU, the mean of the distribution.
        //
        //    Input, double SIGMA, the standard deviation of the distribution.
        //
        //    Output, double NORMAL_MS_MOMENT_CENTRAL_VALUES, the value of the 
        //    central moment.
        //
    {
        double value;

        switch (order)
        {
            case 0:
                value = 1.0;
                break;
            case 1:
                value = 0.0;
                break;
            case 2:
                value = Math.Pow(sigma, 2);
                break;
            case 3:
                value = 0.0;
                break;
            case 4:
                value = 3.0 * Math.Pow(sigma, 4);
                break;
            case 5:
                value = 0.0;
                break;
            case 6:
                value = 15.0 * Math.Pow(sigma, 6);
                break;
            case 7:
                value = 0.0;
                break;
            case 8:
                value = 105.0 * Math.Pow(sigma, 8);
                break;
            case 9:
                value = 0.0;
                break;
            case 10:
                value = 945.0 * Math.Pow(sigma, 10);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_MS_MOMENT_CENTRAL_VALUES - Fatal error!");
                Console.WriteLine("  Only ORDERS 0 through 10 are available.");
                return 1;
        }

        return value;
    }

    public static double normal_ms_moment_values(int order, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_MOMENT_VALUES evaluates moments 0 through 8 of the Normal PDF.
        //
        //  Discussion:
        //
        //    The formula was posted by John D Cook.
        //
        //    Order  Moment
        //    -----  ------
        //      0    1
        //      1    mu
        //      2    mu^2 +         sigma^2
        //      3    mu^3 +  3 mu   sigma^2
        //      4    mu^4 +  6 mu^2 sigma^2 +   3      sigma^4
        //      5    mu^5 + 10 mu^3 sigma^2 +  15 mu   sigma^4
        //      6    mu^6 + 15 mu^4 sigma^2 +  45 mu^2 sigma^4 +  15      sigma^6
        //      7    mu^7 + 21 mu^5 sigma^2 + 105 mu^3 sigma^4 + 105 mu   sigma^6
        //      8    mu^8 + 28 mu^6 sigma^2 + 210 mu^4 sigma^4 + 420 mu^2 sigma^6 
        //           + 105 sigma^8
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the moment.
        //    0 <= ORDER <= 8.
        //
        //    Input, double MU, the mean of the distribution.
        //
        //    Input, double SIGMA, the standard deviation of the distribution.
        //
        //    Output, double NORMAL_MS_MOMENT_VALUES, the value of the central moment.
        //
    {
        double value;

        switch (order)
        {
            case 0:
                value = 1.0;
                break;
            case 1:
                value = mu;
                break;
            case 2:
                value = Math.Pow(mu, 2) + Math.Pow(sigma, 2);
                break;
            case 3:
                value = Math.Pow(mu, 3) + 3.0 * mu * Math.Pow(sigma, 2);
                break;
            case 4:
                value = Math.Pow(mu, 4) + 6.0 * Math.Pow(mu, 2) * Math.Pow(sigma, 2)
                                        + 3.0 * Math.Pow(sigma, 4);
                break;
            case 5:
                value = Math.Pow(mu, 5) + 10.0 * Math.Pow(mu, 3) * Math.Pow(sigma, 2)
                                        + 15.0 * mu * Math.Pow(sigma, 4);
                break;
            case 6:
                value = Math.Pow(mu, 6) + 15.0 * Math.Pow(mu, 4) * Math.Pow(sigma, 2)
                                        + 45.0 * Math.Pow(mu, 2) * Math.Pow(sigma, 4)
                                        + 15.0 * Math.Pow(sigma, 6);
                break;
            case 7:
                value = Math.Pow(mu, 7) + 21.0 * Math.Pow(mu, 5) * Math.Pow(sigma, 2)
                                        + 105.0 * Math.Pow(mu, 3) * Math.Pow(sigma, 4)
                                        + 105.0 * mu * Math.Pow(sigma, 6);
                break;
            case 8:
                value = Math.Pow(mu, 8) + 28.0 * Math.Pow(mu, 6) * Math.Pow(sigma, 2)
                                        + 210.0 * Math.Pow(mu, 4) * Math.Pow(sigma, 4)
                                        + 420.0 * Math.Pow(mu, 2) * Math.Pow(sigma, 6) + 105.0 * Math.Pow(sigma, 8);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_MS_MOMENT_VALUES - Fatal error!");
                Console.WriteLine("  Only ORDERS 0 through 8 are available.");
                return 1;
        }

        return value;
    }

    public static double normal_ms_pdf(double x, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_PDF evaluates the Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(MU,SIGMA;X)
        //      = exp ( - 0.5 * ( ( X - MU ) / SIGMA )^2 )
        //      / ( SIGMA * SQRT ( 2 * PI ) )
        //
        //    The normal PDF is also known as the Gaussian PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_MS_PDF, the value of the PDF.
        //
    {
        double y = (x - mu) / sigma;

        double pdf = Math.Exp(-0.5 * y * y) / (sigma * Math.Sqrt(2.0 * Math.PI));

        return pdf;
    }

    public static double normal_ms_sample(double mu, double sigma, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_SAMPLE samples the Normal PDF.
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
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_MS_SAMPLE, a sample of the PDF.
        //
    {
        double x = normal_01_sample(ref seed);

        x = mu + sigma * x;

        return x;
    }

    public static double normal_ms_variance(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MS_VARIANCE returns the variance of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the parameters of the PDF.
        //    0.0 < SIGMA.
        //
        //    Output, double NORMAL_MS_VARIANCE, the variance of the PDF.
        //
    {
        double variance = sigma * sigma;

        return variance;
    }

    public static bool normal_check(double mu, double sigma)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CHECK checks the parameters of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, bool NORMAL_CHECK, is true if the parameters are legal.
        //
    {
        switch (sigma)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("NORMAL_CHECK - Warning!");
                Console.WriteLine("  SIGMA == 0.");
                return false;
            default:
                return true;
        }
    }

    public static double normal_mean(double mu, double sigma)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_MEAN returns the mean of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
    {
        return mu;
    }

    public static double normal_pdf(double x, double mu, double sigma)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_PDF evaluates the Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(X;MU,SIGMA)
        //      = exp ( - 0.5 * ( ( X - MU ) / SIGMAB )^2 )
        //      / SQRT ( 2 * PI * SIGMA^2 )
        //
        //    The normal PDF is also known as the Gaussian PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
            

        double y = (x - mu) / sigma;

        double pdf = Math.Exp(-0.5 * y * y) / Math.Sqrt(2.0 * Math.PI * sigma * sigma);

        return pdf;
    }

    public static double normal_sample(double mu, double sigma, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SAMPLE samples the Normal PDF.
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
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_SAMPLE, a sample of the PDF.
        //
    {
        double x = normal_01_sample(ref seed);

        x = mu + sigma * x;

        return x;
    }

    public static double[] normal_samples(int n, double mu, double sigma, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_SAMPLES returns multiple samples of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 September 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of samples.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_SAMPLES[N], the samples.
        //
    {
        double[] x = normal_01_samples(n, ref seed);

        for (int i = 0; i < n; i++)
        {
            x[i] = mu + sigma * x[i];
        }

        return x;
    }

    public static double normal_variance(double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_VARIANCE returns the variance of the Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double variance = sigma * sigma;

        return variance;
    }

    public static double[] normal_vector(int n, double mu, double sigma, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_VECTOR samples the normal probability distribution.
        //
        //  Discussion:
        //
        //    The normal probability distribution function (PDF) has
        //    a user-specified mean and standard deviation.
        //
        //    This routine can generate a vector of values on one call.  It
        //    has the feature that it should provide the same results
        //    in the same order no matter how we break up the task.
        //
        //    Before calling this routine, the user may call RANDOM_SEED
        //    in order to set the seed of the random number generator.
        //
        //    The Box-Muller method is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of values desired.  If N is negative,
        //    then the code will flush its internal memory; in particular,
        //    if there is a saved value to be used on the next call, it is
        //    instead discarded.  This is useful if the user has reset the
        //    random number seed, for instance.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double NORMAL_VECTOR[N], a sample of the normal PDF.
        //
    {
        double[] x = normal_01_vector(n, ref seed);

        for (int i = 0; i < n; i++)
        {
            x[i] = mu + sigma * x[i];
        }

        return x;
    }

}