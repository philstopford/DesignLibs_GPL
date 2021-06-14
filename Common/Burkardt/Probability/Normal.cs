using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
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
            double a1 = 0.398942280444;
            double a2 = 0.399903438504;
            double a3 = 5.75885480458;
            double a4 = 29.8213557808;
            double a5 = 2.62433121679;
            double a6 = 48.6959930692;
            double a7 = 5.92885724438;
            double b0 = 0.398942280385;
            double b1 = 3.8052E-08;
            double b2 = 1.00000615302;
            double b3 = 3.98064794E-04;
            double b4 = 1.98615381364;
            double b5 = 0.151679116635;
            double b6 = 5.29330324926;
            double b7 = 4.8385912808;
            double b8 = 15.1508972451;
            double b9 = 0.742380924027;
            double b10 = 30.789933034;
            double b11 = 3.99019417011;
            double cdf;
            double q;
            double y;
            //
            //  |X| <= 1.28.
            //
            if (Math.Abs(x) <= 1.28)
            {
                y = 0.5 * x * x;

                q = 0.5 - Math.Abs(x) * (a1 - a2 * y / (y + a3 - a4 / (y + a5
                                                                         + a6 / (y + a7))));
            //
            //  1.28 < |X| <= 12.7
            //
            }
            else if (Math.Abs(x) <= 12.7)
            {
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
            }
            else
            {
                q = 0.0;
            }

            //
            //  Take account of negative X.
            //
            if (x < 0.0)
            {
                cdf = q;
            }
            else
            {
                cdf = 1.0 - q;
            }

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
            double q;
            double r;
            const double r8_huge = 1.0E+30;
            double split1 = 0.425;
            double split2 = 5.0;
            double value;

            if (p <= 0.0)
            {
                value = -r8_huge;
                return value;
            }

            if (1.0 <= p)
            {
                value = r8_huge;
                return value;
            }

            q = p - 0.5;

            if (Math.Abs(q) <= split1)
            {
                r = const1 - q * q;
                value = q * typeMethods.r8poly_value_horner(7, a, r) / typeMethods.r8poly_value_horner(7, b, r);
            }
            else
            {
                if (q < 0.0)
                {
                    r = p;
                }
                else
                {
                    r = 1.0 - p;
                }

                if (r <= 0.0)
                {
                    value = r8_huge;
                }
                else
                {
                    r = Math.Sqrt(-Math.Log(r));

                    if (r <= split2)
                    {
                        r = r - const2;
                        value = typeMethods.r8poly_value_horner(7, c, r) / typeMethods.r8poly_value_horner(7, d, r);
                    }
                    else
                    {
                        r = r - split2;
                        value = typeMethods.r8poly_value_horner(7, e, r) / typeMethods.r8poly_value_horner(7, f, r);
                    }
                }

                if (q < 0.0)
                {
                    value = -value;
                }

            }

            return value;
        }

        public static void normal_01_cdf_values(ref int n_data, ref double x, ref double fx )
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
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
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
            
            double[] x;

            double[] r1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            double[] r2 = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            x = new double[n];
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
            double variance = 1.0;

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
            int i;
            int m;
            
            double[] r;

            double[] x = new double[n];
            //
            //  Record the range of X we need to fill in.
            //
            int x_lo = 1;
            int x_hi = n;
            //
            //  If we need just one new value, do that here to avoid null arrays.
            //
            if (x_hi - x_lo + 1 == 1)
            {
                r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

                x[x_hi - 1] = Math.Sqrt(-2.0 * Math.Log(r[0])) * Math.Cos(2.0 * Math.PI * r[1]);
            }
            //
            //  If we require an even number of values, that's easy.
            //
            else if ((x_hi - x_lo + 1) % 2 == 0)
            {
                m = (x_hi - x_lo + 1) / 2;

                r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                for (i = 0; i <= 2 * m - 2; i = i + 2)
                {
                    x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                    x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                }
            }
            //
            //  If we require an odd number of values, we generate an even number,
            //  and handle the last pair specially, storing one in X(N).
            //
            else
            {
                x_hi = x_hi - 1;

                m = (x_hi - x_lo + 1) / 2 + 1;

                r = UniformRNG.r8vec_uniform_01_new(2 * m, ref seed);

                for (i = 0; i <= 2 * m - 4; i = i + 2)
                {
                    x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);
                    x[x_lo + i] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Sin(2.0 * Math.PI * r[i + 1]);
                }

                i = 2 * m - 2;

                x[x_lo + i - 1] = Math.Sqrt(-2.0 * Math.Log(r[i])) * Math.Cos(2.0 * Math.PI * r[i + 1]);

            }

            return x;
        }

        public static double normal_cdf(double x, double mu, double sigma)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF evaluates the Normal CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double y = (x - mu) / sigma;

            double cdf = normal_01_cdf(y);

            return cdf;
        }

        public static double normal_cdf_inv(double cdf, double mu, double sigma)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF_INV inverts the Normal CDF.
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
        //  Reference:
        //
        //    Algorithm AS 111,
        //    Applied Statistics,
        //    Volume 26, pages 118-121, 1977.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //    SIGMA should not be zero.
        //
        //    Output, double NORMAL_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            double x2 = normal_01_cdf_inv(cdf);

            double x = mu + sigma * x2;

            return x;
        }

        public static void normal_cdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
            ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_CDF_VALUES returns some values of the Normal CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NormalDistribution [ mu, sigma ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2004
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
        //    Output, double &MU, the mean of the distribution.
        //
        //    Output, double &SIGMA, the variance of the distribution.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 12;

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.9772498680518208E+00,
                0.9999683287581669E+00,
                0.9999999990134124E+00,
                0.6914624612740131E+00,
                0.6305586598182364E+00,
                0.5987063256829237E+00,
                0.5792597094391030E+00,
                0.6914624612740131E+00,
                0.5000000000000000E+00,
                0.3085375387259869E+00,
                0.1586552539314571E+00
            };

            double[] mu_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
            };

            double[] sigma_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01
            };

            double[] x_vec =
            {
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
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
            if (sigma == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_CHECK - Warning!");
                Console.WriteLine("  SIGMA == 0.");
                return false;
            }

            return true;
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

        public static double normal_truncated_ab_cdf(double x, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_CDF evaluates the truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_CDF, the value of the CDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);
            double xi_cdf = normal_01_cdf(xi);

            double cdf = (xi_cdf - alpha_cdf) / (beta_cdf - alpha_cdf);

            return cdf;
        }

        public static double normal_truncated_ab_cdf_inv(double cdf, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_CDF_INV inverts the truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_AB_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);

            double xi_cdf = (beta_cdf - alpha_cdf) * cdf + alpha_cdf;
            double xi = normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_ab_mean(double mu, double s, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_MEAN returns the mean of the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_MEAN, the mean of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);

            double alpha_pdf = normal_01_pdf(alpha);
            double beta_pdf = normal_01_pdf(beta);

            double mean = mu + s * (alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf);

            return mean;
        }

        public static double normal_truncated_ab_pdf(double x, double mu, double s, double a,
            double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_PDF evaluates the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_PDF, the value of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);
            double xi_pdf = normal_01_pdf(xi);

            double pdf = xi_pdf / (beta_cdf - alpha_cdf) / s;

            return pdf;
        }

        public static double normal_truncated_ab_sample(double mu, double s, double a, double b,
            ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_SAMPLE samples the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_AB_SAMPLE, a sample of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);

            double u = UniformRNG.r8_uniform_01(ref seed);
            double xi_cdf = alpha_cdf + u * (beta_cdf - alpha_cdf);
            double xi = normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_ab_variance(double mu, double s, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_AB_VARIANCE returns the variance of the truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //
        //    Output, double NORMAL_TRUNCATED_AB_VARIANCE, the variance of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double beta = (b - mu) / s;

            double alpha_pdf = normal_01_pdf(alpha);
            double beta_pdf = normal_01_pdf(beta);

            double alpha_cdf = normal_01_cdf(alpha);
            double beta_cdf = normal_01_cdf(beta);

            double variance = s * s * (1.0
                                       + (alpha * alpha_pdf - beta * beta_pdf) / (beta_cdf - alpha_cdf)
                                       - Math.Pow((alpha_pdf - beta_pdf) / (beta_cdf - alpha_cdf), 2));

            return variance;
        }

        public static double normal_truncated_a_cdf(double x, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_CDF evaluates the lower truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_CDF, the value of the CDF.
        //
        {
            double alpha = (a - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double xi_cdf = normal_01_cdf(xi);

            double cdf = (xi_cdf - alpha_cdf) / (1.0 - alpha_cdf);

            return cdf;
        }

        public static double normal_truncated_a_cdf_inv(double cdf, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_CDF_INV inverts the lower truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_A_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            double alpha = (a - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);

            double xi_cdf = (1.0 - alpha_cdf) * cdf + alpha_cdf;
            double xi = normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_a_mean(double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_MEAN returns the mean of the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_MEAN, the mean of the PDF.
        //
        {
            double alpha = (a - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);

            double alpha_pdf = normal_01_pdf(alpha);

            double mean = mu + s * alpha_pdf / (1.0 - alpha_cdf);

            return mean;
        }

        public static double normal_truncated_a_pdf(double x, double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_PDF evaluates the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_PDF, the value of the PDF.
        //
        {
            double alpha = (a - mu) / s;
            double xi = (x - mu) / s;

            double alpha_cdf = normal_01_cdf(alpha);
            double xi_pdf = normal_01_pdf(xi);

            double pdf = xi_pdf / (1.0 - alpha_cdf) / s;

            return pdf;
        }

        public static double normal_truncated_a_sample(double mu, double s, double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_SAMPLE samples the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_A_SAMPLE, a sample of the PDF.
        //
        {
            double alpha;
            double alpha_cdf;
            double u;
            double x;
            double xi;
            double xi_cdf;

            alpha = (a - mu) / s;

            alpha_cdf = normal_01_cdf(alpha);

            u = UniformRNG.r8_uniform_01(ref seed);
            xi_cdf = alpha_cdf + u * (1.0 - alpha_cdf);
            xi = normal_01_cdf_inv(xi_cdf);

            x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_a_variance(double mu, double s, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_A_VARIANCE: variance of the lower truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double A, the lower truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_A_VARIANCE, the variance of the PDF.
        //
        {
            double alpha = (a - mu) / s;

            double alpha_pdf = normal_01_pdf(alpha);

            double alpha_cdf = normal_01_cdf(alpha);

            double variance = s * s * (1.0
                                       + alpha * alpha_pdf / (1.0 - alpha_cdf)
                                       - Math.Pow(alpha_pdf / (1.0 - alpha_cdf), 2));

            return variance;
        }

        public static double normal_truncated_b_cdf(double x, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_CDF evaluates the upper truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_CDF, the value of the CDF.
        //
        {
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double beta_cdf = normal_01_cdf(beta);
            double xi_cdf = normal_01_cdf(xi);

            double cdf = xi_cdf / beta_cdf;

            return cdf;
        }

        public static double normal_truncated_b_cdf_inv(double cdf, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_CDF_INV inverts the upper truncated Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("NORMAL_TRUNCATED_B_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            double beta = (b - mu) / s;

            double beta_cdf = normal_01_cdf(beta);

            double xi_cdf = beta_cdf * cdf;
            double xi = normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_b_mean(double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_MEAN returns the mean of the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviatione of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_MEAN, the mean of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_cdf = normal_01_cdf(beta);

            double beta_pdf = normal_01_pdf(beta);

            double mean = mu - s * beta_pdf / beta_cdf;

            return mean;
        }

        public static double normal_truncated_b_pdf(double x, double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_PDF evaluates the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_PDF, the value of the PDF.
        //
        {
            double beta = (b - mu) / s;
            double xi = (x - mu) / s;

            double beta_cdf = normal_01_cdf(beta);
            double xi_pdf = normal_01_pdf(xi);

            double pdf = xi_pdf / beta_cdf / s;

            return pdf;
        }

        public static double normal_truncated_b_sample(double mu, double s, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_SAMPLE samples the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Input/output, int &SEED, a seed for the random number
        //    generator.
        //
        //    Output, double NORMAL_TRUNCATED_B_SAMPLE, a sample of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_cdf = normal_01_cdf(beta);

            double u = UniformRNG.r8_uniform_01(ref seed);
            double xi_cdf = u * beta_cdf;
            double xi = normal_01_cdf_inv(xi_cdf);

            double x = mu + s * xi;

            return x;
        }

        public static double normal_truncated_b_variance(double mu, double s, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMAL_TRUNCATED_B_VARIANCE: variance of the upper truncated Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double MU, S, the mean and standard deviation of the
        //    parent Normal distribution.
        //
        //    Input, double B, the upper truncation limit.
        //
        //    Output, double NORMAL_TRUNCATED_B_VARIANCE, the variance of the PDF.
        //
        {
            double beta = (b - mu) / s;

            double beta_pdf = normal_01_pdf(beta);

            double beta_cdf = normal_01_cdf(beta);

            double variance = s * s * (1.0
                                       - beta * beta_pdf / beta_cdf
                                       - Math.Pow(beta_pdf / beta_cdf, 2));

            return variance;
        }
    }
}