using System;
using Burkardt.Types;

namespace Burkardt.Probability
{
    public static class Student
    {
        public static double student_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_CDF evaluates the central Student T CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Output, double STUDENT_CDF, the value of the CDF.
        //
        {
            double cdf;

            double y = (x - a) / b;

            double a2 = 0.5 * c;
            double b2 = 0.5;
            double c2 = c / (c + y * y);

            if (y <= 0.0)
            {
                cdf = 0.5 * Beta.beta_inc(a2, b2, c2);
            }
            else
            {
                cdf = 1.0 - 0.5 * Beta.beta_inc(a2, b2, c2);
            }

            return cdf;
        }

        public static void student_cdf_values(ref int n_data, ref double c, ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_CDF_VALUES returns some values of the Student CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = StudentTDistribution [ c ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
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
        //    Output, double &C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 13;

            double[] c_vec =
            {
                1.0, 2.0, 3.0, 4.0,
                5.0, 2.0, 5.0, 2.0,
                5.0, 2.0, 3.0, 4.0,
                5.0
            };

            double[] fx_vec =
            {
                0.6000231200328521,
                0.6001080279134390,
                0.6001150934648930,
                0.6000995134721354,
                0.5999341989834830,
                0.7498859393137811,
                0.7500879487671045,
                0.9500004222186464,
                0.9499969138365968,
                0.9900012348724744,
                0.9900017619355059,
                0.9900004567580596,
                0.9900007637471291
            };

            double[] x_vec =
            {
                0.325,
                0.289,
                0.277,
                0.271,
                0.267,
                0.816,
                0.727,
                2.920,
                2.015,
                6.965,
                4.541,
                3.747,
                3.365
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                c = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                c = c_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static bool student_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_CHECK checks the parameter of the central Student T CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        {
            if (b == 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("STUDENT_CHECK - Warning!");
                Console.WriteLine("  B must be nonzero.");
                return false;
            }

            if (c <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("STUDENT_CHECK - Warning!");
                Console.WriteLine("  C must be greater than 0.");
                return false;
            }

            return true;
        }

        public static double student_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_MEAN returns the mean of the central Student T PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Output, double STUDENT_MEAN, the mean of the PDF.
        //
        {
            double mean = a;

            return mean;
        }

        public static double student_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_PDF evaluates the central Student T PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = Gamma ( (C+1)/2 ) /
        //      ( Gamma ( C / 2 ) * Sqrt ( PI * C )
        //      * ( 1 + ((X-A)/B)^2/C )^(C + 1/2 ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Output, double STUDENT_PDF, the value of the PDF.
        //
        {
            

            double y = (x - a) / b;

            double pdf = Helpers.Gamma(0.5 * (c + 1.0)) / (Math.Sqrt(Math.PI * c)
                                                    * Helpers.Gamma(0.5 * c)
                                                    * Math.Sqrt(Math.Pow((1.0 + y * y / c), (2 * c + 1.0))));

            return pdf;
        }
        public static double student_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_SAMPLE samples the central Student T PDF.
        //
        //  Discussion:
        //
        //    For the sampling algorithm, it is necessary that 2 < C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double STUDENT_SAMPLE, a sample of the PDF.
        //
        {
            if (c <= 2.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("STUDENT_SAMPLE - Fatal error!");
                Console.WriteLine("  Sampling fails for C <= 2.");
                return (1);
            }

            double a2 = 0.0;
            double b2 = c / (c - 2.0);

            double x2 = Normal.normal_sample(a2, b2, ref seed);

            double x3 = Chi.chi_square_sample(c, ref seed);
            x3 = x3 * c / (c - 2.0);

            double x = a + b * x2 * Math.Sqrt(c) / x3;

            return x;
        }

        public static double student_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_VARIANCE returns the variance of the central Student T PDF.
        //
        //  Discussion:
        //
        //    The variance is not defined unless 2 < C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, shape parameters of the PDF,
        //    used to transform the argument X to a shifted and scaled
        //    value Y = ( X - A ) / B.  It is required that B be nonzero.
        //    For the standard distribution, A = 0 and B = 1.
        //
        //    Input, double C, is usually called the number of
        //    degrees of freedom of the distribution.  C is typically an
        //    integer, but that is not essential.  It is required that
        //    C be strictly positive.
        //
        //    Output, double STUDENT_VARIANCE, the variance of the PDF.
        //
        {
            if (c <= 2.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("STUDENT_VARIANCE - Fatal error!");
                Console.WriteLine("  Variance not defined for C <= 2.");
                return (1);
            }

            double variance = b * b * c / (c - 2.0);

            return variance;
        }

        public static double student_noncentral_cdf(double x, int idf, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_NONCENTRAL_CDF evaluates the noncentral Student T CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by B E Cooper.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    BE Cooper,
        //    The Integral of the Noncentral T-Distribution,
        //    Algorithm AS 5,
        //    Applied Statistics,
        //    Volume 17, 1968, page 193.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, int IDF, the number of degrees of freedom.
        //
        //    Input, double D, the noncentrality parameter.
        //
        //    Output, double STUDENT_NONCENTRAL_CDF, the value of the CDF.
        //
        {
            double a;
            int a_max = 100;
            double b;
            double cdf;
            double cdf2;
            double drb;
            double emin = 12.5;
            

            double f = (double) idf;

            if (idf == 1)
            {
                a = x / Math.Sqrt(f);
                b = f / (f + x * x);
                drb = d * Math.Sqrt(b);

                cdf2 = Normal.normal_01_cdf(drb);
                cdf = 1.0 - cdf2 + 2.0 * Owen.tfn(drb, a);
            }
            else if (idf <= a_max)
            {
                a = x / Math.Sqrt(f);
                b = f / (f + x * x);
                drb = d * Math.Sqrt(b);
                double sum2 = 0.0;

                double fmkm2 = 0.0;
                if (Math.Abs(drb) < emin)
                {
                    cdf2 = Normal.normal_01_cdf(a * drb);
                    fmkm2 = a * Math.Sqrt(b) * Math.Exp(-0.5 * drb * drb) * cdf2
                            / Math.Sqrt(2.0 * Math.PI);
                }

                double fmkm1 = b * d * a * fmkm2;
                if (Math.Abs(d) < emin)
                {
                    fmkm1 = fmkm1 + 0.5 * b * a * Math.Exp(-0.5 * d * d) / Math.PI;
                }

                if (idf % 2 == 0)
                {
                    sum2 = fmkm2;
                }
                else
                {
                    sum2 = fmkm1;
                }

                double ak = 1.0;

                int k;
                for (k = 2; k <= idf - 2; k = k + 2)
                {
                    double fk = (double) (k);

                    fmkm2 = b * (d * a * ak * fmkm1 + fmkm2) * (fk - 1.0) / fk;

                    ak = 1.0 / (ak * (fk - 1.0));
                    fmkm1 = b * (d * a * ak * fmkm2 + fmkm1) * fk / (fk + 1.0);

                    if (idf % 2 == 0)
                    {
                        sum2 = sum2 + fmkm2;
                    }
                    else
                    {
                        sum2 = sum2 + fmkm1;
                    }

                    ak = 1.0 / (ak * fk);

                }

                if (idf % 2 == 0)
                {
                    cdf2 = Normal.normal_01_cdf(d);
                    cdf = 1.0 - cdf2 + sum2 * Math.Sqrt(2.0 * Math.PI);
                }
                else
                {
                    cdf2 = Normal.normal_01_cdf(drb);
                    cdf = 1.0 - cdf2 + 2.0 * (sum2 + Owen.tfn(drb, a));
                }
            }
            //
            //  Normal approximation.
            //
            else
            {
                a = Math.Sqrt(0.5 * f) * Math.Exp(Helpers.LogGamma(0.5 * (f - 1.0))
                                             - Helpers.LogGamma(0.5 * f)) * d;

                double temp = (x - a) / Math.Sqrt(f * (1.0 + d * d) / (f - 2.0) - a * a);

                cdf2 = Normal.normal_01_cdf(temp);
                cdf = cdf2;
            }

            return cdf;
        }

        public static void student_noncentral_cdf_values(ref int n_data, ref int df, ref double lambda,
        ref double x, ref double fx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = NoncentralStudentTDistribution [ df, lambda ]
        //      CDF [ dist, x ]
        //
        //    Mathematica seems to have some difficulty computing this function
        //    to the desired number of digits.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2004
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
        //    Output, int &DF, double &LAMBDA, the parameters of the
        //    function.
        //
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 30;

            int[] df_vec =
            {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                15, 20, 25,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10
            };

            double[] fx_vec =
            {
                0.8975836176504333E+00,
                0.9522670169E+00,
                0.9711655571887813E+00,
                0.8231218864E+00,
                0.9049021510E+00,
                0.9363471834E+00,
                0.7301025986E+00,
                0.8335594263E+00,
                0.8774010255E+00,
                0.5248571617E+00,
                0.6293856597E+00,
                0.6800271741E+00,
                0.20590131975E+00,
                0.2112148916E+00,
                0.2074730718E+00,
                0.9981130072E+00,
                0.9994873850E+00,
                0.9998391562E+00,
                0.168610566972E+00,
                0.16967950985E+00,
                0.1701041003E+00,
                0.9247683363E+00,
                0.7483139269E+00,
                0.4659802096E+00,
                0.9761872541E+00,
                0.8979689357E+00,
                0.7181904627E+00,
                0.9923658945E+00,
                0.9610341649E+00,
                0.8688007350E+00
            };

            double[] lambda_vec =
            {
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                4.0E+00,
                4.0E+00,
                4.0E+00,
                7.0E+00,
                7.0E+00,
                7.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00
            };

            double[] x_vec =
            {
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                3.00E+00,
                15.00E+00,
                15.00E+00,
                15.00E+00,
                0.05E+00,
                0.05E+00,
                0.05E+00,
                4.00E+00,
                4.00E+00,
                4.00E+00,
                5.00E+00,
                5.00E+00,
                5.00E+00,
                6.00E+00,
                6.00E+00,
                6.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                df = 0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                df = df_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}