using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Zipf
    {
        public static double zipf_cdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_CDF evaluates the Zipf CDF.
        //
        //  Discussion:
        //
        //    Simple summation is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the PDF.
        //    1 <= N
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double c;
            double cdf;
            double pdf;
            int y;

            if (x < 1)
            {
                cdf = 0.0;
            }
            else
            {
                c = typeMethods.r8_zeta(a);

                cdf = 0.0;
                for (y = 1; y <= x; y++)
                {
                    pdf = 1.0 / Math.Pow(y, a) / c;
                    cdf = cdf + pdf;
                }

            }

            return cdf;
        }

        public static int zipf_cdf_inv(double a, double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_CDF_INV inverts the Zipf CDF.
        //
        //  Discussion:
        //
        //    Simple summation is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, int X, the argument such that
        //    CDF(X-1) < CDF <= CDF(X)
        //    1 <= X <= 1000
        //
        {
            double c;
            double cdf2;
            double pdf;
            int x;
            int y;

            if (cdf <= 0.0)
            {
                x = 1;
            }
            else
            {
                c = typeMethods.r8_zeta(a);
                cdf2 = 0.0;

                x = 1000;

                for (y = 1; y <= 1000; y++)
                {
                    pdf = (1.0 / Math.Pow(y, a)) / c;
                    cdf2 = cdf2 + pdf;
                    if (cdf <= cdf2)
                    {
                        x = y;
                        break;
                    }
                }
            }

            return x;
        }

        public static bool zipf_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_CHECK checks the parameter of the Zipf PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, bool ZIPF_CHECK, is true if the parameters are legal.
        //
        {
            if (a <= 1.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("ZIPF_CHECK - Warning!");
                Console.WriteLine("  A <= 1.");
                return false;
            }

            return true;
        }

        public static double zipf_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_MEAN returns the mean of the Zipf PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, double ZIPF_MEAN, the mean of the PDF.
        //    The mean is only defined for 2 < A.
        //
        {
            double mean;

            if (a <= 2.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("ZIPF_MEAN - Fatal error!");
                Console.WriteLine("  No mean defined for A <= 2.");
                return 1.0;
            }

            mean = typeMethods.r8_zeta(a - 1.0) / typeMethods.r8_zeta(a);

            return mean;
        }

        public static double zipf_pdf(int x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_PDF evaluates the Zipf PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = ( 1 / X^A ) / C
        //
        //    where the normalizing constant is chosen so that
        //
        //    C = Sum ( 1 <= I < Infinity ) 1 / I^A.
        //
        //    From observation, the frequency of different words in long
        //    sequences of text seems to follow the Zipf PDF, with
        //    parameter A slightly greater than 1.  The Zipf PDF is sometimes
        //    known as the "discrete Pareto" PDF.
        //
        //    Lotka's law is a version of the Zipf PDF in which A is 2 or approximately
        //    2.  Lotka's law describes the frequency of publications by authors in a
        //    given field, and estimates that the number of authors with X papers is
        //    about 1/X^A of the number of authors with 1 paper.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alfred Lotka,
        //    The frequency distribution of scientific productivity,
        //    Journal of the Washington Academy of Sciences,
        //    Volume 16, Number 12, 1926, pages 317-324.
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the PDF.
        //    1 <= N
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < 1)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = 1.0 / Math.Pow(x, a) / typeMethods.r8_zeta(a);
            }

            return pdf;
        }
        
        public static double[] zipf_probability ( int n, double p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_PROBABILITY sets up a Zipf probability vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    George Zipf,
        //    The Psychobiology of Language,
        //    1935.
        //
        //  Parameters:
        //
        //    Input, unsigned int N, indicates the size of X.
        //
        //    Input, double P, the Zipf parameter.
        //    1.0 < P.
        //
        //    Output, double ZIPF_PROBABILITY[N+2], contains in X[1] through X[N] the
        //    probabilities of outcomes 1 through N.
        //
        {
            int i;
            double[] x;

            x = new double[n+2];

            x[0] = 0.0;
            for ( i = 1; i <= n; i++ )
            {
                x[i] = Math.Pow ( i, - p );
            }
            x[n+1] = 0.0;

            Helpers.normalize ( n, ref x );

            return x;
        }

        public static int zipf_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_SAMPLE samples the Zipf PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Luc Devroye,
        //    Non-Uniform Random Variate Generation,
        //    Springer Verlag, 1986, pages 550-551.
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int ZIPF_SAMPLE, a sample of the PDF.
        //
        {
            double b;
            double t;
            double u;
            double v;
            double w;
            int x;

            b = Math.Pow(2.0, (a - 1.0));

            for (;;)
            {
                u = UniformRNG.r8_uniform_01(ref seed);
                v = UniformRNG.r8_uniform_01(ref seed);
                w = (int) (1.0 / Math.Pow(u, 1.0 / (a - 1.0)));

                t = Math.Pow((w + 1.0) / w, a - 1.0);

                if (v * w * (t - 1.0) * b <= t * (b - 1.0))
                {
                    break;
                }

            }

            x = (int) w;

            return x;
        }

        public static double zipf_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZIPF_VARIANCE returns the variance of the Zipf PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    1.0 < A.
        //
        //    Output, double ZIPF_VARIANCE, the variance of the PDF.
        //    The variance is only defined for 3 < A.
        //
        {
            double mean;
            double variance;

            if (a <= 3.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("ZIPF_VARIANCE - Fatal error!");
                Console.WriteLine("  No variance defined for A <= 3.0.");
                return 1.0;
            }

            mean = zipf_mean(a);

            variance = typeMethods.r8_zeta(a - 2.0) / typeMethods.r8_zeta(a) - mean * mean;

            return variance;
        }
    }
}
