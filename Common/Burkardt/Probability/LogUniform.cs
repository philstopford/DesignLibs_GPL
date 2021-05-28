using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class LogUniform
    {
        public static double log_uniform_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_CDF evaluates the Log Uniform CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= a)
            {
                cdf = 0.0;
            }
            else if (x < b)
            {
                cdf = (Math.Log(x) - Math.Log(a)) / (Math.Log(b) - Math.Log(a));
            }
            else
            {
                cdf = 1.0;
            }

            return cdf;
        }

        public static double log_uniform_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_CDF_INV inverts the Log Uniform CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double LOG_UNIFORM_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine( "");
                Console.WriteLine( "LOG_UNIFORM_CDF_INV - Fatal error!");
                Console.WriteLine( "  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double x = a * Math.Exp((Math.Log(b) - Math.Log(a)) * cdf);

            return x;
        }

        public static bool log_uniform_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_CHECK checks the parameters of the Log Uniform CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1.0 < A < B.
        //
        //    Output, bool LOG_UNIFORM_CHECK, is true if the parameters are legal.
        //
        {
            if (a <= 1.0)
            {
                Console.WriteLine( "");
                Console.WriteLine( "LOG_UNIFORM_CHECK - Warning!");
                Console.WriteLine( "  A <= 1.");
                return false;
            }

            if (b <= a)
            {
                Console.WriteLine( "");
                Console.WriteLine( "LOG_UNIFORM_CHECK - Warning!");
                Console.WriteLine( "  B <= A.");
                return false;
            }

            return true;
        }

        public static double log_uniform_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_MEAN returns the mean of the Log Uniform PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1.0 < A < B.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = (b - a) / (Math.Log(b) - Math.Log(a));

            return mean;
        }

        public static double log_uniform_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_PDF evaluates the Log Uniform PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = ( log ( B ) - log ( A ) ) / X  for A <= X <= B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1.0 < A < B.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x < a)
            {
                pdf = 0.0;
            }
            else if (x <= b)
            {
                pdf = 1.0 / (x * (Math.Log(b) - Math.Log(a)));
            }
            else
            {
                pdf = 0.0;
            }

            return pdf;
        }

        public static double log_uniform_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_SAMPLE samples the Log Uniform PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1.0 < A < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double LOG_UNIFORM_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = log_uniform_cdf_inv(cdf, a, b);

            return x;
        }

        public static double log_uniform_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_UNIFORM_VARIANCE returns the variance of the Log Uniform PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1.0 < A < B.
        //
        //    Output, double LOG_UNIFORM_VARIANCE, the variance of the PDF.
        //
        {
            double mean = log_uniform_mean(a, b);

            double variance = ((0.5 * b * b - 2.0 * mean * b + mean * mean * Math.Log(b))
                               - (0.5 * a * a - 2.0 * mean * a + mean * mean * Math.Log(a)))
                              / (Math.Log(b) - Math.Log(a));

            return variance;
        }
    }
}