using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Arcsin
    {
        static double arcsin_cdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_CDF evaluates the Arcsin CDF.
        //
        //  Modified:
        //
        //    20 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Output, double ARCSIN_CDF, the value of the CDF.
        //
        {
            double cdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= -a)
            {
                cdf = 0.0;
            }
            else if (x < a)
            {
                cdf = 0.5 + Math.Asin(x / a) / r8_pi;
            }
            else
            {
                cdf = 1.0;
            }

            return cdf;
        }

        static double arcsin_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_CDF_INV inverts the Arcsin CDF.
        //
        //  Modified:
        //
        //    20 March 2004
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
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Output, double ARCSIN_CDF_INV, the corresponding argument.
        //
        {
            const double r8_pi = 3.14159265358979323;
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("ARCSIN_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            x = a * Math.Sin(r8_pi * (cdf - 0.5));

            return x;
        }

        static bool arcsin_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_CHECK checks the parameter of the Arcsin CDF.
        //
        //  Modified:
        //
        //    27 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the PDF.
        //    0.0 < A.
        //
        //    Output, bool ARCSIN_CHECK, is TRUE if the data is acceptable.
        {
            if (a <= 0.0)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        static double arcsin_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_MEAN returns the mean of the Arcsin PDF.
        //
        //
        //  Modified:
        //
        //    20 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = 0.0;

            return mean;
        }

        static double arcsin_pdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_PDF evaluates the Arcsin PDF.
        //
        //  Discussion:
        //
        //    The LOGISTIC EQUATION has the form:
        //
        //      X(N+1) = 4.0 * LAMBDA * ( 1.0 - X(N) ).
        //
        //    where 0 < LAMBDA <= 1.  This nonlinear difference equation maps
        //    the unit interval into itself, and is a simple example of a system
        //    exhibiting chaotic behavior.  Ulam and von Neumann studied the
        //    logistic equation with LAMBDA = 1, and showed that iterates of the
        //    function generated a sequence of pseudorandom numbers with
        //    the Arcsin probability density function.
        //
        //    The derived sequence
        //
        //      Y(N) = ( 2 / PI ) * Arcsin ( SQRT ( X(N) ) )
        //
        //    is a pseudorandom sequence with the uniform probability density
        //    function on [0,1].  For certain starting values, such as X(0) = 0, 0.75,
        //    or 1.0, the sequence degenerates into a constant sequence, and for
        //    values very near these, the sequence takes a while before becoming
        //    chaotic.
        //
        //  Formula:
        //
        //    PDF(X) = 1 / ( PI * Sqrt ( A * A - X * X ) )
        //
        //  Reference:
        //
        //    Daniel Zwillinger, Stephen Kokoska,
        //    CRC Standard Probability and Statistics Tables and Formulae,
        //    Chapman and Hall/CRC, 2000, pages 114-115.
        //
        //  Modified:
        //
        //    20 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 < X < 1.0.
        //
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Output, double ARCSIN_PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (a <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("ARCSIN_PDF - Fatal error!");
                Console.WriteLine("  Parameter A <= 0.");
                return 1.0;
            }

            if (x <= -a || a <= x)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = 1.0 / (r8_pi * Math.Sqrt(a * a - x * x));
            }

            return pdf;
        }

        static double arcsin_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_SAMPLE samples the Arcsin PDF.
        //
        //  Modified:
        //
        //    20 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double ARCSIN_SAMPLE, a sample of the PDF.
        //
        {
            double cdf;
            double x;

            cdf = UniformRNG.r8_uniform_01(ref seed);

            x = arcsin_cdf_inv(cdf, a);

            return x;
        }

        static double arcsin_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARCSIN_VARIANCE returns the variance of the Arcsin PDF.
        //
        //  Modified:
        //
        //    20 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the parameter of the CDF.
        //    A must be positive.
        //
        //    Output, double ARCSIN_VARIANCE, the variance of the PDF.
        //
        {
            double variance = a * a / 2.0;

            return variance;
        }
    }
}