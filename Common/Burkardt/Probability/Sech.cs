using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Sech
    {
        public static double sech(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH returns the hyperbolic secant.
        //
        //  Discussion:
        //
        //    SECH ( X ) = 1.0 / COSH ( X ) = 2.0 / ( EXP ( X ) + EXP ( - X ) )
        //
        //    SECH is not a built-in function in FORTRAN, and occasionally it
        //    is handier, or more concise, to be able to refer to it directly
        //    rather than through its definition in terms of the sine function.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double SECH, the hyperbolic secant of X.
        //
        {
            double value = 1.0 / Math.Cosh(x);

            return value;
        }

        public static double sech_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_CDF evaluates the Hyperbolic Secant CDF.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, double SECH_CDF, the value of the CDF.
        //
        {
            double cdf;
            
            double y;

            y = (x - a) / b;

            cdf = 2.0 * Math.Atan(Math.Exp(y)) / Math.PI;

            return cdf;
        }

        public static double sech_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_CDF_INV inverts the Hyperbolic Secant CDF.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double SECH_CDF_INV, the corresponding argument of the CDF.
        //
        {
            const double r8_huge = 1.0E+30;
            
            double x = 0;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("SECH_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return(1);
            }

            if (cdf == 0.0)
            {
                x = -r8_huge;
            }
            else if (cdf < 1.0)
            {
                x = a + b * Math.Log(Math.Tan(0.5 * Math.PI * cdf));
            }
            else if (1.0 == cdf)
            {
                x = r8_huge;
            }

            return x;
        }

        public static bool sech_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_CHECK checks the parameters of the Hyperbolic Secant CDF.
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
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, bool SECH_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("SECH_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
            }

            return true;
        }

        public static double sech_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_MEAN returns the mean of the Hyperbolic Secant PDF.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double SECH_MEAN, the mean of the PDF.
        //
        {
            double mean = a;

            return mean;
        }

        public static double sech_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_PDF evaluates the Hypebolic Secant PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = sech ( ( X - A ) / B ) / ( PI * B )
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double SECH_PDF, the value of the PDF.
        //
        {
            

            double y = (x - a) / b;

            double pdf = sech(y) / (Math.PI * b);

            return pdf;
        }

        public static double sech_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_SAMPLE samples the Hyperbolic Secant PDF.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SECH_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = a + b * Math.Log(Math.Tan(0.5 * Math.PI * cdf));

            return x;
        }

        public static double sech_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECH_VARIANCE returns the variance of the Hyperbolic Secant PDF.
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double SECH_VARIANCE, the variance of the PDF.
        //
        {
            

            double variance = 0.25 * Math.PI * Math.PI * b * b;

            return variance;
        }
    }
}