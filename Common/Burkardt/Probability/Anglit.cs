using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Anglit
    {
        static double anglit_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_CDF evaluates the Anglit CDF.
        //
        //  Modified:
        //
        //    29 December 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double ANGLIT_CDF, the value of the CDF.
        //
        {
            double cdf;
            const double r8_pi = 3.14159265358979323;

            if (x < -0.25 * r8_pi)
            {
                cdf = 0.0;
            }
            else if (x < 0.25 * r8_pi)
            {
                cdf = 0.5 - 0.5 * Math.Cos(2.0 * x + r8_pi / 2.0);
            }
            else
            {
                cdf = 1.0;
            }

            return cdf;
        }

        static double anglit_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_CDF_INV inverts the Anglit CDF.
        //
        //  Modified:
        //
        //    29 December 1999
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
        //    Output, double ANGLIT_CDF_INV, the corresponding argument.
        //
        {
            const double r8_pi = 3.14159265358979323;
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("ANGLIT_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            x = 0.5 * (Math.Acos(1.0 - 2.0 * cdf) - r8_pi / 2.0);

            return x;
        }

        static double anglit_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_MEAN returns the mean of the Anglit PDF.
        //
        //  Modified:
        //
        //    28 December 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double ANGLIT_MEAN, the mean of the PDF.
        //
        {
            return 0.0;
        }

        static double anglit_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_PDF evaluates the Anglit PDF.
        //
        //  Formula:
        //
        //    PDF(X) = SIN ( 2 * X + PI / 2 ) for -PI/4 <= X <= PI/4
        //
        //  Modified:
        //
        //    28 December 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Output, double ANGLIT_PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= -0.25 * r8_pi || 0.25 * r8_pi <= x)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = Math.Sin(2.0 * x + 0.25 * r8_pi);
            }

            return pdf;
        }

        static double anglit_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_SAMPLE samples the Anglit PDF.
        //
        //  Modified:
        //
        //    28 December 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double ANGLIT_SAMPLE, a sample of the PDF.
        //
        {
            double cdf;
            double x;

            cdf = UniformRNG.r8_uniform_01(ref seed);

            x = anglit_cdf_inv(cdf);

            return x;
        }

        static double anglit_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLIT_VARIANCE returns the variance of the Anglit PDF.
        //
        //  Discussion:
        //
        //    Variance =
        //      Integral ( -PI/4 <= X <= PI/4 ) X^2 * SIN ( 2 * X + PI / 2 )
        //
        //    Antiderivative =
        //      0.5 * X * SIN ( 2 * X + PI / 2 )
        //      + ( 0.25 - 0.5 * X^2 ) * COS ( 2 * X + PI / 2 )
        //
        //  Modified:
        //
        //    29 December 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double ANGLIT_VARIANCE, the variance of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;
            double variance;

            variance = 0.0625 * r8_pi * r8_pi - 0.5;

            return variance;
        }

    }
}
