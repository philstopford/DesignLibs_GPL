using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Anglit
    {
        public static double anglit_cdf(double x)
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
            

            if (x < -0.25 * Math.PI)
            {
                cdf = 0.0;
            }
            else if (x < 0.25 * Math.PI)
            {
                cdf = 0.5 - 0.5 * Math.Cos(2.0 * x + Math.PI / 2.0);
            }
            else
            {
                cdf = 1.0;
            }

            return cdf;
        }

        public static double anglit_cdf_inv(double cdf)
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
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("ANGLIT_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            double x = 0.5 * (Math.Acos(1.0 - 2.0 * cdf) - Math.PI / 2.0);

            return x;
        }

        public static double anglit_mean()
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

        public static double anglit_pdf(double x)
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
            

            if (x <= -0.25 * Math.PI || 0.25 * Math.PI <= x)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = Math.Sin(2.0 * x + 0.25 * Math.PI);
            }

            return pdf;
        }

        public static double anglit_sample(ref int seed)
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
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = anglit_cdf_inv(cdf);

            return x;
        }

        public static double anglit_variance()
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
            double variance = 0.0625 * Math.PI * Math.PI - 0.5;

            return variance;
        }

    }
}
