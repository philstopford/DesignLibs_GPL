using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Half
    {
        public static double half_normal_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_CDF evaluates the Half Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Output, double HALF_NORMAL_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= a)
            {
                cdf = 0.0;
            }
            else
            {
                double cdf2 = Normal.normal_cdf(x, a, b);
                cdf = 2.0 * cdf2 - 1.0;
            }

            return cdf;
        }

        public static double half_normal_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_CDF_INV inverts the Half Normal CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Output, double HALF_NORMAL_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("HALF_NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double cdf2 = 0.5 * (cdf + 1.0);

            double x = Normal.normal_cdf_inv(cdf2, a, b);

            return x;
        }

        public static bool half_normal_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_CHECK checks the parameters of the Half Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
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
        //    Output, bool HALF_NORMAL_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("HALF_NORMAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            return true;
        }
        
        public static double half_normal_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_MEAN returns the mean of the Half Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Output, double HALF_NORMAL_MEAN, the mean of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double mean = a + b * Math.Sqrt(2.0 / r8_pi);

            return mean;
        }

        public static double half_normal_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_PDF evaluates the Half Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) =
        //      sqrt ( 2 / PI ) * ( 1 / B ) * exp ( - 0.5 * ( ( X - A ) / B )^2 )
        //
        //    for A <= X
        //
        //    The Half Normal PDF is a special case of both the Chi PDF and the
        //    Folded Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double HALF_NORMAL_PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= a)
            {
                pdf = 0.0;
            }
            else
            {
                double y = (x - a) / b;

                pdf = Math.Sqrt(2.0 / r8_pi) * (1.0 / b) * Math.Exp(-0.5 * y * y);
            }

            return pdf;
        }

        public static double half_normal_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_SAMPLE samples the Half Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Output, double HALF_NORMAL_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = half_normal_cdf_inv(cdf, a, b);

            return x;
        }

        public static double half_normal_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HALF_NORMAL_VARIANCE returns the variance of the Half Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2004
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
        //    Output, double HALF_NORMAL_VARIANCE, the variance of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double variance = b * b * (1.0 - 2.0 / r8_pi);

            return variance;
        }
    }
}