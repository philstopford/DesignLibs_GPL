using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Burr
    {
        public static double burr_cdf(double x, double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_CDF evaluates the Burr CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double BURR_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x <= a)
            {
                cdf = 0.0;
            }
            else
            {
                double y = (x - a) / b;

                cdf = 1.0 - 1.0 / Math.Pow(1.0 + Math.Pow(y, c), d);
            }

            return cdf;
        }

        public static double burr_cdf_inv(double cdf, double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_CDF_INV inverts the Burr CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2016
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
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double BURR_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BURR_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            }

            double y = Math.Pow(Math.Pow(1.0 / (1.0 - cdf), 1.0 / d) - 1.0, 1.0 / c);

            double x = a + b * y;

            return x;
        }

        public static bool burr_check(double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_CHECK checks the parameters of the Burr CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, bool BURR_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BURR_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            if (c <= 0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BURR_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            }

            return true;
        }

        public static double burr_mean(double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_MEAN returns the mean of the Burr PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double BURR_MEAN, the mean of the PDF.
        //
        {
            double ymean = d * typeMethods.r8_beta(d - 1.0 / c, 1.0 + 1.0 / c);

            double mean = a + b * ymean;

            return mean;
        }

        public static double burr_pdf(double x, double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_PDF evaluates the Burr PDF.
        //
        //  Discussion:
        //
        //    Y = ( X - A ) / B;
        //
        //    PDF(X)(A,B,C,D) = ( C * D / B ) * Y ^ ( C - 1 ) / ( 1 + Y ^ C ) ^ ( D + 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    M E Johnson,
        //    Multivariate Statistical Simulation,
        //    Wiley, New York, 1987.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    A <= X
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double BURR_PDF, the value of the PDF.
        //
        {
            double pdf;

            if (x <= a)
            {
                pdf = 0.0;
            }
            else
            {
                double y = (x - a) / b;

                pdf = (c * d / b) * Math.Pow(y, c - 1.0)
                      / Math.Pow(1.0 + Math.Pow(y, c), d + 1.0);
            }

            return pdf;
        }

        public static double burr_sample(double a, double b, double c, double d, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_SAMPLE samples the Burr PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double BURR_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = burr_cdf_inv(cdf, a, b, c, d);

            return x;
        }

        public static double burr_variance(double a, double b, double c, double d)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BURR_VARIANCE returns the variance of the Burr PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, D, the parameters of the PDF.
        //    0 < B,
        //    0 < C.
        //
        //    Output, double BURR_VARIANCE, the variance of the PDF.
        //
        {
            const double r8_huge = 1.0E+30;
            double variance;

            if (c <= 2.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BURR_VARIANCE - Warning!");
                Console.WriteLine("  Variance undefined for C <= 2.");
                variance = r8_huge;
            }
            else
            {
                double mu1 = b * d * typeMethods.r8_beta((c * d - 1.0) / c, (c + 1.0) / c);
                double mu2 = b * b * d * typeMethods.r8_beta((c * d - 2.0) / c, (c + 2.0) / c);
                variance = -mu1 * mu1 + mu2;
            }

            return variance;
        }
    }
}