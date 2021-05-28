using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Frechet
    {
        public static double frechet_cdf(double x, double alpha)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_CDF evaluates the Frechet CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the parameter.
        //    It is required that 0.0 < ALPHA.
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double cdf;

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_CDF - Fatal error!");
                Console.WriteLine("  ALPHA <= 0.0.");
                return (1);
            }

            if (x <= 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = Math.Exp(-1.0 / Math.Pow(x, alpha));
            }

            return cdf;
        }

        public static double frechet_cdf_inv(double cdf, double alpha)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_CDF_INV inverts the Frechet CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
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
        //    Input, double ALPHA, the parameter.
        //    It is required that 0.0 < ALPHA.
        //
        //    Output, double FRECHET_CDF_INV, the corresponding argument of the CDF.
        //
        {
            double x;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_CDF_INV - Fatal error!");
                Console.WriteLine("  ALPHA <= 0.0.");
                return (1);
            }

            if (cdf == 0.0)
            {
                x = 0.0;
            }
            else
            {
                x = Math.Pow(-1.0 / Math.Log(cdf), 1.0 / alpha);
            }

            return x;
        }

        public static double frechet_mean(double alpha)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_MEAN returns the mean of the Frechet PDF.
        //
        //  Discussion:
        //
        //    The distribution does not have a mean value unless 1 < ALPHA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the parameter.
        //    It is required that 1.0 < ALPHA.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            if (alpha <= 1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_MEAN - Fatal error!");
                Console.WriteLine("  Mean does not exist if ALPHA <= 1.");
                return (1);
            }

            double mean = Helpers.Gamma((alpha - 1.0) / alpha);

            return mean;
        }

        public static double frechet_pdf(double x, double alpha)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_PDF evaluates the Frechet PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = ALPHA * exp ( -1 / X^ALPHA ) / X^(ALPHA+1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double ALPHA, the parameter.
        //    It is required that 0.0 < ALPHA.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_PDF - Fatal error!");
                Console.WriteLine("  ALPHA <= 0.0.");
                return (1);
            }

            double pdf = alpha * Math.Exp(-1.0 / Math.Pow(x, alpha)) / Math.Pow(x, alpha + 1.0);

            return pdf;
        }

        public static double frechet_sample(double alpha, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_SAMPLE samples the Frechet PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the parameter.
        //    It is required that 0.0 < ALPHA.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double FRECHET_SAMPLE, a sample of the PDF.
        //
        {
            if (alpha <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_SAMPLE - Fatal error!");
                Console.WriteLine("  ALPHA <= 0.0.");
                return (1);
            }

            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = frechet_cdf_inv(cdf, alpha);

            return x;
        }

        public static double frechet_variance(double alpha)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FRECHET_VARIANCE returns the variance of the Frechet PDF.
        //
        //  Discussion:
        //
        //    The PDF does not have a variance unless 2 < ALPHA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ALPHA, the parameter.
        //    It is required that 2.0 < ALPHA.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            if (alpha <= 2.0)
            {
                Console.WriteLine("");
                Console.WriteLine("FRECHET_VARIANCE - Fatal error!");
                Console.WriteLine("  Variance does not exist if ALPHA <= 2.");
                return (1);
            }

            double mean = Helpers.Gamma((alpha - 1.0) / alpha);

            double variance = Helpers.Gamma((alpha - 2.0) / alpha) - mean * mean;

            return variance;
        }
    }
}
