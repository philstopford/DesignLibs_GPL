using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Gumbel
    {
        public static double gumbel_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_CDF evaluates the Gumbel CDF.
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
        //    Output, double GUMBEL_CDF, the value of the CDF.
        //
        {
            double cdf = Math.Exp(-Math.Exp(-x));

            return cdf;
        }

        public static double gumbel_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_CDF_INV inverts the Gumbel CDF.
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
        //    Output, double GUMBEL_CDF_INV, the corresponding argument of the CDF.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GUMBEL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double x = -Math.Log(-Math.Log(cdf));

            return x;
        }

        public static double gumbel_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_MEAN returns the mean of the Gumbel PDF.
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
        //    Output, double GUMBEL_MEAN, the mean of the PDF.
        //
        {
            double mean = Misc.euler_constant();

            return mean;
        }

        public static double gumbel_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_PDF evaluates the Gumbel PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = EXP ( - X - EXP ( - X  ) ).
        //
        //    GUMBEL_PDF(X) = EXTREME_PDF(0,1;X)
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
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Output, double GUMBEL_PDF, the value of the PDF.
        //
        {
            double pdf = Math.Exp(-x - Math.Exp(-x));

            return pdf;
        }

        public static double gumbel_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_SAMPLE samples the Gumbel PDF.
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double GUMBEL_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = gumbel_cdf_inv(cdf);

            return x;
        }

        public static double gumbel_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUMBEL_VARIANCE returns the variance of the Gumbel PDF.
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
        //    Output, double GUMBEL_VARIANCE, the variance of the PDF.
        //
        {
            

            double variance = Math.PI * Math.PI / 6.0;

            return variance;
        }
    }
}