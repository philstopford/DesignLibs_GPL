﻿using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Chebyshevi
    {
        static double chebyshev1_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_CDF evaluates the Chebyshev1 CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double CHEBYSHEV1_CDF, the value of the CDF.
        //
        {
            double cdf;
            const double r8_pi = 3.14159265358979323;

            if (x < -1.0)
            {
                cdf = 0.0;
            }
            else if (1.0 < x)
            {
                cdf = 1.0;
            }
            else
            {
                cdf = 0.5 + Math.Asin(x) / r8_pi;
            }

            return cdf;
        }

        static double chebyshev1_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_CDF_INV inverts the Chebyshev1 CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
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
        //    Output, double CHEBYSHEV1_CDF_INV, the corresponding argument.
        //
        {
            const double r8_pi = 3.14159265358979323;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("CHEBYSHEV1_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            }

            double x = Math.Sin(r8_pi * (cdf - 0.5));

            return x;
        }

        static double chebyshev1_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_MEAN returns the mean of the Chebyshev1 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CHEBYSHEV1_MEAN, the mean of the PDF.
        //
        {
            double mean = 0.0;

            return mean;
        }

        static double chebyshev1_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_PDF evaluates the Chebyshev1 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= -1.0 || 1.0 <= x)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = 1.0 / r8_pi / Math.Sqrt(1.0 - x * x);
            }

            return pdf;
        }

        static double chebyshev1_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_SAMPLE samples the Chebyshev1 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CHEBYSHEV1_SAMPLE, a random sample.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double value = chebyshev1_cdf_inv(cdf);

            return value;
        }

        static double chebyshev1_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV1_VARIANCE returns the variance of the Chebyshev1 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CHEBYSHEV1_VARIANCE, the variance of the PDF.
        //
        {
            double variance = 0.5;

            return variance;
        }
    }
}