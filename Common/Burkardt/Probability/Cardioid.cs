using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Cardioid
    {
        static double cardioid_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_CDF evaluates the Cardioid CDF.
        //
        //  Discussion:
        //
        //    The angle X is assumed to lie between A - PI and A + PI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
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
        //    0.0 <= B <= 0.5.
        //
        //    Output, double CDF, the value of the PDF.
        //
        {
            double cdf;
            const double r8_pi = 3.14159265358979323;

            if (x <= a - r8_pi)
            {
                cdf = 0.0;
            }
            else if (x < a + r8_pi)
            {
                cdf = (r8_pi + x - a + 2.0 * b * Math.Sin(x - a)) / (2.0 * r8_pi);
            }
            else
            {
                cdf = 1.0;
            }

            return cdf;
        }

        static double cardioid_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_CDF_INV inverts the Cardioid CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0 <= CDF <= 1.
        //
        //    Input, double A, B, the parameters.
        //    0.0 <= B <= 0.5.
        //
        //    Output, double CARDIOD_CDF_INV, the argument with the given CDF.
        //    A - PI <= X <= A + PI.
        //
        {
            const double r8_pi = 3.14159265358979323;
            double tol = 0.000001;
            double x;

            if (cdf <= 0.0)
            {
                x = a - r8_pi;
            }
            else if (cdf < 1.0)
            {
                x = a;

                int it = 0;

                for (;;)
                {
                    double fx = cdf - (r8_pi + x - a + 2.0 * b * Math.Sin(x - a)) / (2.0 * r8_pi);

                    if (Math.Abs(fx) < tol)
                    {
                        break;
                    }

                    if (10 < it)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("CARDIOID_CDF_INV - Fatal error!");
                        Console.WriteLine("  Too many iterations.");
                        return 1;
                    }

                    double fp = -(1.0 + 2.0 * b * Math.Cos(x - a)) / (2.0 * r8_pi);

                    x = x - fx / fp;
                    x = Math.Max(x, a - r8_pi);
                    x = Math.Min(x, a + r8_pi);

                    it = it + 1;
                }
            }
            else
            {
                x = a + r8_pi;
            }

            return x;
        }

        static bool cardioid_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_CHECK checks the parameters of the Cardioid CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    -0.5 <= B <= 0.5.
        //
        //    Output, bool CARDIOID_CHECK, is true if the parameters are legal.
        //
        {
            bool value = true;

            if (b < -0.5 || 0.5 < b)
            {
                Console.WriteLine("");
                Console.WriteLine("CARDIOID_CHECK - Warning!");
                Console.WriteLine("  B < -0.5 or 0.5 < B.");
                value = false;
            }

            return value;
        }

        static double cardioid_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_MEAN returns the mean of the Cardioid PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= B <= 0.5.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = a;

            return mean;
        }

        static double cardioid_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_PDF evaluates the Cardioid PDF.
        //
        //  Discussion:
        //
        //    The cardioid PDF can be thought of as being applied to points on
        //    a circle.  Compare this distribution with the "Cosine PDF".
        //
        //    PDF(A,B;X) = ( 1 / ( 2 * PI ) ) * ( 1 + 2 * B * COS ( X - A ) )
        //    for 0 <= B <= 1/2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    N I Fisher,
        //    Statistical Analysis of Circular Data,
        //    Cambridge, 1993.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= B <= 0.5.
        //
        //    Output, double CARDIOID_PDF, the value of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double pdf = (1.0 + 2.0 * b * Math.Cos(x - a)) / (2.0 * r8_pi);

            return pdf;
        }

        static double cardioid_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_SAMPLE samples the Cardioid PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= B <= 0.5.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CARDIOD_SAMPLE, a sample of the PDF.
        //    A - PI <= X <= A + PI.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = cardioid_cdf_inv(cdf, a, b);

            return x;
        }

        static double cardioid_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDIOID_VARIANCE returns the variance of the Cardioid PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= B <= 0.5.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            double variance = 0.0;

            return variance;
        }
    }
}