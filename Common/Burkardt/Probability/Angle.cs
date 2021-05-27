using System;

namespace Burkardt.Probability
{
    public static class Angle
    {
        static double angle_cdf(double x, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_CDF evaluates the Angle CDF.
        //
        //  Modified:
        //
        //    19 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
        //    Wiley, 1986.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, int N, the spatial dimension.
        //    N must be at least 2.
        //
        //    Output, real CDF, the value of the CDF.
        //
        {
            double cdf;
            const double r8_pi = 3.14159265358979323;
            double zero = 0.0;

            if (n < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("ANGLE_CDF - Fatal error!");
                Console.WriteLine("  N must be at least 2.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1.0;
            }

            if (x < 0.0)
            {
                cdf = 0.0;
            }
            else if (r8_pi < x)
            {
                cdf = 1.0;
            }
            else if (n == 2)
            {
                cdf = x / r8_pi;
            }
            else
            {
                cdf = Misc.sin_power_int(zero, x, n - 2)
                      * Helpers.Gamma((double) (n) / 2.0)
                      / (Math.Sqrt(r8_pi) * Helpers.Gamma((double) (n - 1) / 2.0));
            }

            return cdf;
        }

        static double angle_mean(int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_MEAN returns the mean of the Angle PDF.
        //
        //  Modified:
        //
        //    02 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //    N must be at least 2.
        //
        //    Output, double ANGLE_MEAN, the mean of the PDF.
        //
        {
            double mean;
            const double r8_pi = 3.14159265358979323;

            mean = r8_pi / 2.0;

            return mean;
        }

        static double angle_pdf(double x, int n)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_PDF evaluates the Angle PDF.
        //
        //  Discussion:
        //
        //    X is an angle between 0 and PI, corresponding to the angle
        //    made in an N-dimensional space, between a fixed line passing
        //    through the origin, and an arbitrary line that also passes
        //    through the origin, which is specified by a choosing any point
        //    on the N-dimensional sphere with uniform probability.
        //
        //  Formula:
        //
        //    PDF(X) = ( sin ( X ) )^(N-2) * Gamma ( N / 2 )
        //             / ( sqrt ( PI ) * Gamma ( ( N - 1 ) / 2 ) )
        //
        //    PDF(X) = 1 / PI if N = 2.
        //
        //  Modified:
        //
        //    02 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation and Sensitivity of Queueing Networks,
        //    Wiley, 1986.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, int N, the spatial dimension.
        //    N must be at least 2.
        //
        //    Output, real PDF, the value of the PDF.
        //
        {
            double pdf;
            const double r8_pi = 3.14159265358979323;

            if (n < 2)
            {
                Console.WriteLine("");
                Console.WriteLine("ANGLE_PDF - Fatal error!");
                Console.WriteLine("  N must be at least 2.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1.0;
            }

            if (x < 0.0 || r8_pi < x)
            {
                pdf = 0.0;
            }
            else if (n == 2)
            {
                pdf = 1.0 / r8_pi;
            }
            else
            {
                pdf = Math.Pow((Math.Sin(x)), (n - 2))
                      * Helpers.Gamma((double) (n) / 2.0)
                      / (Math.Sqrt(r8_pi) * Helpers.Gamma((double) (n - 1) / 2.0));
            }
            
            return pdf;
        }
    }
}
