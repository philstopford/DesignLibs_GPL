using System;

namespace Burkardt.Probability;

public static class Angle
{
    public static double angle_cdf(double x, int n)
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
        const double zero = 0.0;

        switch (n)
        {
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("ANGLE_CDF - Fatal error!");
                Console.WriteLine("  N must be at least 2.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1.0;
            default:
                double cdf;
                switch (x)
                {
                    case < 0.0:
                        cdf = 0.0;
                        break;
                    case > Math.PI:
                        cdf = 1.0;
                        break;
                    default:
                    {
                        cdf = n switch
                        {
                            2 => x / Math.PI,
                            _ => Misc.sin_power_int(zero, x, n - 2) * Helpers.Gamma(n / 2.0) /
                                 (Math.Sqrt(Math.PI) * Helpers.Gamma((n - 1) / 2.0))
                        };

                        break;
                    }
                }

                return cdf;
        }
    }

    public static double angle_mean(int n)
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
        const double mean = Math.PI / 2.0;

        return mean;
    }

    public static double angle_pdf(double x, int n)
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

        switch (n)
        {
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("ANGLE_PDF - Fatal error!");
                Console.WriteLine("  N must be at least 2.");
                Console.WriteLine("  The input value of N = " + n + "");
                return 1.0;
        }

        switch (x)
        {
            case < 0.0:
            case > Math.PI:
                pdf = 0.0;
                break;
            default:
            {
                pdf = n switch
                {
                    2 => 1.0 / Math.PI,
                    _ => Math.Pow(Math.Sin(x), n - 2) * Helpers.Gamma(n / 2.0) /
                         (Math.Sqrt(Math.PI) * Helpers.Gamma((n - 1) / 2.0))
                };

                break;
            }
        }
            
        return pdf;
    }
}