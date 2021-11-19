using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class InverseGaussian
{
    public static double inverse_gaussian_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF.
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
        //    0.0 < X.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double INVERSE_GAUSSIAN_CDF, the value of the CDF.
        //
    {
        double cdf;

        switch (x)
        {
            case <= 0.0:
                cdf = 0.0;
                break;
            default:
            {
                double x1 = Math.Sqrt(b / x) * (x - a) / a;
                double cdf1 = Normal.normal_01_cdf(x1);

                double x2 = -Math.Sqrt(b / x) * (x + a) / a;
                double cdf2 = Normal.normal_01_cdf(x2);

                cdf = cdf1 + Math.Exp(2.0 * b / a) * cdf2;
                break;
            }
        }

        return cdf;
    }

    public static bool inverse_gaussian_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_CHECK checks the parameters of the Inverse Gaussian CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, bool INVERSE_GAUSSIAN_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("INVERSE_GAUSSIAN_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("INVERSE_GAUSSIAN_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double inverse_gaussian_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_MEAN returns the mean of the Inverse Gaussian PDF.
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double INVERSE_GAUSSIAN_MEAN, the mean of the PDF.
        //
    {
        double mean = a;

        return mean;
    }

    public static double inverse_gaussian_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF.
        //
        //  Discussion:
        //
        //    The Inverse Gaussian PDF is also known as the Wald PDF
        //    and the Inverse Normal PDF.
        //
        //    PDF(A,B;X)
        //      = sqrt ( B / ( 2 * PI * X^3 ) )
        //        * exp ( - B * ( X - A )^2 / ( 2.0 * A^2 * X ) )
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
        //    0.0 < X
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double INVERSE_GAUSSIAN_PDF, the value of the PDF.
        //
    {
        double pdf = x switch
        {
            <= 0.0 => 0.0,
            _ => Math.Sqrt(b / (2.0 * Math.PI * Math.Pow(x, 3))) * Math.Exp(-b * (x - a) * (x - a) / (2.0 * a * a * x))
        };

        return pdf;
    }

    public static double inverse_gaussian_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian PDF.
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double INVERSE_GAUSSIAN_SAMPLE, a sample of the PDF.
        //
    {
        double phi = b / a;
        double z = Normal.normal_01_sample(ref seed);
        double y = z * z;

        double t = 1.0 + 0.5 * (y - Math.Sqrt(4.0 * phi * y + y * y)) / phi;
        double u = UniformRNG.r8_uniform_01(ref seed);

        double x = (u * (1.0 + t)) switch
        {
            <= 1.0 => a * t,
            _ => a / t
        };

        return x;
    }

    public static double inverse_gaussian_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_GAUSSIAN_VARIANCE returns the variance of the Inverse Gaussian PDF.
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double INVERSE_GAUSSIAN_VARIANCE, the variance of the PDF.
        //
    {
        double variance = a * a * a / b;

        return variance;
    }

}