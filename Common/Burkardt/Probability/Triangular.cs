using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Triangular
{
    public static double triangular_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_CDF evaluates the Triangular CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
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
        //    A < B.
        //
        //    Output, double TRIANGULAR_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x <= a)
        {
            cdf = 0.0;
        }
        else if (x <= 0.5 * (a + b))
        {
            cdf = 2.0 * (x * x - 2.0 * a * x + a * a) / Math.Pow(b - a, 2);
        }
        else if (x <= b)
        {
            cdf = 0.5 + (-2.0 * x * x + 4.0 * b * x + 0.5 * a * a
                         - a * b - 1.5 * b * b) / Math.Pow(b - a, 2);
        }
        else
        {
            cdf = 1.0;
        }

        return cdf;
    }

    public static double triangular_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_CDF_INV inverts the Triangular CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
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
        //    A < B.
        //
        //    Output, double TRIANGULAR_CDF_INV, the corresponding argument.
        //
    {
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("TRIANGULAR_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case <= 0.5:
                x = a + 0.5 * (b - a) * Math.Sqrt(2.0 * cdf);
                break;
            default:
                x = b - 0.5 * (b - a) * Math.Sqrt(2.0 * (1.0 - cdf));
                break;
        }

        return x;
    }

    public static bool triangular_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_CHECK checks the parameters of the Triangular CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, bool TRIANGULAR_CHECK, is true if the parameters are legal.
        //
    {
        if (b <= a)
        {
            Console.WriteLine(" ");
            Console.WriteLine("TRIANGULAR_CHECK - Warning!");
            Console.WriteLine("  B <= A.");
            return false;
        }

        return true;
    }

    public static double triangular_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_MEAN returns the mean of the Triangular PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, double TRIANGULAR_MEAN, the mean of the discrete PDF.
        //
    {
        double mean = 0.5 * (a + b);

        return mean;
    }

    public static double triangular_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_PDF evaluates the Triangular PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = 4 * ( X - A ) / ( B - A )^2 for A <= X <= (A+B)/2
        //                = 4 * ( B - X ) / ( B - A )^2 for (A+B)/2 <= X <= B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
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
        //    A < B.
        //
        //    Output, double TRIANGULAR_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x <= a)
        {
            pdf = 0.0;
        }
        else if (x <= 0.5 * (a + b))
        {
            pdf = 4.0 * (x - a) / (b - a) / (b - a);
        }
        else if (x <= b)
        {
            pdf = 4.0 * (b - x) / (b - a) / (b - a);
        }
        else
        {
            pdf = 0.0;
        }

        return pdf;
    }

    public static double triangular_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_SAMPLE samples the Triangular PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double TRIANGULAR_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = triangular_cdf_inv(cdf, a, b);

        return x;
    }

    public static double triangular_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_VARIANCE returns the variance of the Triangular PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, double TRIANGULAR_VARIANCE, the variance of the PDF.
        //
    {
        double variance = (b - a) * (b - a) / 24.0;

        return variance;
    }
}