using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Power
{
    public static double power_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_CDF evaluates the Power CDF.
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
        //    0.0 < A, 0.0 < B,
        //
        //    Output, double POWER_CDF, the value of the CDF.
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
                if (x <= b)
                {
                    cdf = Math.Pow(x / b, a);
                }
                else
                {
                    cdf = 1.0;
                }

                break;
            }
        }

        return cdf;
    }

    public static double power_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_CDF_INV inverts the Power CDF.
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
        //    0.0 < A, 0.0 < B.
        //
        //    Output, double POWER_CDF_INV, the argument of the CDF.
        //
    {
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("POWER_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = 0.0;
                break;
            case < 1.0:
                x = b * Math.Exp(Math.Log(cdf) / a);
                break;
            default:
                x = b;
                break;
        }

        return x;
    }

    public static bool power_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_CHECK checks the parameter of the Power PDF.
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
        //    0.0 < A, 0.0 < B.
        //
        //    Output, bool POWER_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POWER_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("POWER_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double power_mean(double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_MEAN returns the mean of the Power PDF.
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
        //    0.0 < A, 0.0 < B.
        //
        //    Output, double POWER_MEAN, the mean of the PDF.
        //
    {
        double mean = a * b / (a + 1.0);

        return mean;
    }

    public static double power_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_PDF evaluates the Power PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = (A/B) * (X/B)^(A-1)
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
        //    Daniel Zwillinger and Stephen Kokoska,
        //    CRC Standard Probability and Statistics Tables and Formulae,
        //    Chapman and Hall/CRC, 2000, pages 152-153.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X <= B.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A, 0.0 < B.
        //
        //    Output, double POWER_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x < 0.0 || b < x)
        {
            pdf = 0.0;
        }
        else
        {
            pdf = a / b * Math.Pow(x / b, a - 1.0);
        }

        return pdf;
    }

    public static double power_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_SAMPLE samples the Power PDF.
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
        //    0.0 < A, 0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double POWER_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = power_cdf_inv(cdf, a, b);

        return x;
    }

    public static double power_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_VARIANCE returns the variance of the Power PDF.
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
        //    0.0 < A, 0.0 < B.
        //
        //    Output, double POWER_VARIANCE, the variance of the PDF.
        //
    {
        double variance = b * b * a / ((a + 1.0) * (a + 1.0) * (a + 2.0));

        return variance;
    }
}