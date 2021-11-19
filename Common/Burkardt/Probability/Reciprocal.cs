using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Reciprocal
{
    public static double reciprocal_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_CDF evaluates the Reciprocal CDF.
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
        //    0.0 < A <= B.
        //
        //    Output, double RECIPROCAL_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            <= 0.0 => 0.0,
            > 0.0 => Math.Log(a / x) / Math.Log(a / b),
            _ => 0
        };

        return cdf;
    }

    public static double reciprocal_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_CDF_INV inverts the Reciprocal CDF.
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
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A <= B.
        //
        //    Output, double RECIPROCAL_CDF_INV, the corresponding argument of the CDF.
        //
    {
        double x = 0;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("RECIPROCAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = 0.0;
                break;
            case > 0.0:
                x = Math.Pow(b, cdf) / Math.Pow(a, cdf - 1.0);
                break;
        }

        return x;
    }

    public static bool reciprocal_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_CHECK checks the parameters of the Reciprocal CDF.
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
        //    0.0 < A <= B.
        //
        //    Output, bool RECIPROCAL_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("RECIPROCAL_CHECK - Warning!");
                Console.WriteLine("  A <= 0.0");
                return false;
        }

        if (!(b < a))
        {
            return true;
        }

        Console.WriteLine(" ");
        Console.WriteLine("RECIPROCAL_CHECK - Warning!");
        Console.WriteLine("  B < A");
        return false;

    }

    public static double reciprocal_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_MEAN returns the mean of the Reciprocal PDF.
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
        //    0.0 < A <= B.
        //
        //    Output, double RECIPROCAL_MEAN, the mean of the PDF.
        //
    {
        double mean = (a - b) / Math.Log(a / b);

        return mean;
    }

    public static double reciprocal_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_PDF evaluates the Reciprocal PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = 1.0 / ( X * LOG ( B / A ) )
        //    for 0.0 <= X
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
        //    0.0 < A <= B.
        //
        //    Output, double RECIPROCAL_PDF, the value of the PDF.
        //
    {
        double pdf = x switch
        {
            <= 0.0 => 0.0,
            > 0.0 => 1.0 / (x * Math.Log(b / a)),
            _ => 0
        };

        return pdf;
    }

    public static double reciprocal_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_SAMPLE samples the Reciprocal PDF.
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
        //    0.0 < A <= B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double RECIPROCAL_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = Math.Pow(b, cdf) / Math.Pow(a, cdf - 1.0);

        return x;
    }

    public static double reciprocal_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RECIPROCAL_VARIANCE returns the variance of the Reciprocal PDF.
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
        //    0.0 < A <= B.
        //
        //    Output, double RECIPROCAL_VARIANCE, the variance of the PDF.
        //
    {
        double d = Math.Log(a / b);

        double variance = (a - b) * (a * (d - 2.0) + b * (d + 2.0))
                          / (2.0 * d * d);

        return variance;
    }
}