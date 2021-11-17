using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Gompertz
{
    public static double gompertz_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOMPERTZ_CDF evaluates the Gompertz CDF.
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
        //    Johnson, Kotz, and Balakrishnan,
        //    Continuous Univariate Distributions, Volume 2, second edition,
        //    Wiley, 1994, pages 25-26.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1 < A, 0 < B.
        //
        //    Output, double GOMPERTZ_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            <= 0.0 => 0.0,
            _ => 1.0 - Math.Exp(-b * (Math.Pow(a, x) - 1.0) / Math.Log(a))
        };

        return cdf;
    }

    public static double gompertz_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOMPERTZ_CDF_INV inverts the Gompertz CDF.
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
        //    Johnson, Kotz, and Balakrishnan,
        //    Continuous Univariate Distributions, Volume 2, second edition,
        //    Wiley, 1994, pages 25-26.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1 < A, 0 < B.
        //
        //    Output, double GOMPERTZ_CDF_INV, the corresponding argument.
        //
    {
        const double r8_huge = 1.0E+30;
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("GOMPERTZ_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case < 1.0:
                x = Math.Log(1.0 - Math.Log(1.0 - cdf) * Math.Log(a) / b) / Math.Log(a);
                break;
            default:
                x = r8_huge;
                break;
        }

        return x;
    }

    public static bool gompertz_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOMPERTZ_CHECK checks the parameters of the Gompertz PDF.
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
        //  Reference:
        //
        //    Johnson, Kotz, and Balakrishnan,
        //    Continuous Univariate Distributions, Volume 2, second edition,
        //    Wiley, 1994, pages 25-26.
        //
        //  Parameters:
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1 < A, 0 < B.
        //
        //    Output, bool GOMPERTZ_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("GOMPERTZ_CHECK - Warning!");
                Console.WriteLine("  A <= 1.0!");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("GOMPERTZ_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0!");
                return false;
            default:
                return true;
        }
    }

    public static double gompertz_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOMPERTZ_PDF evaluates the Gompertz PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = B * A**X / exp ( B * ( A**X - 1 ) / log ( A ) )
        //
        //    for
        //
        //      0.0 <= X
        //      1.0 <  A
        //      0.0 <  B
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
        //    Johnson, Kotz, and Balakrishnan,
        //    Continuous Univariate Distributions, Volume 2, second edition,
        //    Wiley, 1994, pages 25-26.
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    1 < A, 0 < B.
        //
        //    Output, double GOMPERTZ_PDF, the value of the PDF.
        //
    {
        double pdf = 0;

        switch (x)
        {
            case < 0.0:
                pdf = 0.0;
                break;
            default:
            {
                pdf = a switch
                {
                    > 1.0 => Math.Exp(Math.Log(b) + x * Math.Log(a) - b / Math.Log(a) * (Math.Pow(a, x) - 1.0)),
                    _ => pdf
                };

                break;
            }
        }

        return pdf;
    }

    public static double gompertz_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOMPERTZ_SAMPLE samples the Gompertz PDF.
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
        //    1 < A, 0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double GOMPERTZ_SAMPLE, a sample of the PDF.
        //
    {
        double cdf;
        double x;

        cdf = UniformRNG.r8_uniform_01(ref seed);

        x = gompertz_cdf_inv(cdf, a, b);

        return x;
    }
}