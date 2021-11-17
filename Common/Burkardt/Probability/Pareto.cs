using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Pareto
{
    public static double pareto_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_CDF evaluates the Pareto CDF.
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double PARETO_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x < a)
        {
            cdf = 0.0;
        }
        else
        {
            cdf = 1.0 - Math.Pow(a / x, b);
        }

        return cdf;
    }

    public static double pareto_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_CDF_INV inverts the Pareto CDF.
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
        //    0.0 < A,
        //    0.0 < B.
        //
        //    Output, double PARETO_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("PARETO_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                double x = a / Math.Pow(1.0 - cdf, 1.0 / b);

                return x;
            }
        }
    }

    public static bool pareto_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_CHECK checks the parameters of the Pareto CDF.
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
        //    Output, bool PARETO_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("PARETO_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("PARETO_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double pareto_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_MEAN returns the mean of the Pareto PDF.
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
        //    Output, double PARETO_MEAN, the mean of the PDF.
        //
    {
        double mean;

        switch (b)
        {
            case <= 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("PARETO_MEAN - Fatal error!");
                Console.WriteLine("  For B <= 1, the mean does not exist.");
                mean = 0.0;
                return mean;
            default:
                mean = b * a / (b - 1.0);

                return mean;
        }
    }

    public static double pareto_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_PDF evaluates the Pareto PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = B * A^B / X^(B+1).
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
        //    A <= X
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < A.
        //    0.0 < B.
        //
        //    Output, double PARETO_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x < a)
        {
            pdf = 0.0;
        }
        else
        {
            pdf = b * Math.Pow(a, b) / Math.Pow(x, b + 1.0);
        }

        return pdf;
    }

    public static double pareto_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_SAMPLE samples the Pareto PDF.
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
        //    0.0 < A.
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double PARETO_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = pareto_cdf_inv(cdf, a, b);

        return x;
    }

    public static double pareto_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARETO_VARIANCE returns the variance of the Pareto PDF.
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
        //    Output, double PARETO_VARIANCE, the variance of the PDF.
        //
    {
        double variance;

        switch (b)
        {
            case <= 2.0:
                Console.WriteLine(" ");
                Console.WriteLine("PARETO_VARIANCE - Warning!");
                Console.WriteLine("  For B <= 2, the variance does not exist.");
                variance = 0.0;
                return variance;
            default:
                variance = a * a * b / (Math.Pow(b - 1.0, 2) * (b - 2.0));

                return variance;
        }
    }
}