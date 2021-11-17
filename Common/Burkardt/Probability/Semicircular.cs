using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Semicircular
{
    public static double semicircular_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_CDF evaluates the Semicircular CDF.
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
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, double SEMICIRCULAR_CDF, the value of the CDF.
        //
    {
        double cdf = 0;
            
        double y;

        if (x <= a - b)
        {
            cdf = 0.0;
        }
        else if (x <= a + b)
        {
            y = (x - a) / b;

            cdf = 0.5 + (y * Math.Sqrt(1.0 - y * y) + Math.Asin(y)) / Math.PI;
        }
        else if (a + b < x)
        {
            cdf = 1.0;
        }

        return cdf;
    }

    public static double semicircular_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used on the interval [ A - B, A + B ].
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
        //    0.0 < B.
        //
        //    Output, double SEMICIRCULAR_CDF_INV, the corresponding argument
        //    of the CDF.
        //
    {
        double cdf1;
        double cdf3;
        int it;
        int it_max = 100;
        double tol = 0.0001;
        double x;
        double x1;
        double x2;
        double x3;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("SEMICIRCULAR_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = a - b;
                return x;
            case 1.0:
                x = a + b;
                return x;
        }

        x1 = a - b;
        cdf1 = 0.0;

        x2 = a + b;
        //
        //  Now use bisection.
        //
        it = 0;

        for (;;)
        {
            it += 1;

            x3 = 0.5 * (x1 + x2);
            cdf3 = semicircular_cdf(x3, a, b);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("SEMICIRCULAR_CDF_INV - Fatal error!");
                Console.WriteLine("  Iteration limit exceeded.");
                return 1;
            }

            if (cdf3 <= cdf && cdf1 <= cdf || cdf <= cdf3 && cdf <= cdf1)
            {
                x1 = x3;
                cdf1 = cdf3;
            }
            else
            {
                x2 = x3;
            }

        }

        return x;
    }

    public static bool semicircular_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_CHECK checks the parameters of the Semicircular CDF.
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
        //    Input, double A, B, the parameter of the PDF.
        //    0.0 < B.
        //
        //    Output, bool SEMICIRCULAR_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("SEMICIRCULAR_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
            default:
                return true;
        }
    }

    public static double semicircular_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_MEAN returns the mean of the Semicircular PDF.
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
        //    0.0 < B.
        //
        //    Output, double SEMICIRCULAR_MEAN, the mean of the PDF.
        //
    {
        double mean = a;

        return mean;
    }

    public static double semicircular_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_PDF evaluates the Semicircular PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = ( 2 / ( B * PI ) ) * SQRT ( 1 - ( ( X - A ) / B )^2 )
        //    for A - B <= X <= A + B
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
        //    0.0 < B.
        //
        //    Output, double SEMICIRCULAR_PDF, the value of the PDF.
        //
    {
        double pdf = 0;
            
        double y;

        if (x < a - b)
        {
            pdf = 0.0;
        }
        else if (x <= a + b)
        {
            y = (x - a) / b;

            pdf = 2.0 / (b * Math.PI) * Math.Sqrt(1.0 - y * y);
        }
        else if (a + b < x)
        {
            pdf = 0.0;
        }

        return pdf;
    }

    public static double semicircular_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_SAMPLE samples the Semicircular PDF.
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
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double SEMICIRCULAR_SAMPLE, a sample of the PDF.
        //
    {
        double radius = UniformRNG.r8_uniform_01(ref seed);
        radius = b * Math.Sqrt(radius);
        double angle = Math.PI * UniformRNG.r8_uniform_01(ref seed);
        double x = a + radius * Math.Cos(angle);

        return x;
    }

    public static double semicircular_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEMICIRCULAR_VARIANCE returns the variance of the Semicircular PDF.
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
        //    0.0 < B.
        //
        //    Output, double SEMICIRCULAR_VARIANCE, the variance of the PDF.
        //
    {
        double variance = b * b / 4.0;

        return variance;
    }
}