using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Cosine
{
    public static double cosine_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_CDF evaluates the Cosine CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Output, double CDF, the value of the CDF.
        //
    {
        double cdf = 0;

        if (x <= a - Math.PI * b)
        {
            cdf = 0.0;
        }
        else if (x <= a + Math.PI * b)
        {
            double y = (x - a) / b;

            cdf = 0.5 + (y + Math.Sin(y)) / (2.0 * Math.PI);
        }
        else if (a + Math.PI * b < x)
        {
            cdf = 1.0;
        }

        return cdf;
    }

    public static double cosine_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_CDF_INV inverts the Cosine CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used on the interval
        //    [ A - PI * B, A + PI * B ].
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
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 < B.
        //
        //    Output, double COSINE_CDF_INV, the corresponding argument of the CDF.
        //
    {
        const int it_max = 100;
            
        const double tol = 0.0001;
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("COSINE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = a - Math.PI * b;
                return x;
            case 1.0:
                x = a + Math.PI * b;
                return x;
        }

        double x1 = a - Math.PI * b;
        double cdf1 = 0.0;

        double x2 = a + Math.PI * b;
        //
        //  Now use bisection.
        //
        int it = 0;

        for (it = 1; it <= it_max; it++)
        {
            double x3 = 0.5 * (x1 + x2);
            double cdf3 = cosine_cdf(x3, a, b);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                return x;
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

        Console.WriteLine(" ");
        Console.WriteLine("COSINE_CDF_INV - Fatal error!");
        Console.WriteLine("  Iteration limit exceeded.");
        return 1;
    }

    public static bool cosine_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_CHECK checks the parameters of the Cosine CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Output, bool COSINE_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("COSINE_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
            default:
                return true;
        }
    }

    public static double cosine_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_MEAN returns the mean of the Cosine PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Output, double MEAN, the mean of the PDF.
        //
    {
        return a;
    }

    public static double cosine_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_PDF evaluates the Cosine PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = ( 1 / ( 2 * PI * B ) ) * COS ( ( X - A ) / B )
        //    for A - PI * B <= X <= A + PI * B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf = 0;

        if (x < a - Math.PI * b)
        {
            pdf = 0.0;
        }
        else if (x <= a + Math.PI * b)
        {
            double y = (x - a) / b;

            pdf = 1.0 / (2.0 * Math.PI * b) * Math.Cos(y);
        }
        else if (a + Math.PI * b < x)
        {
            pdf = 0.0;
        }

        return pdf;
    }

    public static double cosine_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_SAMPLE samples the Cosine PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
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
        //    Output, double COSINE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = cosine_cdf_inv(cdf, a, b);

        return x;
    }

    public static double cosine_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COSINE_VARIANCE returns the variance of the Cosine PDF.
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
        //    0.0 < B.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
            

        double variance = (Math.PI * Math.PI / 3.0 - 2.0) * b * b;

        return variance;
    }
}