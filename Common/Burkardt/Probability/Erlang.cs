using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Erlang
{
    public static double erlang_cdf(double x, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_CDF evaluates the Erlang CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, double ERLANG_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x < a)
        {
            cdf = 0.0;
        }
        else
        {
            double x2 = (x - a) / b;
            double p2 = c;

            cdf = typeMethods.r8_gamma_inc(p2, x2);
        }

        return cdf;
    }

    public static double erlang_cdf_inv(double cdf, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_CDF_INV inverts the Erlang CDF.
        //
        //  Discussion:
        //
        //    A simple bisection method is used.
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
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, double ERLANG_CDF_INV, the corresponding argument of the CDF.
        //
    {
        double cdf2;
        int it_max = 100;
        const double r8_huge = 1.0E+30;
        double tol = 0.0001;

        double x = 0.0;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("ERLANG_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = a;
                return x;
            case 1.0:
                x = r8_huge;
                return x;
        }

        double x1 = a;
        double cdf1 = 0.0;

        double x2 = a + 1.0;

        for (;;)
        {
            cdf2 = erlang_cdf(x2, a, b, c);

            if (cdf < cdf2)
            {
                break;
            }

            x2 = a + 2.0 * (x2 - a);
        }

        //
        //  Now use bisection.
        //
        int it = 0;

        for (;;)
        {
            it += 1;

            double x3 = 0.5 * (x1 + x2);
            double cdf3 = erlang_cdf(x3, a, b, c);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("ERLANG_CDF_INV - Warning!");
                Console.WriteLine("  Iteration limit exceeded.");
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
                cdf2 = cdf3;
            }
        }

        return x;
    }

    public static bool erlang_check(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_CHECK checks the parameters of the Erlang PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, bool ERLANG_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("ERLANG_CHECK - Warning!");
                Console.WriteLine("  B <= 0.0");
                return false;
        }

        switch (c)
        {
            case <= 0:
                Console.WriteLine(" ");
                Console.WriteLine("ERLANG_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double erlang_mean(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_MEAN returns the mean of the Erlang PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, double ERLANG_MEAN, the mean of the PDF.
        //
    {
        double mean = a + b * c;

        return mean;
    }

    public static double erlang_pdf(double x, double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_PDF evaluates the Erlang PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = ( ( X - A ) / B )^( C - 1 )
        //      / ( B * Gamma ( C ) * EXP ( ( X - A ) / B ) )
        //
        //    for 0 < B, 0 < C integer, A <= X.
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
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, double ERLANG_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x <= a)
        {
            pdf = 0.0;
        }
        else
        {
            double y = (x - a) / b;

            pdf = Math.Pow(y, c - 1) / (b * typeMethods.r8_factorial(c - 1) * Math.Exp(y));
        }

        return pdf;
    }

    public static double erlang_sample(double a, double b, int c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_SAMPLE samples the Erlang PDF.
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
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double ERLANG_SAMPLE, a sample of the PDF.
        //
    {
        double a2 = 0.0;
        double b2 = b;
        double x = a;

        for (int i = 1; i <= c; i++)
        {
            double x2 = Exponential.exponential_sample(a2, b2, ref seed);
            x += x2;
        }

        return x;
    }

    public static double erlang_variance(double a, double b, int c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERLANG_VARIANCE returns the variance of the Erlang PDF.
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
        //    Input, double A, B, int C, the parameters of the PDF.
        //    0.0 < B.
        //    0 < C.
        //
        //    Output, double ERLANG_VARIANCE, the variance of the PDF.
        //
    {
        double variance = b * b * c;

        return variance;
    }
}