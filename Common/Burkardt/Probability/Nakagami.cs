using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Nakagami
{
    public static double nakagami_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_CDF evaluates the Nakagami CDF.
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
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B
        //    0.0 < C.
        //
        //    Output, double NAKAGAMI_CDF, the value of the CDF.
        //
    {
        double cdf = 0;

        switch (x)
        {
            case <= 0.0:
                cdf = 0.0;
                break;
            case > 0.0:
            {
                double y = (x - a) / b;
                double x2 = c * y * y;
                double p2 = c;

                cdf = typeMethods.r8_gamma_inc(p2, x2);
                break;
            }
        }

        return cdf;
    }

    public static double nakagami_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_CDF_INV inverts the Nakagami CDF.
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
        //    02 August 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B.
        //    0.0 < C.
        //
        //    Output, double NAKAGAMI_CDF_INV, the corresponding argument of the CDF.
        //
    {
        double cdf2;
        int it_max = 100;
        const double r8_huge = 1.0E+30;
        double tol = 0.000000001;

        double x = 0.0;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("NAKAGAMI_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = c * a * a;
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
            cdf2 = nakagami_cdf(x2, a, b, c);

            if (cdf < cdf2)
            {
                break;
            }

            x2 = a + 2.0 * (x2 - a);
        }

        /*
        Now use bisection.
        */
        int it = 0;

        for (;;)
        {
            it += 1;

            double x3 = 0.5 * (x1 + x2);
            double cdf3 = nakagami_cdf(x3, a, b, c);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine("");
                Console.WriteLine("NAKAGAMI_CDF_INV - Warning!");
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

    public static bool nakagami_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_CHECK checks the parameters of the Nakagami PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, bool NAKAGAMI_CHECK, is true if the parameters are legal.
        //
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("NAKAGAMI_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
        }

        switch (c)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("NAKAGAMI_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double nakagami_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_MEAN returns the mean of the Nakagami PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B
        //    0.0 < C
        //
        //    Output, double NAKAGAMI_MEAN, the mean of the PDF.
        //
    {
        double mean;

        mean = a + b * Helpers.Gamma(c + 0.5) / (Math.Sqrt(c) * Helpers.Gamma(c));

        return mean;
    }

    public static double nakagami_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_PDF evaluates the Nakagami PDF.
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
        //    Input, double X, the argument of the PDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B
        //    0.0 < C.
        //
        //    Output, double NAKAGAMI_PDF, the value of the PDF.
        //
    {
        double pdf = 0;
        double y;

        switch (x)
        {
            case <= 0.0:
                pdf = 0.0;
                break;
            case > 0.0:
                y = (x - a) / b;

                pdf = 2.0 * Math.Pow(c, c) / (b * Helpers.Gamma(c))
                      * Math.Pow(y, 2.0 * c - 1.0)
                      * Math.Exp(-c * y * y);
                break;
        }

        return pdf;
    }

    public static double nakagami_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NAKAGAMI_VARIANCE returns the variance of the Nakagami PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B
        //    0.0 < C
        //
        //    Output, double NAKAGAMI_VARIANCE, the variance of the PDF.
        //
    {
        double t1 = Helpers.Gamma(c + 0.5);
        double t2 = Helpers.Gamma(c);

        double variance = b * b * (1.0 - t1 * t1 / (c * t2 * t2));

        return variance;
    }

}