using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Fisk
{
    public static double fisk_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_CDF evaluates the Fisk CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double FISK_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x <= a)
        {
            cdf = 0.0;
        }
        else
        {
            cdf = 1.0 / (1.0 + Math.Pow(b / (x - a), c));
        }

        return cdf;
    }

    public static double fisk_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_CDF_INV inverts the Fisk CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double FISK_CDF_INV, the corresponding argument of the CDF.
        //
    {
        const double r8_huge = 1.0E+30;
        double x = 0;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("FISK_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case <= 0.0:
                x = a;
                break;
            case < 1.0:
                x = a + b * Math.Pow(cdf / (1.0 - cdf), 1.0 / c);
                break;
            case >= 1.0:
                x = r8_huge;
                break;
        }

        return x;
    }
//****************************************************************************80

    public static bool fisk_check(double a, double b, double c)

//****************************************************************************80
//
//  Purpose:
//
//    FISK_CHECK checks the parameters of the Fisk PDF.
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
//    Output, bool FISK_CHECK, is true if the parameters are legal.
//
    {
        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("FISK_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
        }

        switch (c)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("FISK_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            default:
                return true;
        }
    }

    public static double fisk_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_MEAN returns the mean of the Fisk PDF.
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
        //    Output, double FISK_MEAN, the mean of the PDF.
        //
    {
            

        switch (c)
        {
            case <= 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("FISK_MEAN - Fatal error!");
                Console.WriteLine("  No mean defined for C <= 1.0");
                return 1;
            default:
            {
                double mean = a + Math.PI * (b / c) * typeMethods.r8_csc(Math.PI / c);

                return mean;
            }
        }
    }

    public static double fisk_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_PDF evaluates the Fisk PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) =
        //      ( C / B ) * ( ( X - A ) / B )^( C - 1 ) /
        //      ( 1 + ( ( X - A ) / B )^C )^2
        //
        //    The Fisk PDF is also known as the Log Logistic PDF.
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
        //    A <= X
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double FISK_PDF, the value of the PDF.
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

            pdf = c / b * Math.Pow(y, c - 1.0) / Math.Pow(1.0 + Math.Pow(y, c), 2);
        }

        return pdf;
    }

    public static double fisk_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_SAMPLE samples the Fisk PDF.
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double FISK_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = fisk_cdf_inv(cdf, a, b, c);

        return x;
    }

    public static double fisk_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FISK_VARIANCE returns the variance of the Fisk PDF.
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
        //    Output, double FISK_VARIANCE, the variance of the PDF.
        //
    {
        double g;
            
        double variance;

        switch (c)
        {
            case <= 2.0:
                Console.WriteLine(" ");
                Console.WriteLine("FISK_VARIANCE - Fatal error!");
                Console.WriteLine("  No variance defined for C <= 2.0");
                return 1;
        }

        g = Math.PI / c;

        variance = b * b * (2.0 * g * typeMethods.r8_csc(2.0 * g)
                            - Math.Pow(g * typeMethods.r8_csc(g), 2));

        return variance;
    }

}