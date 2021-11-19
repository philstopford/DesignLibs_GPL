using System;

namespace Burkardt.Probability;

public static class Pearson
{
    public static bool pearson_05_check(double a, double b, double c)

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_CHECK checks the parameters of the Pearson 5 PDF.
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
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, bool PEARSON_05_CHECK, is true if the parameters are legal.
//
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("PEARSON_05_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("PEARSON_05_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }
//****************************************************************************80

    public static double pearson_05_mean(double a, double b, double c)

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_MEAN evaluates the mean of the Pearson 5 PDF.
//
//  Discussion:
//
//    The mean is undefined for B <= 1.
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
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double PEARSON_05_MEAN, the mean of the PDF.
//
    {
        double mean;

        switch (b)
        {
            case <= 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("PEARSON_05_MEAN - Warning!");
                Console.WriteLine("  MEAN undefined for B <= 1.");
                mean = c;
                return mean;
            default:
                mean = c + a / (b - 1.0);

                return mean;
        }
    }
//****************************************************************************80

    public static double pearson_05_pdf(double x, double a, double b, double c)

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_PDF evaluates the Pearson 5 PDF.
//
//  Discussion:
//
//    PDF(A,B;X) = A^B * ( X - C )^(-B-1)
//      * exp ( - A / ( X - C ) ) / Gamma ( B )
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
//    C < X
//
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Output, double PDF, the value of the PDF.
//
    {
        double pdf;

        if (x <= c)
        {
            pdf = 0.0;
        }
        else
        {
            pdf = Math.Pow(a, b) * Math.Pow(x - c, -b - 1.0)
                                 * Math.Exp(-a / (x - c)) / Helpers.Gamma(b);
        }

        return pdf;
    }
//****************************************************************************80

    public static double pearson_05_sample(double a, double b, double c, ref int seed)

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_SAMPLE samples the Pearson 5 PDF.
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
//    Input, double A, B, C, the parameters of the PDF.
//    0.0 < A, 0.0 < B.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double PEARSON_05_SAMPLE, a sample of the PDF.
//
    {
        const double a2 = 0.0;
        double c2 = 1.0 / a;

        double x2 = Gamma.gamma_sample(a2, b, c2, ref seed);

        double x = c + 1.0 / x2;

        return x;
    }
}