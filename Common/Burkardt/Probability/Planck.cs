using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Planck
{
    public static bool planck_check(double a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_CHECK checks the parameters of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
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
//    Output, bool PLANCK_CHECK, is true if the parameters are legal.
//
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("PLANCK_CHECK - Warning!");
                Console.WriteLine("  A <= 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("PLANCK_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            default:
                return true;
        }
    }
//****************************************************************************80

    public static double planck_mean(double a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_MEAN returns the mean of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
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
//    Output, double PLANCK_MEAN, the mean of the PDF.
//
    {
        double mean = (b + 1.0) * typeMethods.r8_zeta(b + 2.0) / typeMethods.r8_zeta(b + 1.0);

        return mean;
    }
//****************************************************************************80

    public static double planck_pdf(double x, double a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_PDF evaluates the Planck PDF.
//
//  Discussion:
//
//    The Planck PDF has the form
//
//      PDF(A,B;X) = A^(B+1) * X^B / ( exp ( A * X ) - 1 ) / K
//
//    where K is the normalization constant, and has the value
//
//      K = Gamma ( B + 1 ) * Zeta ( B + 1 ).
//
//    The original Planck distribution governed the frequencies in
//    blackbody radiation at a given temperature T, and has the form
//
//      PDF(A;X) = K * X^3 / ( exp ( A * X ) - 1 )
//
//    where
//
//      K = 15 / PI^4.
//
//    Thus, in terms of the Planck PDF, the original Planck distribution
//    has A = 1, B = 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//    0.0 <= X
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Output, double PLANCK_PDF, the value of the PDF.
//
    {
        double pdf;

        switch (x)
        {
            case < 0.0:
                pdf = 0.0;
                break;
            default:
                double k = Helpers.Gamma(b + 1.0) * typeMethods.r8_zeta(b + 1.0);
                pdf = Math.Pow(a, b + 1.0) * Math.Pow(x, b) / (Math.Exp(a * x) - 1.0) / k;
                break;
        }

        return pdf;
    }
//****************************************************************************80

    public static double planck_sample(double a, double b, ref int seed)

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_SAMPLE samples the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Luc Devroye,
//    Non-Uniform Random Variate Generation,
//    Springer Verlag, 1986, pages 552.
//
//  Parameters:
//
//    Input, double A, B, the parameters of the PDF.
//    0.0 < A,
//    0.0 < B.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double PLANCK_SAMPLE, a sample of the PDF.
//
    {
        const double a2 = 0.0;
        const double b2 = 1.0;
        double c2 = b + 1.0;

        double g = Gamma.gamma_sample(a2, b2, c2, ref seed);

        int z = Zipf.zipf_sample(c2, ref seed);

        double x = g / (a * z);

        return x;
    }
//****************************************************************************80

    public static double planck_variance(double a, double b)

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_VARIANCE returns the variance of the Planck PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2004
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
//    Output, double PLANCK_VARIANCE, the variance of the PDF.
//
    {
        const double variance = 0.0;

        return variance;
    }

}