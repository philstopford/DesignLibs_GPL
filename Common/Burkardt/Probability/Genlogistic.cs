using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Genlogistic
    {
        public static double genlogistic_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_CDF evaluates the Generalized Logistic CDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double CDF, the value of the CDF.
        //
        {
            double y = (x - a) / b;

            double cdf = 1.0 / Math.Pow((1.0 + Math.Exp(-y)), c);

            return cdf;
        }

        public static double genlogistic_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_CDF_INV inverts the Generalized Logistic CDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double GENLOGISTIC_CDF_INV, the corresponding argument.
        //
        {
            const double r8_huge = 1.0E+30;
            double x = 0;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GENLOGISTIC_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            if (cdf == 0.0)
            {
                x = -r8_huge;
            }
            else if (cdf < 1.0)
            {
                x = a - b * Math.Log(Math.Pow(cdf, -1.0 / c) - 1.0);
            }
            else if (1.0 == cdf)
            {
                x = r8_huge;
            }

            return x;
        }

        public static bool genlogistic_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_CHECK checks the parameters of the Generalized Logistic CDF.
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
        //    Output, bool GENLOGISTIC_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GENLOGISTIC_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            if (c <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("GENLOGISTIC_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            }

            return true;
        }

        public static double genlogistic_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_MEAN returns the mean of the Generalized Logistic PDF.
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
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double MEAN, the mean of the PDF.
        //
        {
            double mean = a + b * (Misc.euler_constant() + Digamma.digamma(c));

            return mean;
        }

        public static double genlogistic_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_PDF evaluates the Generalized Logistic PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = ( C / B ) * EXP ( ( A - X ) / B ) /
        //      ( ( 1 + EXP ( ( A - X ) / B ) )^(C+1) )
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
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double PDF, the value of the PDF.
        //
        {
            double y = (x - a) / b;

            double pdf = (c / b) * Math.Exp(-y) / Math.Pow(1.0 + Math.Exp(-y), c + 1.0);

            return pdf;
        }

        public static double genlogistic_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_SAMPLE samples the Generalized Logistic PDF.
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
        //    Input, double A, B, C, the parameters of the PDF.
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double GENLOGISTIC_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = genlogistic_cdf_inv(cdf, a, b, c);

            return x;
        }

        public static double genlogistic_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENLOGISTIC_VARIANCE returns the variance of the Generalized Logistic PDF.
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
        //    0.0 < B,
        //    0.0 < C.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double variance = b * b * (r8_pi * r8_pi / 6.0 + trigamma(c));

            return variance;
        }
    }
}