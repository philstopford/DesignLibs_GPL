using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Lerch
    {
        public static double lerch(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LERCH estimates the Lerch transcendent function.
        //
        //  Discussion:
        //
        //    The Lerch transcendent function is defined as:
        //
        //      LERCH ( A, B, C ) = Sum ( 0 <= K < Infinity ) A^K / ( C + K )^B
        //
        //    excluding any term with ( C + K ) = 0.
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
        //  Reference:
        //
        //    Eric Weisstein, editor,
        //    CRC Concise Encylopedia of Mathematics,
        //    CRC Press, 1998.
        //
        //  Thanks:
        //
        //    Oscar van Vlijmen
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the function.
        //
        //    Output, double LERCH, an approximation to the Lerch
        //    transcendent function.
        //
        {
            double sum2 = 0.0;
            int k = 0;
            double a_k = 1.0;

            for (;;)
            {
                double sum2_old = sum2;

                if (c + (double) (k) == 0.0)
                {
                    k = k + 1;
                    a_k = a_k * a;
                    continue;
                }

                sum2 = sum2 + a_k / Math.Pow(c + (double) (k), b);

                if (sum2 <= sum2_old)
                {
                    break;
                }

                k = k + 1;
                a_k = a_k * a;
            }

            return sum2;
        }

        public static double levy_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVY_CDF evaluates the Levy CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the CDF.
        //    Normally, A <= X.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0 < B.
        //
        //    Output, double LEVY_CDF, the value of the PDF.
        //
        {
            double cdf;

            if (b <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVY_CDF - Fatal error!");
                Console.WriteLine("  B <= 0.0.");
                return (1);
            }

            if (x <= a)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = 1.0 - typeMethods.r8_error_f(Math.Sqrt(b / (2.0 * (x - a))));
            }

            return cdf;
        }

        public static double levy_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVY_CDF_INV inverts the Levy CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2007
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
        //    0 < B.
        //
        //    Output, double LEVY_CDF_INV, the corresponding argument.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("LEVY_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            if (b <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVY_CDF_INV - Fatal error!");
                Console.WriteLine("  B <= 0.0.");
                return (1);
            }

            double cdf1 = 1.0 - 0.5 * cdf;
            double x1 = Normal.normal_01_cdf_inv(cdf1);
            double x = a + b / (x1 * x1);

            return x;
        }

        public static double levy_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVY_PDF evaluates the Levy PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = sqrt ( B / ( 2 * PI ) )
        //               * exp ( - B / ( 2 * ( X - A ) )
        //               / ( X - A )^(3/2)
        //
        //    for A <= X.
        //
        //    Note that the Levy PDF does not have a finite mean or variance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    Normally, A <= X.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0 < B.
        //
        //    Output, double LEVY_PDF, the value of the PDF.
        //
        {
            double pdf;
            

            if (b <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("LEVY_PDF - Fatal error!");
                Console.WriteLine("  B <= 0.0.");
                return (1);
            }

            if (x <= a)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = Math.Sqrt(b / (2.0 * Math.PI))
                      * Math.Exp(-b / (2.0 * (x - a)))
                      / Math.Sqrt(Math.Pow(x - a, 3));
            }

            return pdf;
        }

        public static double levy_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVY_SAMPLE samples the Levy PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2007
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
        //    Output, double LEVY_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = levy_cdf_inv(cdf, a, b);

            return x;
        }
    }
}