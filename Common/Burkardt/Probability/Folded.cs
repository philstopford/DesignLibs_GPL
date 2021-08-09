using System;
using Burkardt.CDFLib;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Folded
    {
        public static double folded_normal_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_CDF evaluates the Folded Normal CDF.
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
        //    0.0 <= X.
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, double FOLDED_NORMAL_CDF, the value of the CDF.
        //
        {
            double cdf;

            if (x < 0.0)
            {
                cdf = 0.0;
            }
            else
            {
                double x1 = (x - a) / b;
                double cdf1 = Normal.normal_01_cdf(x1);
                double x2 = (-x - a) / b;
                double cdf2 = Normal.normal_01_cdf(x2);
                cdf = cdf1 - cdf2;
            }

            return cdf;
        }

        public static double folded_normal_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.
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
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, double FOLDED_NORMAL_CDF_INV, the argument of the CDF.
        //    0.0 <= X.
        //
        {
            int it_max = 100;
            const double r8_huge = 1.0E+30;
            double tol = 0.0001;
            double x;
            double x1;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("FOLDED_NORMAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            if (cdf == 0.0)
            {
                x = 0.0;
                return x;
            }
            else if (1.0 == cdf)
            {
                x = r8_huge;
                return x;
            }

            //
            //  Find X1, for which the value of CDF will be too small.
            //
            if (0.0 <= a)
            {
                x1 = CDF.normal_cdf_inv(cdf, a, b);
            }
            else
            {
                x1 = CDF.normal_cdf_inv(cdf, -a, b);
            }

            x1 = Math.Max(x1, 0.0);
            double cdf1 = folded_normal_cdf(x1, a, b);
            //
            //  Find X2, for which the value of CDF will be too big.
            //
            double cdf2 = (1.0 - cdf) / 2.0;

            double xa = CDF.normal_cdf_inv(cdf2, a, b);
            double xb = CDF.normal_cdf_inv(cdf2, -a, b);
            double x2 = Math.Max(Math.Abs(xa), Math.Abs(xb));
            cdf2 = folded_normal_cdf(x2, a, b);
            //
            //  Now use bisection.
            //
            int it = 0;

            for (;;)
            {
                it = it + 1;

                double x3 = 0.5 * (x1 + x2);
                double cdf3 = folded_normal_cdf(x3, a, b);

                if (Math.Abs(cdf3 - cdf) < tol)
                {
                    x = x3;
                    break;
                }

                if (it_max < it)
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("FOLDED_NORMAL_CDF_INV - Fatal error!");
                    Console.WriteLine("  Iteration limit IT_MAX = " + it_max + " exceeded.");
                    return (1);
                }

                if ((cdf3 <= cdf && cdf1 <= cdf) || (cdf <= cdf3 && cdf <= cdf1))
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

        public static bool folded_normal_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_CHECK checks the parameters of the Folded Normal CDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, bool FOLDED_NORMAL_CHECK, is true if the parameters are legal.
        //
        {
            if (a < 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("FOLDED_NORMAL_CHECK - Warning!");
                Console.WriteLine("  A < 0.");
                return false;
            }

            if (b <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("FOLDED_NORMAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0.");
                return false;
            }

            return true;
        }

        public static double folded_normal_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_MEAN returns the mean of the Folded Normal PDF.
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
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, double FOLDED_NORMAL_MEAN, the mean of the PDF.
        //
        {
            

            double a2 = a / b;

            double cdf = Normal.normal_01_cdf(a2);

            double mean = b * Math.Sqrt(2.0 / Math.PI) * Math.Exp(-0.5 * a2 * a2)
                          - a * (1.0 - 2.0 * cdf);

            return mean;
        }

        public static double folded_normal_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_PDF evaluates the Folded Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = sqrt ( 2 / PI ) * ( 1 / B ) * COSH ( A * X / B^2 )
        //      * EXP ( - 0.5 * ( X^2 + A^2 ) / B^2 )
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
        //    0.0 <= X
        //
        //    Input, double A, B, the parameters of the PDF.
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, double FOLDED_NORMAL_PDF, the value of the PDF.
        //
        {
            double pdf;
            

            if (x < 0.0)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = Math.Sqrt(2.0 / Math.PI) * (1.0 / b) * Math.Cosh(a * x / b / b)
                      * Math.Exp(-0.5 * (x * x + a * a) / b / b);
            }

            return pdf;
        }

        public static double folded_normal_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_SAMPLE samples the Folded Normal PDF.
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
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double FOLDED_NORMAL_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = folded_normal_cdf_inv(cdf, a, b);

            return x;
        }

        public static double folded_normal_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FOLDED_NORMAL_VARIANCE returns the variance of the Folded Normal PDF.
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
        //    0.0 <= A,
        //    0.0 < B.
        //
        //    Output, double FOLDED_NORMAL_VARIANCE, the variance of the PDF.
        //
        {
            double mean = folded_normal_mean(a, b);

            double variance = a * a + b * b - mean * mean;

            return variance;
        }
    }
}