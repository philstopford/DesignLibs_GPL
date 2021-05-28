using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Bradford
    {
        static double bradford_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_CDF evaluates the Bradford CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
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
        //    A < B,
        //    0.0 < C.
        //
        //    Output, double BRADFORD_CDF, the value of the CDF.
        //
        {
            double cdf = 0;

            if (x <= a)
            {
                cdf = 0.0;
            }
            else if (x <= b)
            {
                cdf = Math.Log(1.0 + c * (x - a) / (b - a)) / Math.Log(c + 1.0);
            }
            else if (b < x)
            {
                cdf = 1.0;
            }

            return cdf;
        }

        static double bradford_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_CDF_INV inverts the Bradford CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
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
        //    A < B,
        //    0.0 < C.
        //
        //    Output, double BRADFORD_CDF_INV, the corresponding argument of the CDF.
        //
        {
            double x = 0;

            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BRADFORD_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1.0;
            }

            if (cdf <= 0.0)
            {
                x = a;
            }
            else if (cdf < 1.0)
            {
                x = a + (b - a) * (Math.Pow((c + 1.0), cdf) - 1.0) / c;
            }
            else if (1.0 <= cdf)
            {
                x = b;
            }

            return x;
        }

        static bool bradford_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_CHECK checks the parameters of the Bradford PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A < B,
        //    0.0 < C.
        //
        //    Output, bool BRADFORD_CHECK, is true if the parameters are legal.
        //
        {
            if (b <= a)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BRADFORD_CHECK - Warning!");
                Console.WriteLine("  B <= A.");
                return false;
            }

            if (c <= 0.0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("BRADFORD_CHECK - Warning!");
                Console.WriteLine("  C <= 0.");
                return false;
            }

            return true;
        }

        static double bradford_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_MEAN returns the mean of the Bradford PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A < B,
        //    0.0 < C.
        //
        //    Output, double BRADFORD_MEAN, the mean of the PDF.
        //
        {
            double mean = (c * (b - a) + Math.Log(c + 1.0) * (a * (c + 1.0) - b))
                          / (c * Math.Log(c + 1.0));

            return mean;
        }

        static double bradford_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_PDF evaluates the Bradford PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B,C;X) = C / ( ( C * ( X - A ) + B - A ) * log ( C + 1 ) )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
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
        //    A < B,
        //    0.0 < C.
        //
        //    Output, double BRADFORD_PDF, the value of the PDF.
        //
        {
            double pdf = 0;

            if (x <= a)
            {
                pdf = 0.0;
            }
            else if (x <= b)
            {
                pdf = c / ((c * (x - a) + b - a) * Math.Log(c + 1.0));
            }
            else if (b < x)
            {
                pdf = 0.0;
            }

            return pdf;
        }

        static double bradford_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_SAMPLE samples the Bradford PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A < B,
        //    0.0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double BRADFORD_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = a + (b - a) * (Math.Pow((c + 1.0), cdf) - 1.0) / c;

            return x;
        }

        static double bradford_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BRADFORD_VARIANCE returns the variance of the Bradford PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A < B,
        //    0.0 < C.
        //
        //    Output, double BRADFORD_VARIANCE, the variance of the PDF.
        //
        {
            double variance = (b - a) * (b - a) *
                              (c * (Math.Log(c + 1.0) - 2.0) + 2.0 * Math.Log(c + 1.0))
                              / (2.0 * c * Math.Pow((Math.Log(c + 1.0)), 2));

            return variance;
        }
    }
}