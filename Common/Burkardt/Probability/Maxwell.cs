using System;
using Burkardt.Types;

namespace Burkardt.Probability;

public static class Maxwell
{
    public static double maxwell_cdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_CDF evaluates the Maxwell CDF.
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
        //    0.0 <= X
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double MAXWELL_CDF, the value of the CDF.
        //
    {
        double cdf;

        switch (x)
        {
            case <= 0.0:
                cdf = 0.0;
                break;
            default:
            {
                double x2 = x / a;
                double p2 = 1.5;

                cdf = typeMethods.r8_gamma_inc(p2, x2);
                break;
            }
        }

        return cdf;
    }

    public static double maxwell_cdf_inv(double cdf, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_CDF_INV inverts the Maxwell CDF.
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
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double MAXWELL_CDF_INV, the corresponding argument of the CDF.
        //
    {
        int it;
        int it_max = 100;
        const double r8_huge = 1.0E+30;
        double tol = 0.0001;
        double x;
        double x1;
        double x2;
        double x3;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("MAXWELL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            case 0.0:
                x = 0.0;
                return x;
            case 1.0:
                x = r8_huge;
                return x;
        }

        x1 = 0.0;
        double cdf1 = 0.0;

        x2 = 1.0;

        for (;;)
        {
            double cdf2 = maxwell_cdf(x2, a);

            if (cdf < cdf2)
            {
                break;
            }

            x2 = 2.0 * x2;

            switch (x2)
            {
                case > 1000000.0:
                    Console.WriteLine(" ");
                    Console.WriteLine("MAXWELL_CDF_INV - Fatal error!");
                    Console.WriteLine("  Initial bracketing effort fails.");
                    return 1;
            }
        }

        //
        //  Now use bisection.
        //
        it = 0;

        for (;;)
        {
            it += 1;

            x3 = 0.5 * (x1 + x2);
            double cdf3 = maxwell_cdf(x3, a);

            if (Math.Abs(cdf3 - cdf) < tol)
            {
                x = x3;
                break;
            }

            if (it_max < it)
            {
                Console.WriteLine(" ");
                Console.WriteLine("MAXWELL_CDF_INV - Fatal error!");
                Console.WriteLine("  Iteration limit exceeded.");
                return 1;
            }

            if (cdf3 <= cdf && cdf1 <= cdf || cdf <= cdf3 && cdf <= cdf1)
            {
                x1 = x3;
                cdf1 = cdf3;
            }
            else
            {
                x2 = x3;
                double cdf2 = cdf3;
            }
        }

        return x;
    }

    public static bool maxwell_check(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_CHECK checks the parameters of the Maxwell CDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, bool MAXWELL_CHECK, is true if the parameters are legal.
        //
    {
        switch (a)
        {
            case <= 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("MAXWELL_CHECK - Warning!");
                Console.WriteLine("  A <= 0.0.");
                return false;
            default:
                return true;
        }
    }

    public static double maxwell_mean(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_MEAN returns the mean of the Maxwell PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double MAXWELL_MEAN, the mean value.
        //
    {
        double mean = Math.Sqrt(2.0) * a * Helpers.Gamma(2.0) / Helpers.Gamma(1.5);

        return mean;
    }

    public static double maxwell_pdf(double x, double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_PDF evaluates the Maxwell PDF.
        //
        //  Discussion:
        //
        //    PDF(A;X) = EXP ( - 0.5 * ( X / A )^2 ) * ( X / A )^2 /
        //      ( sqrt ( 2 ) * A * GAMMA ( 1.5 ) )
        //
        //    MAXWELL_PDF(A;X) = CHI_PDF(0,A,3;X)
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
        //    0 < X
        //
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double MAXWELL_PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case <= 0.0:
                pdf = 0.0;
                break;
            default:
            {
                double y = x / a;

                pdf = Math.Exp(-0.5 * y * y) * y * y / (Math.Sqrt(2.0) * a * Helpers.Gamma(1.5));
                break;
            }
        }

        return pdf;
    }

    public static double maxwell_sample(double a, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_SAMPLE samples the Maxwell PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double MAXWELL_SAMPLE, a sample of the PDF.
        //
    {
        double a2 = 3.0;
        double x = Chi.chi_square_sample(a2, ref seed);

        x = a * Math.Sqrt(x);

        return x;
    }

    public static double maxwell_variance(double a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAXWELL_VARIANCE returns the variance of the Maxwell PDF.
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
        //    Input, double A, the parameter of the PDF.
        //    0 < A.
        //
        //    Output, double MAXWELL_VARIANCE, the variance of the PDF.
        //
    {
        double temp = Helpers.Gamma(2.0) / Helpers.Gamma(1.5);

        double variance = a * a * (3.0 - 2.0 * temp * temp);

        return variance;
    }
}