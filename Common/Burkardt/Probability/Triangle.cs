using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Triangle
{
    public static double triangle_cdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CDF evaluates the Triangle CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
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
        //    A <= B <= C and A < C.
        //
        //    Output, double TRIANGLE_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x <= a)
        {
            cdf = 0.0;
        }
        else if (x <= b)
        {
            if (Math.Abs(a - b) <= double.Epsilon)
            {
                cdf = 0.0;
            }
            else
            {
                cdf = (x - a) * (x - a) / (b - a) / (c - a);
            }
        }
        else if (x <= c)
        {
            cdf = (b - a) / (c - a)
                  + (2.0 * c - b - x) * (x - b) / (c - b) / (c - a);
        }
        else
        {
            cdf = 1.0;
        }

        return cdf;
    }

    public static double triangle_cdf_inv(double cdf, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CDF_INV inverts the Triangle CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2004
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
        //    A <= B <= C and A < C.
        //
        //    Output, double TRIANGLE_CDF_INV, the corresponding argument.
        //
    {
        double x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double d = 2.0 / (c - a);
        double cdf_mid = 0.5 * d * (b - a);

        if (cdf <= cdf_mid)
        {
            x = a + Math.Sqrt(cdf * (b - a) * (c - a));
        }
        else
        {
            x = c - Math.Sqrt((c - b) * (c - b - (cdf - cdf_mid) * (c - a)));
        }

        return x;
    }

    public static bool triangle_check(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CHECK checks the parameters of the Triangle CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A <= B <= C and A < C.
        //
        //    Output, bool TRIANGLE_CHECK, is true if the parameters are legal.
        //
    {
        if (b < a)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CHECK - Warning!");
            Console.WriteLine("  B < A.");
            return false;
        }

        if (c < b)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CHECK - Warning!");
            Console.WriteLine("  C < B.");
            return false;
        }

        if (Math.Abs(a - c) <= double.Epsilon)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CHECK - Warning!");
            Console.WriteLine("  A == C.");
            return false;
        }

        return true;
    }

    public static double triangle_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_MEAN returns the mean of the Triangle PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A <= B <= C and A < C.
        //
        //    Output, double TRIANGLE_MEAN, the mean of the discrete uniform PDF.
        //
    {
        double mean = a + (c + b - 2.0 * a) / 3.0;

        return mean;
    }

    public static double triangle_pdf(double x, double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_PDF evaluates the Triangle PDF.
        //
        //  Discussion:
        //
        //    Given points A <= B <= C, the probability is 0 to the left of A,
        //    rises linearly to a maximum of 2/(C-A) at B, drops linearly to zero
        //    at C, and is zero for all values greater than C.
        //
        //    PDF(A,B,C;X)
        //      = 2 * ( X - A ) / ( B - A ) / ( C - A ) for A <= X <= B
        //      = 2 * ( C - X ) / ( C - B ) / ( C - A ) for B <= X <= C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
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
        //    A <= B <= C and A < C.
        //
        //    Output, double TRIANGLE_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x <= a)
        {
            pdf = 0.0;
        }
        else if (x <= b)
        {
            if (Math.Abs(a - b) <= double.Epsilon)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = 2.0 * (x - a) / (b - a) / (c - a);
            }
        }
        else if (x <= c)
        {
            if (Math.Abs(b - c) <= double.Epsilon)
            {
                pdf = 0.0;
            }
            else
            {
                pdf = 2.0 * (c - x) / (c - b) / (c - a);
            }
        }
        else
        {
            pdf = 0.0;
        }

        return pdf;
    }

    public static double triangle_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_SAMPLE samples the Triangle PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A <= B <= C and A < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double TRIANGLE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = triangle_cdf_inv(cdf, a, b, c);

        return x;
    }

    public static double triangle_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_VARIANCE returns the variance of the Triangle PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the PDF.
        //    A <= B <= C and A < C.
        //
        //    Output, double VARIANCE, the variance of the PDF.
        //
    {
        double variance = ((c - a) * (c - a)
                           - (c - a) * (b - a)
                           + (b - a) * (b - a)) / 18.0;

        return variance;
    }
}