﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Discrete
{
    public static double discrete_cdf(int x, int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_CDF evaluates the Discrete CDF.
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
        //    Input, int X, the item whose probability is desired.
        //
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of outcomes
        //    1 through A.  Each entry must be nonnegative.
        //
        //    Output, double DISCRETE_CDF, the value of the CDF.
        //
    {
        double cdf = 0;

        switch (x)
        {
            case < 1:
                cdf = 0.0;
                break;
            default:
            {
                if (x < a)
                {
                    cdf = typeMethods.r8vec_sum(x, b) / typeMethods.r8vec_sum(a, b);
                }
                else if (a <= x)
                {
                    cdf = 1.0;
                }

                break;
            }
        }

        return cdf;
    }

    public static int discrete_cdf_inv(double cdf, int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_CDF_INV inverts the Discrete CDF.
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
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of outcomes
        //    1 through A.  Each entry must be nonnegative.
        //
        //    Output, int DISCRETE_CDF_INV, the corresponding argument for which
        //    CDF(X-1) < CDF <= CDF(X)
        //
    {
        int x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("DISCRETE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double b_sum = typeMethods.r8vec_sum(a, b);

        double cum = 0.0;

        for (int j = 1; j <= a; j++)
        {
            cum += b[j - 1] / b_sum;

            if (!(cdf <= cum))
            {
                continue;
            }

            x = j;
            return x;
        }

        x = a;

        return x;
    }

    public static bool discrete_check(int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_CHECK checks the parameters of the Discrete CDF.
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
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of
        //    outcomes 1 through A.  Each entry must be nonnegative.
        //
        //    Output, bool DISCRETE_CHECK, is true if the parameters are legal.
        //
    {
        for (int j = 0; j < a; j++)
        {
            switch (b[j])
            {
                case < 0.0:
                    Console.WriteLine(" ");
                    Console.WriteLine("DISCRETE_CHECK - Warning!");
                    Console.WriteLine("  Negative probabilities not allowed.");
                    return false;
            }
        }

        double b_sum = typeMethods.r8vec_sum(a, b);

        switch (b_sum)
        {
            case 0.0:
                Console.WriteLine(" ");
                Console.WriteLine("DISCRETE_CHECK - Warning!");
                Console.WriteLine("  Total probablity is zero.");
                return false;
            default:
                return true;
        }
    }

    public static double discrete_mean(int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_MEAN evaluates the mean of the Discrete PDF.
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
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of
        //    outcomes 1 through A.  Each entry must be nonnegative.
        //
        //    Output, double DISCRETE_MEAN, the mean of the PDF.
        //
    {
        double b_sum = typeMethods.r8vec_sum(a, b);

        double mean = 0.0;
        for (int j = 0; j < a; j++)
        {
            mean += (j + 1) * b[j];
        }

        mean /= b_sum;

        return mean;
    }

    public static double discrete_pdf(int x, int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_PDF evaluates the Discrete PDF.
        //
        //  Discusion:
        //
        //    PDF(A,B;X) = B(X) if 1 <= X <= A
        //                = 0    otherwise
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
        //    Input, int X, the item whose probability is desired.
        //
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of
        //    outcomes 1 through A.  Each entry must be nonnegative.
        //
        //    Output, double DISCRETE_PDF, the value of the PDF.
        //
    {
        double b_sum = typeMethods.r8vec_sum(a, b);

        double pdf = x switch
        {
            >= 1 when x <= a => b[x - 1] / b_sum,
            _ => 0.0
        };

        return pdf;
    }

    public static int discrete_sample(int a, double[] b, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_SAMPLE samples the Discrete PDF.
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
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of
        //    outcomes 1 through A.  Each entry must be nonnegative.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int DISCRETE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = discrete_cdf_inv(cdf, a, b);

        return x;
    }

    public static double discrete_variance(int a, double[] b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISCRETE_VARIANCE evaluates the variance of the Discrete PDF.
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
        //    Input, int A, the number of probabilities assigned.
        //
        //    Input, double B[A], the relative probabilities of
        //    outcomes 1 through A.  Each entry must be nonnegative.
        //
        //    Output, double DISCRETE_VARIANCE, the variance of the PDF.
        //
    {
        double b_sum = typeMethods.r8vec_sum(a, b);

        double mean = 0.0;
        for (int j = 1; j <= a; j++)
        {
            mean += j * b[j - 1];
        }

        mean /= b_sum;

        double variance = 0.0;
        for (int j = 1; j <= a; j++)
        {
            variance += b[j - 1] * Math.Pow(j - mean, 2);
        }

        variance /= b_sum;

        return variance;
    }
}