using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class EmpiricalDiscrete
    {
        static double empirical_discrete_cdf(double x, int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_CDF evaluates the Empirical Discrete CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B[A], the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C[A], the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, double EMPIRICAL_DISCRETE_CDF, the value of the CDF.
        //
        {
            double cdf = 0.0;

            double bsum = typeMethods.r8vec_sum(a, b);

            for (int i = 1; i <= a; i++)
            {
                if (x < c[i - 1])
                {
                    return cdf;
                }

                cdf = cdf + b[i - 1] / bsum;
            }

            return cdf;
        }

        static double empirical_discrete_cdf_inv(double cdf, int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_CDF_INV inverts the Empirical Discrete CDF.
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
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B(A), the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C(A), the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, double EMPIRICAL_DISCRETE_CDF_INV, the smallest argument
        //    whose CDF is greater than or equal to CDF.
        //
        {
            if (cdf < 0.0 || 1.0 < cdf)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EMPIRICAL_DISCRETE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return (1);
            }

            double bsum = typeMethods.r8vec_sum(a, b);

            double x = c[0];
            double cdf2 = b[0] / bsum;

            for (int i = 1; i < a; i++)
            {
                if (cdf <= cdf2)
                {
                    return x;
                }

                x = c[i];
                cdf2 = cdf2 + b[i] / bsum;
            }

            return x;
        }

        static bool empirical_discrete_check(int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_CHECK checks the parameters of the Empirical Discrete CDF.
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
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B[A], the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C[A], the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, bool EMPIRICAL_DISCRETE_CHECK, is true if the parameters
        //    are legal.
        //
        {

            if (a <= 0)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EMPIRICAL_DISCRETE_CHECK - Warning!");
                Console.WriteLine("  A must be positive.");
                Console.WriteLine("  Input A = " + a + "");
                Console.WriteLine("  A is the number of weights.");
                return false;
            }

            for (int i = 0; i < a; i++)
            {
                if (b[i] < 0.0)
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("EMPIRICAL_DISCRETE_CHECK - Warning!");
                    Console.WriteLine("  B[" + i + "] < 0.");
                    Console.WriteLine("  But all B values must be nonnegative.");
                    return false;
                }
            }

            bool positive = false;

            for (int i = 0; i < a; i++)
            {
                if (0.0 < b[i])
                {
                    positive = true;
                }
            }

            if (!positive)
            {
                Console.WriteLine(" ");
                Console.WriteLine("EMPIRICAL_DISCRETE_CHECK - Warning!");
                Console.WriteLine("  All B(*) = 0.");
                Console.WriteLine("  But at least one B values must be nonzero.");
                return false;
            }

            for (int i = 0; i < a; i++)
            {
                for (int j = i + 1; j < a; j++)
                {
                    if (c[i] == c[j])
                    {
                        Console.WriteLine(" ");
                        Console.WriteLine("EMPIRICAL_DISCRETE_CHECK - Warning!");
                        Console.WriteLine("  All values C must be unique.");
                        Console.WriteLine("  But at least two values are identical.");
                        return false;
                    }
                }
            }

            for (int i = 0; i < a - 1; i++)
            {
                if (c[i + 1] < c[i])
                {
                    Console.WriteLine(" ");
                    Console.WriteLine("EMPIRICAL_DISCRETE_CHECK - Warning!");
                    Console.WriteLine("  The values in C must be in ascending order.");
                    return false;
                }
            }

            return true;
        }

        static double empirical_discrete_mean(int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_MEAN returns the mean of the Empirical Discrete PDF.
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
        //  Parameters:
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B(A), the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C(A), the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, double EMPIRICAL_DISCRETE_MEAN, the mean of the PDF.
        //
        {
            int i;
            double mean;

            mean = 0.0;
            for (i = 0; i < a; i++)
            {
                mean = mean + b[i] * c[i];
            }

            mean = mean / typeMethods.r8vec_sum(a, b);

            return mean;
        }

        static double empirical_discrete_pdf(double x, int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_PDF evaluates the Empirical Discrete PDF.
        //
        //  Discussion:
        //
        //    A set of A values C(1:A) are assigned nonnegative weights B(1:A),
        //    with at least one B nonzero.  The probability of C(I) is the
        //    value of B(I) divided by the sum of the weights.
        //
        //    The C's must be distinct, and given in ascending order.
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
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B(A), the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C(A), the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, double EMPIRICAL_DISCRETE_PDF, the value of the PDF.
        //
        {
            double pdf;

            for (int i = 0; i <= a; i++)
            {
                if (x == c[i])
                {
                    pdf = b[i] / typeMethods.r8vec_sum(a, b);
                    return pdf;
                }
            }

            pdf = 0.0;

            return pdf;
        }

        static double empirical_discrete_sample(int a, double[] b, double[] c, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_SAMPLE samples the Empirical Discrete PDF.
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
        //  Parameters:
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B(A), the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C(A), the values.
        //    The values must be distinct and in ascending order.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double EMPIRICAL_DISCRETE_SAMPLE, a sample of the PDF.
        //
        {
            double cdf = UniformRNG.r8_uniform_01(ref seed);

            double x = empirical_discrete_cdf_inv(cdf, a, b, c);

            return x;
        }

        static double empirical_discrete_variance(int a, double[] b, double[] c )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMPIRICAL_DISCRETE_VARIANCE returns the variance of the Empirical Discrete PDF.
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
        //  Parameters:
        //
        //    Input, int A, the number of values.
        //    0 < A.
        //
        //    Input, double B(A), the weights of each value.
        //    0 <= B(1:A) and at least one value is nonzero.
        //
        //    Input, double C(A), the values.
        //    The values must be distinct and in ascending order.
        //
        //    Output, double EMPIRICAL_DISCRETE_VARIANCE, the variance of the PDF.
        //
        {
            double bsum = typeMethods.r8vec_sum(a, b);

            double mean = empirical_discrete_mean(a, b, c);

            double variance = 0.0;

            for (int i = 0; i < a; i++)
            {
                variance = variance + (b[i] / bsum) * Math.Pow(c[i] - mean, 2);
            }

            return variance;
        }
    }
}