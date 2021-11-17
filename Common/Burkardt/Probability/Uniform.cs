using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class Uniform
{
    public static double uniform_01_cdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_CDF evaluates the Uniform 01 CDF.
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
        //    Input, double X, the argument of the CDF.
        //
        //    Output, double UNIFORM_01_CDF, the value of the CDF.
        //
    {
        double cdf = x switch
        {
            < 0.0 => 0.0,
            > 1.0 => 1.0,
            _ => x
        };

        return cdf;
    }

    public static double uniform_01_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_CDF_INV inverts the Uniform 01 CDF.
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
        //    Output, double UNIFORM_01_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("UNIFORM_01_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                double x = cdf;

                return x;
            }
        }
    }

    public static double uniform_01_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_MEAN returns the mean of the Uniform 01 PDF.
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
        //    Output, double UNIFORM_01_MEAN, the mean of the discrete uniform PDF.
        //
    {
        double mean = 0.5;

        return mean;
    }

    public static double uniform_01_pdf(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_PDF evaluates the Uniform 01 PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = 1 for 0 <= X <= 1
        //           = 0 otherwise
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the PDF.
        //    0.0 <= X <= 1.0.
        //
        //    Output, double PDF, the value of the PDF.
        //
    {
        double pdf;

        switch (x)
        {
            case < 0.0:
            case > 1.0:
                pdf = 0.0;
                break;
            default:
                pdf = 1.0;
                break;
        }

        return pdf;
    }

    public static double uniform_01_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_SAMPLE is a random number generator.
        //
        //  Discussion:
        //
        //    SEED = SEED * (7^5) mod (2^31 - 1)
        //    UNIFORM_01_SAMPLE = SEED * / ( 2^31 - 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, the integer "seed" used to generate
        //    the output random number, and updated in preparation for the
        //    next one.  SEED should not be zero.
        //
        //    Output, double UNIFORM_01_SAMPLE, a random value between 0 and 1.
        //
        //  Local Parameters:
        //
        //    Local, int IA = 7^5.
        //
        //    Local, int IB = 2^15.
        //
        //    Local, int IB16 = 2^16.
        //
        //    Local, int IP = 2^31 - 1.
        //
    {
        int ia = 16807;
        int ib15 = 32768;
        int ib16 = 65536;
        int ip = 2147483647;
        int iprhi;
        int ixhi;
        int k;
        int leftlo;
        int loxa;
        double temp;
        seed = seed switch
        {
            //
            //  Don't let SEED be 0.
            //
            0 => ip,
            _ => seed
        };

        //
        //  Get the 15 high order bits of SEED.
        //
        ixhi = seed / ib16;
        //
        //  Get the 16 low bits of SEED and form the low product.
        //
        loxa = (seed - ixhi * ib16) * ia;
        //
        //  Get the 15 high order bits of the low product.
        //
        leftlo = loxa / ib16;
        //
        //  Form the 31 highest bits of the full product.
        //
        iprhi = ixhi * ia + leftlo;
        //
        //  Get overflow past the 31st bit of full product.
        //
        k = iprhi / ib15;
        //
        //  Assemble all the parts and presubtract IP.  The parentheses are essential.
        //
        seed = loxa - leftlo * ib16 - ip +
               (iprhi - k * ib15) * ib16 + k;
        switch (seed)
        {
            //
            //  Add IP back in if necessary.
            //
            case < 0:
                seed += ip;
                break;
        }

        //
        //  Multiply by 1 / (2^31-1).
        //
        temp = seed * 4.656612875E-10;

        return temp;
    }

    public static double uniform_01_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_VARIANCE returns the variance of the Uniform 01 PDF.
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
        //    Output, double UNIFORM_01_VARIANCE, the variance of the PDF.
        //
    {
        double variance = 1.0 / 12.0;

        return variance;
    }

    public static double[] uniform_01_order_sample(int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_01_ORDER_SAMPLE samples the Uniform 01 Order PDF.
        //
        //  Discussion:
        //
        //    In effect, this routine simply generates N samples of the
        //    Uniform 01 PDF; but it generates them in order.  (Actually,
        //    it generates them in descending order, but stores them in
        //    the array in ascending order).  This saves the work of
        //    sorting the results.  Moreover, if the order statistics
        //    for another PDF are desired, and the inverse CDF is available,
        //    then the desired values may be generated, presorted, by
        //    calling this routine and using the results as input to the
        //    inverse CDF routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 168.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the sample.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_01_ORDER_SAMPLE[N], N samples of the Uniform 01 PDF, in
        //    ascending order.
        //
    {
        double[] x = new double[n];

        double v = 1.0;

        for (int i = n - 1; 0 <= i; i--)
        {
            double u = UniformRNG.r8_uniform_01(ref seed);
            v *= Math.Pow(u, 1.0 / (i + 1));
            x[i] = v;
        }

        return x;
    }

    public static double uniform_cdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_CDF evaluates the Uniform CDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, double UNIFORM_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x < a)
        {
            cdf = 0.0;
        }
        else if (b < x)
        {
            cdf = 1.0;
        }
        else
        {
            cdf = (x - a) / (b - a);
        }

        return cdf;
    }

    public static double uniform_cdf_inv(double cdf, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_CDF_INV inverts the Uniform CDF.
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
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, double UNIFORM_CDF_INV, the corresponding argument.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("UNIFORM_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
            default:
            {
                double x = a + (b - a) * cdf;

                return x;
            }
        }
    }

    public static bool uniform_check(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_CHECK checks the parameters of the Uniform CDF.
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
        //    A < B.
        //
        //    Output, bool UNIFORM_CHECK, is true if the parameters are legal.
        //
    {
        if (b <= a)
        {
            Console.WriteLine(" ");
            Console.WriteLine("UNIFORM_CHECK - Warning!");
            Console.WriteLine("  B <= A.");
            return false;
        }

        return true;
    }

    public static double uniform_mean(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_MEAN returns the mean of the Uniform PDF.
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
        //    A < B.
        //
        //    Output, double UNIFORM_MEAN, the mean of the discrete uniform PDF.
        //
    {
        double mean = 0.5 * (a + b);

        return mean;
    }

    public static double uniform_pdf(double x, double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_PDF evaluates the Uniform PDF.
        //
        //  Discussion:
        //
        //    The Uniform PDF is also known as the "Rectangular" or "de Moivre" PDF.
        //
        //    PDF(A,B;X) = 1 / ( B - A ) for A <= X <= B
        //               = 0 otherwise
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
        //    Input, double A, B, the parameters of the PDF.
        //    A < B.
        //
        //    Output, double UNIFORM_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x < a || b < x)
        {
            pdf = 0.0;
        }
        else
        {
            pdf = 1.0 / (b - a);
        }

        return pdf;
    }

    public static double uniform_sample(double a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_SAMPLE samples the Uniform PDF.
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
        //    A < B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        double x = uniform_cdf_inv(cdf, a, b);

        return x;
    }

    public static double uniform_variance(double a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_VARIANCE returns the variance of the Uniform PDF.
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
        //    A < B.
        //
        //    Output, double UNIFORM_VARIANCE, the variance of the PDF.
        //
    {
        double variance = (b - a) * (b - a) / 12.0;

        return variance;
    }

    public static double uniform_discrete_cdf(int x, int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_CDF evaluates the Uniform Discrete CDF.
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
        //    Input, int X, the argument of the CDF.
        //
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, double UNIFORM_DISCRETE_CDF, the value of the CDF.
        //
    {
        double cdf;

        if (x < a)
        {
            cdf = 0.0;
        }
        else if (b < x)
        {
            cdf = 1.0;
        }
        else
        {
            cdf = (x + 1 - a) / (double) (b + 1 - a);
        }

        return cdf;
    }

    public static int uniform_discrete_cdf_inv(double cdf, int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_CDF_INV inverts the Uniform Discrete CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 November 2004
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
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, int UNIFORM_DISCRETE_CDF_INV, the smallest argument whose
        //    CDF is greater than or equal to CDF.
        //
    {
        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("UNIFORM_DISCRETE_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double a2 = a - 0.5;
        double b2 = b + 0.5;
        double x2 = a + cdf * (b2 - a2);

        int x = (int)typeMethods.r8_nint(x2);

        x = Math.Max(x, a);
        x = Math.Min(x, b);

        return x;
    }

    public static bool uniform_discrete_check(int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_CHECK checks the parameters of the Uniform discrete CDF.
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
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, bool UNIFORM_DISCRETE_CHECK, is true if the parameters
        //    are legal.
        //
    {
        if (b < a)
        {
            Console.WriteLine(" ");
            Console.WriteLine("UNIFORM_DISCRETE_CHECK - Warning!");
            Console.WriteLine("  B < A.");
            return false;
        }

        return true;
    }

    public static double uniform_discrete_mean(int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_MEAN returns the mean of the Uniform discrete PDF.
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
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, double UNIFORM_DISCRETE_MEAN, the mean of the PDF.
        //
    {
        double mean = 0.5 * (a + b);

        return mean;
    }

    public static double uniform_discrete_pdf(int x, int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_PDF evaluates the Uniform discrete PDF.
        //
        //  Discussion:
        //
        //    The Uniform Discrete PDF is also known as the "Rectangular"
        //    Discrete PDF.
        //
        //    PDF(A,B;X) = 1 / ( B + 1 - A ) for A <= X <= B.
        //
        //    The parameters define the interval of integers
        //    for which the PDF is nonzero.
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
        //    Input, int X, the argument of the PDF.
        //
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, double UNIFORM_DISCRETE_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x < a || b < x)
        {
            pdf = 0.0;
        }
        else
        {
            pdf = 1.0 / (b + 1 - a);
        }

        return pdf;
    }

    public static int uniform_discrete_sample(int a, int b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_SAMPLE samples the Uniform discrete PDF.
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
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int UNIFORM_DISCRETE_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = uniform_discrete_cdf_inv(cdf, a, b);

        return x;
    }

    public static double uniform_discrete_variance(int a, int b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_DISCRETE_VARIANCE returns the variance of the Uniform discrete PDF.
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
        //    Input, int A, B, the parameters of the PDF.
        //    A <= B.
        //
        //    Output, double UNIFORM_DISCRETE_VARIANCE, the variance of the PDF.
        //
    {
        double variance = ((b + 1 - a) * (b + 1 - a) - 1) / 12.0;

        return variance;
    }

    public static double[] uniform_nsphere_sample(int n, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_NSPHERE_SAMPLE samples the Uniform Unit Sphere PDF.
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
        //    Jerry Banks, editor,
        //    Handbook of Simulation,
        //    Engineering and Management Press Books, 1998, page 168.
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the sphere.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_NSPHERE_SAMPLE[N], a point on the unit N
        //    sphere, chosen with a uniform probability.
        //
    {
        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            x[i] = Normal.normal_01_sample(ref seed);
        }

        double sum = 0.0;
        for (int i = 0; i < n; i++)
        {
            sum += x[i] * x[i];
        }

        sum = Math.Sqrt(sum);

        for (int i = 0; i < n; i++)
        {
            x[i] /= sum;
        }

        return x;
    }

}