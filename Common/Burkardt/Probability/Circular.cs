using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Circular
    {
        static double[] circular_normal_01_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_01_MEAN returns the mean of the Circular Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CIRCULAR_01_MEAN[2], the mean of the PDF.
        //
        {
                double[] mean = new double[2];

            mean[0] = 0.0;
            mean[1] = 0.0;

            return mean;
        }

        static double circular_normal_01_pdf(double[] x )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_01_PDF evaluates the Circular Normal 01 PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = EXP ( - 0.5 * ( X(1)^2 + X(2)^2 ) ) / ( 2 * PI )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[2], the argument of the PDF.
        //
        //    Output, double CIRCULAR_NORMAL_01_PDF, the value of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double pdf = Math.Exp(-0.5 * (x[0] * x[0] + x[1] * x[1])) / (2.0 * r8_pi);

            return pdf;
        }

        static double[] circular_normal_01_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_01_SAMPLE samples the Circular Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CIRCULAR_NORMAL_01_SAMPLE[2], a sample of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double[] x = new double[2];

            double v1 = UniformRNG.r8_uniform_01(ref seed);
            double v2 = UniformRNG.r8_uniform_01(ref seed);

            x[0] = Math.Sqrt(-2.0 *Math.Log(v1)) * Math.Cos(2.0 * r8_pi * v2);
            x[1] = Math.Sqrt(-2.0 * Math.Log(v1)) * Math.Sin(2.0 * r8_pi * v2);

            return x;
        }

        static double[] circular_normal_01_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_01_VARIANCE returns the variance of the Circular Normal 01 PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, double CIRCULAR_NORMAL_01_VARIANCE[2], the variance of the PDF.
        //
        {
            double[] variance = new double[2];

            variance[0] = 1.0;
            variance[1] = 1.0;

            return variance;
        }

        static double[] circular_normal_mean(double[] a, double b )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_MEAN returns the mean of the Circular Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], a parameter of the PDF, the mean value.
        //
        //    Input, double B, a parameter of the PDF, the standard deviation.
        //
        //    Output, double CIRCULAR_MEAN[2], the mean of the PDF.
        //
        {
            double[] mean = new double[2];

            mean[0] = a[0];
            mean[1] = a[1];

            return mean;
        }

        static double circular_normal_pdf(double[] x, double[] a, double b )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_PDF evaluates the Circular Normal PDF.
        //
        //  Discussion:
        //
        //    PDF(X) = EXP ( - 0.5D+00 * ( ( (X(1)-A(1))^2 + (X(2)-A(2))^2 ) / B^2 ) 
        //      / ( 2 * PI * B^2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X[2], the argument of the PDF.
        //
        //    Input, double A[2], a parameter of the PDF, the mean value.
        //
        //    Input, double B, a parameter of the PDF, the standard deviation.
        //
        //    Output, double CIRCULAR_NORMAL_PDF, the value of the PDF.
        //
        {
            double d;
            double pdf;
            const double r8_pi = 3.14159265358979323;

            d = (Math.Pow(x[0] - a[0], 2)
                 + Math.Pow(x[1] - a[1], 2)) / Math.Pow(b, 2);

            pdf = Math.Exp(-0.5 * d) / (2.0 * b * b * r8_pi);

            return pdf;
        }

        static double[] circular_normal_sample(double[] a, double b, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_SAMPLE samples the Circular Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], a parameter of the PDF, the mean value.
        //
        //    Input, double B, a parameter of the PDF, the standard deviation.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double CIRCULAR_NORMAL_SAMPLE[2], a sample of the PDF.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double[] x = new double[2];

            double v1 = UniformRNG.r8_uniform_01(ref seed);
            double v2 = UniformRNG.r8_uniform_01(ref seed);

            double r = Math.Sqrt(-2.0 * Math.Log(v1));

            x[0] = a[0] + b * r * Math.Cos(2.0 * r8_pi * v2);
            x[1] = a[1] + b * r * Math.Sin(2.0 * r8_pi * v2);

            return x;
        }

        static double[] circular_normal_variance(double[] a, double b )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCULAR_NORMAL_VARIANCE returns the variance of the Circular Normal PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], a parameter of the PDF, the mean value.
        //
        //    Input, double B, a parameter of the PDF, the standard deviation.
        //
        //    Output, double CIRCULAR_NORMAL_VARIANCE[2], the variance of the PDF.
        //
        {
            double[] variance = new double[2];

            variance[0] = b;
            variance[1] = b;

            return variance;
        }

    }
}