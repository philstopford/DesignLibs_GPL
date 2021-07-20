using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8vec_variance(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_variance() returns the variance of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N, the number of entries in the vector.
            //
            //    double X[N], the vector whose variance is desired.
            //
            //  Output:
            //
            //    double R8VEC_VARIANCE, the variance of the vector entries.
            //
        {
            int i;
            double mean;
            double variance;

            mean = r8vec_mean(n, x);

            variance = 0.0;
            for (i = 0; i < n; i++)
            {
                variance = variance + (x[i] - mean) * (x[i] - mean);
            }

            if (1 < n)
            {
                variance = variance / (double) (n);
            }
            else
            {
                variance = 0.0;
            }

            return variance;
        }

        public static double r8vec_variance_circular(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_VARIANCE_CIRCULAR returns the circular variance of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 December 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double X[N], the vector whose variance is desired.
            //
            //    Output, double R8VEC_VARIANCE_CIRCULAR, the circular variance
            //    of the vector entries.
            //
        {
            int i;
            double mean;
            double sum_c;
            double sum_s;
            double value;

            mean = r8vec_mean(n, x);

            sum_c = 0.0;
            for (i = 0; i < n; i++)
            {
                sum_c = sum_c + Math.Cos(x[i] - mean);
            }

            sum_s = 0.0;
            for (i = 0; i < n; i++)
            {
                sum_s = sum_s + Math.Sin(x[i] - mean);
            }

            value = Math.Sqrt(sum_c * sum_c + sum_s * sum_s) / (double) n;

            value = 1.0 - value;

            return value;
        }

        public static double r8vec_variance_sample(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_VARIANCE_SAMPLE returns the sample variance of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N, the number of entries in the vector.
            //
            //    double X[N], the vector whose variance is desired.
            //
            //  Output:
            //
            //    double R8VEC_VARIANCE_SAMPLE, the sample variance.
            //
        {
            int i;
            double mean;
            double variance;

            mean = r8vec_mean(n, x);

            variance = 0.0;
            for (i = 0; i < n; i++)
            {
                variance = variance + (x[i] - mean) * (x[i] - mean);
            }

            if (1 < n)
            {
                variance = variance / (double) (n - 1);
            }
            else
            {
                variance = 0.0;
            }

            return variance;
        }

        public static void r8vec_variance_sample_update(int nm1, double mean_nm1,
                double variance_nm1, double xn, ref int n, ref double mean_n, ref double variance_n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_variance_sample_update updates sample variance with one new value.
            //
            //  Discussion:
            //
            //    On first call:
            //      nm1 = 0
            //      mean_nm1 = 0.0
            //      variance_nm1 = 0.0
            //      xn = first value to be handled.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John D Cook,
            //    Accurately computing running variance,
            //    https://www.johndcook.com/blog/standard_deviation/
            //
            //  Input:
            //
            //    int NM1, the number of entries in the old vector.
            //
            //    double MEAN_NM1, the mean of the old vector.
            //
            //    double VARIANCE_NM1, the variance of the old vector.
            //
            //    double XN, the new N-th entry of the vector.
            //
            //  Output:
            //
            //    int &N, the number of entries in the new vector.
            //
            //    double &MEAN_N, the mean of the new vector.
            //
            //    double &VARIANCE_N, the variance of the new vector.
            //
        {
            if (nm1 <= 0)
            {
                n = 1;
                mean_n = xn;
                variance_n = 0.0;
            }
            else
            {
                n = nm1 + 1;
                mean_n = mean_nm1 + (xn - mean_nm1) / (double) (n);
                variance_n = (variance_nm1 * (double) (nm1 - 1)
                              + (xn - mean_nm1) * (xn - mean_n)) / (double) (n - 1);
            }

            return;
        }

        public static void r8vec_variance_update(int nm1, double mean_nm1, double variance_nm1,
                double xn, ref int n, ref double mean_n, ref double variance_n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_variance_update() updates the variance with one new value.
            //
            //  Discussion:
            //
            //    On first call:
            //      nm1 = 0
            //      mean_nm1 = 0.0
            //      variance_nm1 = 0.0
            //      xn = first value to be handled.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int NM1, the number of entries in the old vector.
            //
            //    double MEAN_NM1, the mean of the old vector.
            //
            //    double VARIANCE_NM1, the variance of the old vector.
            //
            //    double XN, the new N-th entry of the vector.
            //
            //  Output:
            //
            //    int &N, the number of entries in the new vector.
            //
            //    double &MEAN_N, the mean of the new vector.
            //
            //    double &VARIANCE_N, the variance of the new vector.
            //
        {
            if (nm1 <= 0)
            {
                n = 1;
                mean_n = xn;
                variance_n = 0.0;
            }
            else
            {
                n = nm1 + 1;
                mean_n = mean_nm1 + (xn - mean_nm1) / (double) (n);
                variance_n = (variance_nm1 * (double) (nm1)
                              + (xn - mean_nm1) * (xn - mean_n)) / (double) (n);
            }
        }


    }
}