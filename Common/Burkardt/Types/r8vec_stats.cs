using System;
using System.Linq;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8vec_correlation(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_correlation returns the correlation of two R8VEC's.
            //
            //  Discussion:
            //
            //    The correlation coefficient is also known as Pearson's r coefficient.
            //
            //    It must be the case that -1 <= r <= +1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 August 2019
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int n: the dimension of the vectors.
            //
            //    double x[n], y[n]: the vectors.
            //
            //  Output:
            //
            //    double r8vec_correlation, the correlation of X and Y.
            //
        {
            double dot;
            int i;
            double r;
            double x_mean;
            double x_std;
            double y_mean;
            double y_std;

            if (n <= 1)
            {
                r = 0.0;
            }
            else
            {
                x_mean = r8vec_mean(n, x);
                x_std = r8vec_std_sample(n, x);

                y_mean = r8vec_mean(n, y);
                y_std = r8vec_std_sample(n, y);

                dot = 0.0;
                for (i = 0; i < n; i++)
                {
                    dot = dot + (x[i] - x_mean) * (y[i] - y_mean);
                }

                r = dot / x_std / y_std / (double) (n - 1);
            }

            return r;
        }

        public static double r8vec_covariance(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_COVARIANCE computes the covariance of two vectors.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 April 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the two vectors.
            //
            //    Input, double X[N], Y[N], the two vectors.
            //
            //    Output, double R8VEC_COVARIANCE, the covariance of the two vectors.
            //
        {
            int i;
            double value;
            double x_average;
            double y_average;

            x_average = 0.0;
            for (i = 0; i < n; i++)
            {
                x_average = x_average + x[i];
            }

            x_average = x_average / (double) (n);

            y_average = 0.0;
            for (i = 0; i < n; i++)
            {
                y_average = y_average + x[i];
            }

            y_average = y_average / (double) (n);

            value = 0.0;
            for (i = 0; i < n; i++)
            {
                value = value + (x[i] - x_average) * (y[i] - y_average);
            }

            value = value / (double) (n - 1);

            return value;
        }


        public static double r8vec_mean(int n, double[] x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MEAN returns the mean of an R8VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 May 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double X[N], the vector whose mean is desired.
            //
            //    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
            //
        {
            if (x.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, x.Length);

            if (n == x.Length)
            {
                return x.Average();
            }

            return x.Take(n).Average();
        }
        
        public static double r8vec_mean_geometric ( int n, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MEAN_GEOMETRIC returns the geometric mean of an R8VEC.
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
            //    27 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double X[N], the vector whose geometric mean is desired.
            //    All entries should be nonnegative.
            //
            //    Output, double R8VEC_MEAN_GEOMETRIC, the geometric mean of the 
            //    vector entries.
            //
        {
            int i;
            double mean;

            mean = 0.0;
            for ( i = 0; i < n; i++ )
            {
                mean = mean + Math.Log ( x[i] );
            }

            mean = mean / ( double ) n;
            mean = Math.Exp ( mean );

            return mean;
        }

        public static double[] r8vec_mean_running ( int n, double[] v )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MEAN_RUNNING computes the running means of an R8VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 February 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items.
            //
            //    Input, double V(N), the data.
            //
            //    Output, double R8VEC_MEAN_RUNNING[N+1], the running averages.  A[i] is 
            //    the average value of the first I-1 values in V.
            //
        {
            double[] a;
            int i;

            a = new double[n+1];
            //
            //  Sum.
            //
            a[0] = 0.0;
            for ( i = 1; i < n + 1; i++ )
            {
                a[i] = a[i-1] + v[i-1];
            }
            //
            //  Average.
            //
            for ( i = 1; i < n + 1; i++ )
            {
                a[i] = a[i] / ( double ) ( i );
            }

            return a;
        }
        
        public static double r8vec_median ( int n, ref double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MEDIAN returns the median of an unsorted R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    Hoare's algorithm is used.  The values of the vector are
            //    rearranged by this routine.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input/output, double A[N], the array to search.  On output,
            //    the order of the elements of A has been somewhat changed.
            //
            //    Output, double R8VEC_MEDIAN, the value of the median of A.
            //
        {
            int k;
            double median;

            k = ( n + 1 ) / 2;

            median = r8vec_frac ( n, ref a, k );

            return median;
        }

        public static double r8vec_circular_variance(int n, double[] x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC
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
            //    Input, double X(N), the vector whose variance is desired.
            //
            //    Output, double R8VEC_CIRCULAR VARIANCE, the circular variance
            //    of the vector entries.
            //
        {
            double mean = r8vec_mean(n, x);

            double sum_c = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum_c = sum_c + Math.Cos(x[i] - mean);
            }

            double sum_s = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum_s = sum_s + Math.Sin(x[i] - mean);
            }

            double value = Math.Sqrt(sum_c * sum_c + sum_s * sum_s) / (double) n;

            value = 1.0 - value;

            return value;
        }

        public static double r8vec_std(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_std() returns the standard deviation of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The standard deviation of a vector X of length N is defined as
            //
            //      mean ( X(1:n) ) = sum ( X(1:n) ) / n
            //
            //      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N, the number of entries in the vector.
            //    N should be at least 2.
            //
            //    double A[N], the vector.
            //
            //  Output:
            //
            //    double R8VEC_STD, the standard deviation of the vector.
            //
        {
            int i;
            double mean;
            double std;

            if (n < 2)
            {
                std = 0.0;
            }
            else
            {
                mean = 0.0;
                for (i = 0; i < n; i++)
                {
                    mean = mean + a[i];
                }

                mean = mean / ((double) n);

                std = 0.0;
                for (i = 0; i < n; i++)
                {
                    std = std + (a[i] - mean) * (a[i] - mean);
                }

                std = Math.Sqrt(std / ((double) (n)));
            }

            return std;
        }

        public static double r8vec_std_sample(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_std_sample() returns the sample standard deviation of an R8VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    The standard deviation of a vector X of length N is defined as
            //
            //      mean ( X(1:n) ) = sum ( X(1:n) ) / n
            //
            //      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n - 1 ) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N, the number of entries in the vector.
            //    N should be at least 2.
            //
            //    double A[N], the vector.
            //
            //  Output:
            //
            //    double R8VEC_STD_SAMPLE, the sample standard deviation of the vector.
            //
        {
            int i;
            double mean;
            double std;

            if (n < 2)
            {
                std = 0.0;
            }
            else
            {
                mean = 0.0;
                for (i = 0; i < n; i++)
                {
                    mean = mean + a[i];
                }

                mean = mean / ((double) n);

                std = 0.0;
                for (i = 0; i < n; i++)
                {
                    std = std + (a[i] - mean) * (a[i] - mean);
                }

                std = Math.Sqrt(std / ((double) (n - 1)));
            }

            return std;
        }

        public static void r8vec_std_sample_update(int nm1, double mean_nm1, double std_nm1,
                double xn, ref int n, ref double mean_n, ref double std_n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_std_sample_update() updates sample standard deviation with one new value.
            //
            //  Discussion:
            //
            //    On first call:
            //      nm1 = 0
            //      mean_nm1 = 0.0
            //      std_nm1 = 0.0
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
            //    double STD_NM1, the sample standard deviation of the old vector.
            //
            //    double XN, the new N-th entry of the vector.
            //
            //  Output:
            //
            //    int &N, the number of entries in the new vector.
            //
            //    double &MEAN_N, the mean of the new vector.
            //
            //    double &STD_N, the sample standard deviation of the new vector.
            //
        {
            if (nm1 <= 0)
            {
                n = 1;
                mean_n = xn;
                std_n = 0.0;
            }
            else
            {
                n = nm1 + 1;
                mean_n = mean_nm1 + (xn - mean_nm1) / (double) (n);
                std_n = Math.Sqrt((std_nm1 * std_nm1 * (double) (nm1 - 1)
                              + (xn - mean_nm1) * (xn - mean_n)) / (double) (n - 1));
            }

            return;
        }

        public static void r8vec_std_update(int nm1, double mean_nm1, double std_nm1,
                double xn, ref int n, ref double mean_n, ref double std_n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8vec_std_update() updates standard deviation with one new value.
            //
            //  Discussion:
            //
            //    On first call:
            //      nm1 = 0
            //      mean_nm1 = 0.0
            //      std_nm1 = 0.0
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
            //    double STD_NM1, the standard deviation of the old vector.
            //
            //    double XN, the new N-th entry of the vector.
            //
            //  Output:
            //
            //    int &N, the number of entries in the new vector.
            //
            //    double &MEAN_N, the mean of the new vector.
            //
            //    double &STD_N, the standard deviation of the new vector.
            //
        {
            if (nm1 <= 0)
            {
                n = 1;
                mean_n = xn;
                std_n = 0.0;
            }
            else
            {
                n = nm1 + 1;
                mean_n = mean_nm1 + (xn - mean_nm1) / (double) (n);
                std_n = Math.Sqrt((std_nm1 * std_nm1 * (double) (nm1)
                              + (xn - mean_nm1) * (xn - mean_n)) / (double) (n));
            }

            return;
        }


    }
}