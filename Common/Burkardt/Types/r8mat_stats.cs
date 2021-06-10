using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8mat_covariance(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_COVARIANCE computes the sample covariance of a set of vector data.
            //
            //  Discussion:
            //
            //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 June 2013
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int M, the size of a single data vectors.
            //
            //    Input, int N, the number of data vectors.
            //    N should be greater than 1.
            //
            //    Input, double X[M*N], an array of N data vectors, each
            //    of length M.
            //
            //    Output, double C[M*M], the covariance matrix for the data.
            //
        {
            double[] c;
            int i;
            int j;
            int k;
            double[] x_mean;

            c = new double[m * m];
            for (j = 0; j < m; j++)
            {
                for (i = 0; i < m; i++)
                {
                    c[i + j * m] = 0.0;
                }
            }

            //
            //  Special case of N = 1.
            //
            if (n == 1)
            {
                for (i = 0; i < m; i++)
                {
                    c[i + i * m] = 1.0;
                }

                return c;
            }

            //
            //  Determine the sample means.
            //
            x_mean = new double[m];
            for (i = 0; i < m; i++)
            {
                x_mean[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    x_mean[i] = x_mean[i] + x[i + j * m];
                }

                x_mean[i] = x_mean[i] / (double) (n);
            }

            //
            //  Determine the sample covariance.
            //
            for (j = 0; j < m; j++)
            {
                for (i = 0; i < m; i++)
                {
                    for (k = 0; k < n; k++)
                    {
                        c[i + j * m] = c[i + j * m]
                                       + (x[i + k * m] - x_mean[i]) * (x[j + k * m] - x_mean[j]);
                    }
                }
            }

            for (j = 0; j < m; j++)
            {
                for (i = 0; i < m; i++)
                {
                    c[i + j * m] = c[i + j * m] / (double) (n - 1);
                }
            }

            return c;
        }

        public static double[] r8mat_standardize(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_STANDARDIZE standardizes an R8MAT.
            //
            //  Discussion:
            //
            //    The output array will have columns of 0 mean and unit standard deviation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double X[M*N], the array to be standardized.
            //
            //    Output, double R8MAT_STANDARDIZE[M*N], the standardized array.
            //
        {
            int i;
            int j;
            double[] mu;
            double[] sigma;
            double[] xs;

            mu = r8mat_mean_columns(m, n, x);
            sigma = r8mat_std_columns(m, n, x);

            xs = new double[m * n];

            for (j = 0; j < n; j++)
            {
                if (sigma[j] != 0.0)
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = (x[i + j * m] - mu[j]) / sigma[j];
                    }
                }
                else
                {
                    for (i = 0; i < m; i++)
                    {
                        xs[i + j * m] = 0.0;
                    }
                }
            }

            return xs;
        }

        public static double[] r8mat_std_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_STD_COLUMNS returns the column standard deviation of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_STD_COLUMNS[N], the column stds.
            //
        {
            int i;
            int j;
            double mean;
            double std;
            double[] std_columns;

            std_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                mean = 0.0;
                for (i = 0; i < m; i++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) m;
                std = 0.0;
                for (i = 0; i < m; i++)
                {
                    std = std + Math.Pow(a[i + j * m] - mean, 2);
                }

                std = Math.Sqrt(std / (double) (m - 1));
                std_columns[j] = std;
            }

            return std_columns;
        }

        public static double[] r8mat_std_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_STD_ROWS returns the row standard deviations of an R8MAT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_STD_ROWS[M], the row stds.
            //
        {
            int i;
            int j;
            double mean;
            double[] std_rows;
            double std;

            std_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                mean = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) n;
                std = 0.0;
                for (j = 0; j < n; j++)
                {
                    std = std + Math.Pow(a[i + j * m] - mean, 2);
                }

                std = Math.Sqrt(std / (double) (n - 1));
                std_rows[i] = std;
            }

            return std_rows;
        }

        public static double[] r8mat_variance_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_VARIANCE_COLUMNS returns the column variances of an R8MAT.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_VARIANCE_COLUMNS[N], the column variances.
            //
        {
            int i;
            int j;
            double mean;
            double variance;
            double[] variance_columns;

            variance_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                mean = 0.0;
                for (i = 0; i < m; i++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) m;
                variance = 0.0;
                for (i = 0; i < m; i++)
                {
                    variance = variance + Math.Pow(a[i + j * m] - mean, 2);
                }

                variance = variance / (double) (m - 1);
                variance_columns[j] = variance;
            }

            return variance_columns;
        }

        public static double[] r8mat_variance_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_VARIANCE_ROWS returns the row variances of an R8MAT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_VARIANCE_ROWS[M], the row variances.
            //
        {
            int i;
            int j;
            double mean;
            double[] variance_rows;
            double variance;

            variance_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                mean = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) n;
                variance = 0.0;
                for (j = 0; j < n; j++)
                {
                    variance = variance + Math.Pow(a[i + j * m] - mean, 2);
                }

                variance = variance / (double) (n - 1);
                variance_rows[i] = variance;
            }

            return variance_rows;
        }

    }
}