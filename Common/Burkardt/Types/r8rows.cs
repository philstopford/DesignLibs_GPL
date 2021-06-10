using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8row_max(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ROW_MAX returns the maximums of an R8ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MAX =
            //      3
            //      7
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double A[M*N], the array to be examined.
            //
            //    Output, double R8ROW_MAX[M], the maximums of the rows.
            //
        {
            int i;
            int j;
            double[] amax = new double[m];

            for (i = 0; i < m; i++)
            {
                amax[i] = a[i + 0 * m];

                for (j = 1; j < n; j++)
                {
                    if (amax[i] < a[i + j * m])
                    {
                        amax[i] = a[i + j * m];
                    }
                }
            }

            return amax;
        }

        public static double[] r8row_mean(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ROW_MEAN returns the means of an R8ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MEAN =
            //      2
            //      5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double A[M*N], the array to be examined.
            //
            //    Output, double R8ROW_MEAN[M], the means, or averages, of the rows.
            //
        {
            int i;
            int j;
            double[] mean = new double[m];

            for (i = 0; i < m; i++)
            {
                mean[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean[i] = mean[i] + a[i + j * m];
                }

                mean[i] = mean[i] / (double) (n);
            }

            return mean;
        }

        public static double[] r8row_min(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ROW_MIN returns the minimums of an R8ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MIN =
            //      1
            //      2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double A[M*N], the array to be examined.
            //
            //    Output, double R8ROW_MIN[M], the minimums of the rows.
            //
        {
            int i;
            int j;
            double[] amin = new double[m];

            for (i = 0; i < m; i++)
            {
                amin[i] = a[i + 0 * m];
                for (j = 1; j < n; j++)
                {
                    if (a[i + j * m] < amin[i])
                    {
                        amin[i] = a[i + j * m];
                    }
                }
            }

            return amin;
        }

        public static double[] r8row_variance(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ROW_VARIANCE returns the variances of an R8ROW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double A[M*N], the array whose variances are desired.
            //
            //    Output, double R8ROW_VARIANCE[M], the variances of the rows.
            //
        {
            int i;
            int j;
            double mean;
            double[] variance = new double[m];

            for (i = 0; i < m; i++)
            {
                mean = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) (n);

                variance[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    variance[i] = variance[i] + Math.Pow((a[i + j * m] - mean), 2);
                }

                if (1 < n)
                {
                    variance[i] = variance[i] / (double) (n - 1);
                }
                else
                {
                    variance[i] = 0.0;
                }

            }

            return variance;
        }
        
        public static double[] r8rows_to_r8mat ( int m, int n, double[] r8rows )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8ROWS_TO_R8MAT converts a row-major vector to an R8MAT.
            //
            //  Discussion:
            //
            //    This function allows me to declare a vector of the right type and length,
            //    fill it with data that I can display row-wise, and then have the
            //    data copied into a column-wise doubly dimensioned array array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 September 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, double R8ROWS[M*N], the rowwise data.
            //
            //    Output, double R8ROWS_TO_R8MAT[M*N], the doubly-dimensioned columnwise data.
            //
        {
            int i;
            double[] r8mat;
            int j;
            int k;

            r8mat = new double[m*n];

            k = 0;
            for ( i = 0; i < m; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    r8mat[i+j*m] = r8rows[k];
                    k = k + 1;
                }
            }

            return r8mat;
        }
    }
}