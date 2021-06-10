using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int[] i4row_max(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MAX returns the maximums of an I4ROW.
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
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, int I4ROW_AMAX[M], the maximums of the rows.
            //
        {
            int i;
            int j;
            int[] amax = new int[m];

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

        public static double[] i4row_mean(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MEAN returns the means of an I4ROW.
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
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, double I4ROW_MEAN[M], the means, or averages, of the rows.
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
                    mean[i] = mean[i] + (double) a[i + j * m];
                }

                mean[i] = mean[i] / (double) (n);
            }

            return mean;
        }

        public static int[] i4row_min(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MIN returns the minimums of an I4ROW.
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
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, int I4ROW_AMIN[M], the minimums of the rows.
            //
        {
            int i;
            int j;
            int[] amin = new int[m];

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

        public static double[] i4row_variance(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_VARIANCE returns the variances of an I4ROW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array whose variances are desired.
            //
            //    Output, double I4ROW_VARIANCE[M], the variances of the rows.
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
                    mean = mean + (double) a[i + j * m];
                }

                mean = mean / (double) (n);

                variance[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    variance[i] = variance[i] + Math.Pow(((double) a[i + j * m] - mean), 2);
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
        
    }
}