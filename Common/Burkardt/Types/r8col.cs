using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8col_variance(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_VARIANCE returns the variances of an R8COL.
            //
            //  Discussion:
            //
            //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
            //    each of length M.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 September 2005
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
            //    Output, double R8COL_VARIANCE[N], the variances of the rows.
            //
        {
            double[] variance = new double[n];

            for (int j = 0; j < n; j++)
            {
                double mean = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean = mean + a[i + j * m];
                }

                mean = mean / (double) (m);

                variance[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    variance[j] = variance[j] + Math.Pow(a[i + j * m] - mean, 2);
                }

                if (1 < m)
                {
                    variance[j] = variance[j] / (double) (m - 1);
                }
                else
                {
                    variance[j] = 0.0;
                }
            }

            return variance;
        }

        public static double[] r8col_mean(int m, int n, double[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_MEAN returns the column means of an R8COL.
            //
            //  Discussion:
            //
            //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
            //    each of length M.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    R8COL_MEAN =
            //      1.5  4.0  5.0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A[M*N], the array to be examined.
            //
            //    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
            //
        {
            double[] mean = new double[n];

            for (int j = 0; j < n; j++)
            {
                mean[j] = 0.0;
                for (int i = 0; i < m; i++)
                {
                    mean[j] = mean[j] + a[i + j * m];
                }

                mean[j] = mean[j] / (double) (m);
            }

            return mean;
        }

        public static void r8col_separation(int m, int n, double[] a, ref double d_min, ref double d_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8COL_SEPARATION returns the "separation" of an R8COL.
            //
            //  Discussion:
            //
            //    D_MIN is the minimum distance between two columns,
            //    D_MAX is the maximum distance between two columns.
            //
            //    The distances are measured using the Loo norm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns 
            //    in the array.  If N < 2, it does not make sense to call this routine.
            //
            //    Input, double A[M*N], the array whose variances are desired.
            //
            //    Output, double &D_MIN, &D_MAX, the minimum and maximum distances.
            //
        {
            double d;
            int i;
            int j1;
            int j2;

            d_min = r8_huge();
            d_max = 0.0;

            for (j1 = 0; j1 < n; j1++)
            {
                for (j2 = j1 + 1; j2 < n; j2++)
                {
                    d = 0.0;
                    for (i = 0; i < m; i++)
                    {
                        d = Math.Max(d, Math.Abs(a[i + j1 * m] - a[i + j2 * m]));
                    }

                    d_min = Math.Min(d_min, d);
                    d_max = Math.Max(d_max, d);
                }
            }
        }
    }
}