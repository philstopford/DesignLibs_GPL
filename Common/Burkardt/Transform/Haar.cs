using System;
using Burkardt.Types;

namespace Burkardt.Transform
{
    public static class Haar
    {
        public static void haar(int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAAR performs a Haar transform.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2011
            //
            //  Author:
            //
            //    Ken Beauchamp
            //
            //  Reference:
            //
            //    Ken Beauchamp,
            //    Walsh functions and their applications,
            //    Academic Press, 1975,
            //    ISBN: 0-12-084050-2,
            //    LC: QA404.5.B33.
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in X.
            //    N must be a power of 2.
            //
            //    Input/output, double X[N], the data to be transformed.
            //
        {
            int i;
            int j;
            int jj;
            int k;
            int l;
            int l2;
            int l3;
            double[] y;

            y = new double[n];

            k = (int)Math.Log2(n);

            for (i = 1; i <= k; i++)
            {
                l = k + 1 - i;
                l2 = (int)Math.Pow(2, l - 1);

                typeMethods.r8vec_copy(2 * l2, x, ref y);

                for (j = 1; j <= l2; j++)
                {
                    l3 = l2 + j;
                    jj = 2 * j - 1;
                    x[j - 1] = y[jj - 1] + y[jj];
                    x[l3 - 1] = y[jj - 1] - y[jj];
                }
            }
        }

        public static void haarin(int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAARIN inverts a Haar transform.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2011
            //
            //  Author:
            //
            //    Ken Beauchamp
            //
            //  Reference:
            //
            //    Ken Beauchamp,
            //    Walsh functions and their applications,
            //    Academic Press, 1975,
            //    ISBN: 0-12-084050-2,
            //    LC: QA404.5.B33.
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in X.
            //    N must be a power of 2.
            //
            //    Input/output, double X[N], the data to be transformed.
            //
        {
            int i;
            int j;
            int jj;
            int jj1;
            int k;
            int l;
            int lj;
            double[] y;

            y = new double[n];

            k = (int) Math.Log2(n);

            for (i = 1; i <= k; i++)
            {
                l = (int) Math.Pow(2, i - 1);
                typeMethods.r8vec_copy(2 * l, x, ref y);
                for (j = 1; j <= l; j++)
                {
                    lj = l + j;
                    jj = 2 * j;
                    jj1 = jj - 1;
                    x[jj - 1] = y[j - 1] - y[lj - 1];
                    x[jj1 - 1] = y[j - 1] + y[lj - 1];
                }
            }
        }

        public static void hnorm(int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HNORM computes normalization factors for a forward or inverse Haar transform.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2011
            //
            //  Author:
            //
            //    Ken Beauchamp
            //
            //  Reference:
            //
            //    Ken Beauchamp,
            //    Walsh functions and their applications,
            //    Academic Press, 1975,
            //    ISBN: 0-12-084050-2,
            //    LC: QA404.5.B33.
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in X.
            //    N must be a power of 2.
            //
            //    Input/output, double X[N], the data to be transformed.
            //
        {
            int i;
            int ii;
            int j;
            int jmax;
            int jmin;
            int k;
            double wlk;

            k = (int) Math.Log2(n);

            x[0] = x[0] / Math.Pow(2.0, k);

            if (1 <= k)
            {
                x[1] = x[1] / Math.Pow(2.0, k);
            }

            for (ii = 2; ii <= k; ii++)
            {
                i = ii - 1;
                wlk = 1.0 / Math.Pow(2.0, k - i);
                jmin = (int) Math.Pow(2, i);
                jmax = (int) Math.Pow(2, ii) - 1;
                for (j = jmin; j <= jmax; j++)
                {
                    x[j] = x[j] * wlk;
                }
            }
        }

        public static void haar_1d(int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAAR_1D computes the Haar transform of a vector.
            //
            //  Discussion:
            //
            //    For the classical Haar transform, N should be a power of 2.
            //    However, this is not required here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vector.
            //
            //    Input/output, double X[N], on input, the vector to be transformed.
            //    On output, the transformed vector.
            //
        {
            int i;
            int k;
            double s;
            double[] y;

            s = Math.Sqrt(2.0);

            y = new double[n];
            //
            //  Initialize.
            //
            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            //
            //  Determine K, the largest power of 2 such that K <= N.
            //
            k = 1;
            while (k * 2 <= n)
            {
                k = k * 2;
            }

            while (1 < k)
            {
                k = k / 2;
                for (i = 0; i < k; i++)
                {
                    y[i] = (x[2 * i] + x[2 * i + 1]) / s;
                    y[i + k] = (x[2 * i] - x[2 * i + 1]) / s;
                }

                for (i = 0; i < k * 2; i++)
                {
                    x[i] = y[i];
                }
            }
        }

        public static void haar_1d_inverse(int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAAR_1D_INVERSE computes the inverse Haar transform of a vector.
            //
            //  Discussion:
            //
            //    For the classical Haar transform, N should be a power of 2.
            //    However, this is not required here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vector.  
            //
            //    Input/output, double X[N], on input, the vector to be transformed.
            //    On output, the transformed vector.
            //
        {
            int i;
            int k;
            double s;
            double[] y;

            s = Math.Sqrt(2.0);

            y = new double[n];
            //
            //  Initialize.
            //
            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            k = 1;
            while (k * 2 <= n)
            {
                for (i = 0; i < k; i++)
                {
                    y[2 * i] = (x[i] + x[i + k]) / s;
                    y[2 * i + 1] = (x[i] - x[i + k]) / s;
                }

                for (i = 0; i < k * 2; i++)
                {
                    x[i] = y[i];
                }

                k = k * 2;
            }
        }

        public static void haar_2d(int m, int n, ref double[] u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAAR_2D computes the Haar transform of an array.
            //
            //  Discussion:
            //
            //    For the classical Haar transform, M and N should be a power of 2.
            //    However, this is not required here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimensions of the array.
            //
            //    Input/output, double U[M*N], the array to be transformed.
            //
        {
            int i;
            int j;
            int k;
            double s;
            double[] v;

            s = Math.Sqrt(2.0);

            v = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = u[i + j * m];
                }
            }

            //
            //  Determine K, the largest power of 2 such that K <= M.
            //
            k = 1;
            while (k * 2 <= m)
            {
                k = k * 2;
            }

            //
            //  Transform all columns.
            //
            while (1 < k)
            {
                k = k / 2;

                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < k; i++)
                    {
                        v[i + j * m] = (u[2 * i + j * m] + u[2 * i + 1 + j * m]) / s;
                        v[k + i + j * m] = (u[2 * i + j * m] - u[2 * i + 1 + j * m]) / s;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < 2 * k; i++)
                    {
                        u[i + j * m] = v[i + j * m];
                    }
                }
            }

            //
            //  Determine K, the largest power of 2 such that K <= N.
            //
            k = 1;
            while (k * 2 <= n)
            {
                k = k * 2;
            }

            //
            //  Transform all rows.
            //
            while (1 < k)
            {
                k = k / 2;

                for (j = 0; j < k; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        v[i + (j) * m] = (u[i + 2 * j * m] + u[i + (2 * j + 1) * m]) / s;
                        v[i + (k + j) * m] = (u[i + 2 * j * m] - u[i + (2 * j + 1) * m]) / s;
                    }
                }

                for (j = 0; j < 2 * k; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        u[i + j * m] = v[i + j * m];
                    }
                }
            }
        }

        public static void haar_2d_inverse(int m, int n, ref double[] u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HAAR_2D_INVERSE inverts the Haar transform of an array.
            //
            //  Discussion:
            //
            //    For the classical Haar transform, M and N should be a power of 2.
            //    However, this is not required here.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 March 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the dimensions of the array.
            //
            //    Input/output, double U[M*N], the array to be transformed.
            //
        {
            int i;
            int j;
            int k;
            double s;
            double[] v;

            s = Math.Sqrt(2.0);

            v = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = u[i + j * m];
                }
            }

            //
            //  Inverse transform of all rows.
            //
            k = 1;

            while (k * 2 <= n)
            {
                for (j = 0; j < k; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        v[i + (2 * j) * m] = (u[i + j * m] + u[i + (k + j) * m]) / s;
                        v[i + (2 * j + 1) * m] = (u[i + j * m] - u[i + (k + j) * m]) / s;
                    }
                }

                for (j = 0; j < 2 * k; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        u[i + j * m] = v[i + j * m];
                    }
                }

                k = k * 2;
            }

            //
            //  Inverse transform of all columns.
            //
            k = 1;

            while (k * 2 <= m)
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < k; i++)
                    {
                        v[2 * i + j * m] = (u[i + j * m] + u[k + i + j * m]) / s;
                        v[2 * i + 1 + j * m] = (u[i + j * m] - u[k + i + j * m]) / s;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < 2 * k; i++)
                    {
                        u[i + j * m] = v[i + j * m];
                    }
                }

                k = k * 2;
            }
        }
    }
}