using System;

namespace Burkardt.Transform
{
    public static class Haar
    {
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