using System;
using Burkardt.Types;

namespace Burkardt.Transform;

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

        double[] y = new double[n];

        int k = (int)Math.Log2(n);

        for (i = 1; i <= k; i++)
        {
            int l = k + 1 - i;
            int l2 = (int)Math.Pow(2, l - 1);

            typeMethods.r8vec_copy(2 * l2, x, ref y);

            int j;
            for (j = 1; j <= l2; j++)
            {
                int l3 = l2 + j;
                int jj = 2 * j - 1;
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

        double[] y = new double[n];

        int k = (int) Math.Log2(n);

        for (i = 1; i <= k; i++)
        {
            int l = (int) Math.Pow(2, i - 1);
            typeMethods.r8vec_copy(2 * l, x, ref y);
            int j;
            for (j = 1; j <= l; j++)
            {
                int lj = l + j;
                int jj = 2 * j;
                int jj1 = jj - 1;
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
        int ii;

        int k = (int) Math.Log2(n);

        x[0] /= Math.Pow(2.0, k);

        switch (k)
        {
            case >= 1:
                x[1] /= Math.Pow(2.0, k);
                break;
        }

        for (ii = 2; ii <= k; ii++)
        {
            int i = ii - 1;
            double wlk = 1.0 / Math.Pow(2.0, k - i);
            int jmin = (int) Math.Pow(2, i);
            int jmax = (int) Math.Pow(2, ii) - 1;
            int j;
            for (j = jmin; j <= jmax; j++)
            {
                x[j] *= wlk;
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

        double s = Math.Sqrt(2.0);

        double[] y = new double[n];
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
        int k = 1;
        while (k * 2 <= n)
        {
            k *= 2;
        }

        while (1 < k)
        {
            k /= 2;
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

        double s = Math.Sqrt(2.0);

        double[] y = new double[n];
        //
        //  Initialize.
        //
        for (i = 0; i < n; i++)
        {
            y[i] = 0.0;
        }

        int k = 1;
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

            k *= 2;
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

        double s = Math.Sqrt(2.0);

        double[] v = new double[m * n];

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
        int k = 1;
        while (k * 2 <= m)
        {
            k *= 2;
        }

        //
        //  Transform all columns.
        //
        while (1 < k)
        {
            k /= 2;

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
            k *= 2;
        }

        //
        //  Transform all rows.
        //
        while (1 < k)
        {
            k /= 2;

            for (j = 0; j < k; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + j * m] = (u[i + 2 * j * m] + u[i + (2 * j + 1) * m]) / s;
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

        double s = Math.Sqrt(2.0);

        double[] v = new double[m * n];

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
        int k = 1;

        while (k * 2 <= n)
        {
            for (j = 0; j < k; j++)
            {
                for (i = 0; i < m; i++)
                {
                    v[i + 2 * j * m] = (u[i + j * m] + u[i + (k + j) * m]) / s;
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

            k *= 2;
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

            k *= 2;
        }
    }
}