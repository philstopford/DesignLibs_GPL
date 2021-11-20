using System;

namespace Burkardt.MatrixNS;

public static class ToeplitzMatrix
{
    public static double[] t_cholesky_lower(int n, double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The first row of the Toeplitz matrix A is supplied.
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = L * L'.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double T[N], the first row.
        //
        //    Output, double T_CHOLESKY_LOWER[N*N], the lower Cholesky factor.
        //
    {
        int i;
        int j;

        double[] g = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            g[0 + j * 2] = t[j];
        }

        g[1 + 0 * 2] = 0.0;
        for (j = 1; j < n; j++)
        {
            g[1 + j * 2] = t[j];
        }

        double[] l = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                l[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            l[i + 0 * n] = g[0 + i * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * 2] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));
            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                l[j + i * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }

        return l;
    }

    public static double[] t_cholesky_upper(int n, double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The first row of the Toeplitz matrix A is supplied.
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = R' * R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double T[N}, the first row.
        //
        //    Output, double T_CHOLESKY_UPPER[N*N], the upper Cholesky factor.
        //
    {
        int i;
        int j;

        double[] g = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            g[0 + j * 2] = t[j];
        }

        g[1 + 0 * 2] = 0.0;
        for (j = 1; j < n; j++)
        {
            g[1 + j * 2] = t[j];
        }

        double[] r = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                r[i + j * n] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            r[0 + j * n] = g[0 + j * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * n] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));

            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                r[i + j * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }

        return r;
    }

    public static double[] toep_cholesky_lower(int n, double[] g)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEP_CHOLESKY_LOWER: lower Cholesky factor of a compressed Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The Toeplitz matrix A is supplied in a compressed form G.
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = L * L'.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double G[2*N], the compressed Toeplitz matrix.
        //    G(1,1:N) contains the first row.
        //    G(2,2:N) contains the first column.
        //
        //    Output, double TOEP_CHOLESKY_LOWER[N*N], the lower Cholesky factor.
        //
    {
        int i;
        int j;

        double[] l = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                l[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            l[i + 0 * n] = g[0 + i * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * 2] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));
            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                l[j + i * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }

        return l;
    }

    public static double[] toep_cholesky_upper(int n, double[] g)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEP_CHOLESKY_UPPER: upper Cholesky factor of a compressed Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The Toeplitz matrix A is supplied in a compressed form G.
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = R' * R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double G[2*N}, the compressed Toeplitz matrix.
        //    G(1,1:N) contains the first row.
        //    G(2,2:N) contains the first column.
        //
        //    Output, double TOEP_CHOLESKY_UPPER[N*N], the upper Cholesky factor.
        //
    {
        int i;
        int j;

        double[] r = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                r[i + j * n] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            r[0 + j * n] = g[0 + j * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * n] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));

            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                r[i + j * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }

        return r;
    }

    public static double[] toeplitz_cholesky_lower(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEPLITZ_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = L * L'.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the Toeplitz matrix.
        //
        //    Output, double TOEPLITZ_CHOLESKY_LOWER[N*N], the lower Cholesky factor.
        //
    {
        int i;
        int j;

        double[] l = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                l[i + j * n] = 0.0;
            }
        }

        double[] g = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            g[0 + j * 2] = a[0 + j * n];
        }

        g[1 + 0 * 2] = 0.0;
        for (j = 1; j < n; j++)
        {
            g[1 + j * 2] = a[j + 0 * n];
        }

        for (i = 0; i < n; i++)
        {
            l[i + 0 * n] = g[0 + i * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * 2] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));

            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                l[j + i * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }
            
        return l;
    }

    public static double[] toeplitz_cholesky_upper(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TOEPLITZ_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
        //
        //  Discussion:
        //
        //    The Toeplitz matrix must be positive semi-definite.
        //
        //    After factorization, A = R' * R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Michael Stewart,
        //    Cholesky factorization of semi-definite Toeplitz matrices.
        //    Linear Algebra and its Applications,
        //    Volume 254, pages 497-525, 1997.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the Toeplitz matrix.
        //
        //    Output, double TOEPLITZ_CHOLESKY_UPPER[N*N], the upper Cholesky factor.
        //
    {
        int i;
        int j;

        double[] r = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                r[i + j * n] = 0.0;
            }
        }

        double[] g = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            g[0 + j * 2] = a[0 + j * n];
        }

        g[1 + 0 * 2] = 0.0;
        for (j = 1; j < n; j++)
        {
            g[1 + j * 2] = a[j + 0 * n];
        }

        for (j = 0; j < n; j++)
        {
            r[0 + j * n] = g[0 + j * 2];
        }

        for (j = n - 1; 1 <= j; j--)
        {
            g[0 + j * 2] = g[0 + (j - 1) * 2];
        }

        g[0 + 0 * 2] = 0.0;

        for (i = 1; i < n; i++)
        {
            double rho = -g[1 + i * 2] / g[0 + i * 2];
            double div = Math.Sqrt((1.0 - rho) * (1.0 + rho));
            for (j = i; j < n; j++)
            {
                double g1j = g[0 + j * 2];
                double g2j = g[1 + j * 2];
                g[0 + j * 2] = (g1j + rho * g2j) / div;
                g[1 + j * 2] = (rho * g1j + g2j) / div;
            }

            for (j = i; j < n; j++)
            {
                r[i + j * n] = g[0 + j * 2];
            }

            for (j = n - 1; i < j; j--)
            {
                g[0 + j * 2] = g[0 + (j - 1) * 2];
            }

            g[0 + i * 2] = 0.0;
        }
            
        return r;
    }

}