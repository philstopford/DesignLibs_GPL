using System;
using Burkardt.Probability;
using Burkardt.Sequence;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class Moments
{
    public static void moment_method(int n, double[] moment, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENT_METHOD computes a quadrature rule by the method of moments.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Gene Golub, John Welsch,
        //    Calculation of Gaussian Quadrature Rules,
        //    Mathematics of Computation,
        //    Volume 23, Number 106, April 1969, pages 221-230.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the quadrature rule.
        //
        //    Input, double MOMENT[2*N+1], moments 0 through 2*N.
        //
        //    Output, double X[N], W[N], the points and weights of the quadrature rule.
        //
    {
        int flag = 0;
        int i;
        int it_num = 0;
        int j;
        int rot_num = 0;

        const bool debug = false;

        switch (debug)
        {
            case true:
                typeMethods.r8vec_print(2 * n + 1, moment, "  Moments:");
                break;
        }

        //
        //  Define the N+1 by N+1 Hankel matrix H(I,J) = moment(I+J).
        //
        double[] h = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                h[i + j * (n + 1)] = moment[i + j];
            }
        }

        switch (debug)
        {
            case true:
                typeMethods.r8mat_print(n + 1, n + 1, h, "  Hankel matrix:");
                break;
        }

        //
        //  Compute R, the upper triangular Cholesky factor of H.
        //
        double[] r = typeMethods.r8mat_cholesky_factor_upper(n + 1, h, ref flag);

        if (flag != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("MOMENT_METHOD - Fatal error!");
            Console.WriteLine("  R8MAT_CHOLESKY_FACTOR_UPPER returned FLAG = " + flag + "");
            return;
        }

        switch (debug)
        {
            case true:
                typeMethods.r8mat_print(n + 1, n + 1, r, "  Cholesky factor:");
                break;
        }

        //
        //  Compute ALPHA and BETA from R, using Golub and Welsch's formula.
        //
        double[] alpha = new double[n];

        alpha[0] = r[0 + 1 * (n + 1)] / r[0 + 0 * (n + 1)];
        for (i = 1; i < n; i++)
        {
            alpha[i] = r[i + (i + 1) * (n + 1)] / r[i + i * (n + 1)]
                       - r[i - 1 + i * (n + 1)] / r[i - 1 + (i - 1) * (n + 1)];
        }

        double[] beta = new double[n - 1];

        for (i = 0; i < n - 1; i++)
        {
            beta[i] = r[i + 1 + (i + 1) * (n + 1)] / r[i + i * (n + 1)];
        }

        //
        //  Compute the points and weights from the moments.
        //
        double[] jacobi = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                jacobi[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            jacobi[i + i * n] = alpha[i];
        }

        for (i = 0; i < n - 1; i++)
        {
            jacobi[i + (i + 1) * n] = beta[i];
            jacobi[i + 1 + i * n] = beta[i];
        }

        switch (debug)
        {
            case true:
                typeMethods.r8mat_print(n, n, jacobi, "  The Jacobi matrix:");
                break;
        }

        //
        //  Get the eigendecomposition of the Jacobi matrix.
        //
        const int it_max = 100;
        double[] v = new double[n * n];

        Jacobi.jacobi_eigenvalue(n, jacobi, it_max, ref v, ref x, ref it_num, ref rot_num);

        switch (debug)
        {
            case true:
                typeMethods.r8mat_print(n, n, v, "  Eigenvector");
                break;
        }

        for (i = 0; i < n; i++)
        {
            w[i] = moment[0] * Math.Pow(v[0 + i * n], 2);
        }
    }

    public static double[] moments_normal(int m, double mu, double sigma)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENTS_NORMAL returns moments of the standard Normal distribution.
        //
        //  Discussion:
        //
        //    pdf(x) = exp ( -((x-mu)/sigma)^2/2 ) / sigma / sqrt ( Math.PI * 2 )
        //    mu(k) = integral ( -oo < x < +oo ) x^k pdf(x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments desired.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //
        //    Output, double W(0:M-1), the weighted integrals of X^0 
        //    through X^(M-1).
        //
    {
        int k;

        double[] w = new double[m];

        for (k = 0; k < m; k++)
        {
            double t = 0.0;
            int j_hi = k / 2;
            int j;
            for (j = 0; j <= j_hi; j++)
            {
                t += typeMethods.r8_choose(k, 2 * j) * typeMethods.r8_factorial2(2 * j - 1)
                                                     * Math.Pow(sigma, 2 * j) * Math.Pow(mu, k - 2 * j);
            }

            w[k] = t;
        }

        return w;
    }

    public static double[] moments_truncated_normal_ab(int m, double mu, double sigma,
            double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENTS_TRUNCATED_NORMAL_AB: moments of truncated Normal distribution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments desired.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //
        //    Input, double A, B, the lower and upper truncation limits.
        //    A < B.
        //
        //    Output, double W(0:M-1), the weighted integrals of X^0 
        //    through X^(M-1).
        //
    {
        int order;

        double[] w = new double[m];

        for (order = 0; order < m; order++)
        {
            w[order] = Truncated.truncated_normal_ab_moment(order, mu, sigma, a, b);
        }

        return w;
    }

    public static double[] moments_truncated_normal_a(int m, double mu, double sigma,
            double a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENTS_TRUNCATED_NORMAL_A: moments of lower truncated Normal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments desired.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //
        //    Input, double A, the lower truncation limit.
        //    A < B.
        //
        //    Output, double W(0:M-1), the weighted integrals of X^0 
        //    through X^(M-1).
        //
    {
        int order;

        double[] w = new double[m];

        for (order = 0; order < m; order++)
        {
            w[order] = Truncated.truncated_normal_a_moment(order, mu, sigma, a);
        }

        return w;
    }

    public static double[] moments_truncated_normal_b(int m, double mu, double sigma,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOMENTS_TRUNCATED_NORMAL_B: moments of upper truncated Normal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments desired.
        //
        //    Input, double MU, SIGMA, the mean and standard deviation.
        //
        //    Input, double B, the upper truncation limit.
        //    A < B.
        //
        //    Output, double W(0:M-1), the weighted integrals of X^0 
        //    through X^(M-1).
        //
    {
        int order;

        double[] w = new double[m];

        for (order = 0; order < m; order++)
        {
            w[order] = Truncated.truncated_normal_b_moment(order, mu, sigma, b);
        }

        return w;
    }
}