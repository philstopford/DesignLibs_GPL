using System;
using Burkardt.Types;

namespace Burkardt.Weight;

public static class WM
{
    public static double[] wm(int m, int kind, double alpha, double beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WM evaluates the first M moments of classical weight functions.
        //
        //  Discussion:
        //
        //    W(K) = Integral ( A <= X <= B ) X**(K-1) * W(X) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Sylvan Elhay, Jaroslav Kautsky,
        //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
        //    Interpolatory Quadrature,
        //    ACM Transactions on Mathematical Software,
        //    Volume 13, Number 4, December 1987, pages 399-415.
        //
        //  Parameters:
        //
        //    Input, int M, the number of moments to evaluate.
        //
        //    Input, int KIND, the rule.
        //    1, Legendre,             (a,b)       1.0
        //    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
        //    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
        //    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
        //    5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
        //    6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
        //    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
        //    8, Rational,             (a,+oo)     (x-a)^alpha*(x+b)^beta
        //    9, Chebyshev Type 2,     (a,b)       ((b-x)*(x-a))^(+0.5)
        //
        //    Input, double ALPHA, the value of Alpha, if needed.
        //
        //    Input, double BETA, the value of Beta, if needed.
        //
        //    Output, double WM[M], the first M moments.
        //
    {
        double als;
        int k;
            
        double rk;

        PARCHK.parchk(kind, m, alpha, beta);

        double[] w = new double[m];

        for (k = 2; k <= m; k += 2)
        {
            w[k - 1] = 0.0;
        }

        switch (kind)
        {
            case 1:
            {
                for (k = 1; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = 2.0 / rk;
                }

                break;
            }
            case 2:
            {
                w[0] = Math.PI;
                for (k = 3; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = w[k - 3] * (rk - 2.0) / (rk - 1.0);
                }

                break;
            }
            case 3:
            {
                w[0] = Math.Sqrt(Math.PI) * typeMethods.r8_gamma(alpha + 1.0)
                       / typeMethods.r8_gamma(alpha + 3.0 / 2.0);

                for (k = 3; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = w[k - 3] * (rk - 2.0) / (2.0 * alpha + rk);
                }

                break;
            }
            case 4:
            {
                als = alpha + beta + 1.0;
                w[0] = Math.Pow(2.0, als) * typeMethods.r8_gamma(alpha + 1.0)
                    / typeMethods.r8_gamma(als + 1.0) * typeMethods.r8_gamma(beta + 1.0);

                for (k = 2; k <= m; k++)
                {
                    double sum = 0.0;
                    double trm = 1.0;
                    rk = k;

                    int i;
                    for (i = 0; i <= (k - 2) / 2; i++)
                    {
                        double tmpa = trm;
                        int ja;
                        for (ja = 1; ja <= 2 * i; ja++)
                        {
                            tmpa = tmpa * (alpha + ja) / (als + ja);
                        }

                        int jb;
                        for (jb = 1; jb <= k - 2 * i - 1; jb++)
                        {
                            tmpa = tmpa * (beta + jb) / (als + 2 * i + jb);
                        }

                        tmpa = tmpa / (2 * i + 1.0) *
                               (2 * i * (beta + alpha) + beta - (rk - 1.0) * alpha)
                               / (beta + rk - 2 * i - 1.0);
                        sum += tmpa;

                        trm = trm * (rk - 2 * i - 1.0)
                            / (2 * i + 1.0) * (rk - 2 * i - 2.0) / (2 * i + 2.0);
                    }

                    if (k % 2 != 0)
                    {
                        double tmpb = 1.0;
                        for (i = 1; i <= k - 1; i++)
                        {
                            tmpb = tmpb * (alpha + i) / (als + i);
                        }

                        sum += tmpb;
                    }

                    w[k - 1] = sum * w[0];
                }

                break;
            }
            case 5:
            {
                w[0] = typeMethods.r8_gamma(alpha + 1.0);

                for (k = 2; k <= m; k++)
                {
                    rk = k;
                    w[k - 1] = (alpha + rk - 1.0) * w[k - 2];
                }

                break;
            }
            case 6:
            {
                w[0] = typeMethods.r8_gamma((alpha + 1.0) / 2.0);

                for (k = 3; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = w[k - 3] * (alpha + rk - 2.0) / 2.0;
                }

                break;
            }
            case 7:
            {
                als = alpha;
                for (k = 1; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = 2.0 / (rk + als);
                }

                break;
            }
            case 8:
            {
                w[0] = typeMethods.r8_gamma(alpha + 1.0)
                       * typeMethods.r8_gamma(-alpha - beta - 1.0)
                       / typeMethods.r8_gamma(-beta);

                for (k = 2; k <= m; k++)
                {
                    rk = k;
                    w[k - 1] = -w[k - 2] * (alpha + rk - 1.0) / (alpha + beta + rk);
                }

                break;
            }
            case 9:
            {
                w[0] = Math.PI / 2.0;
                for (k = 3; k <= m; k += 2)
                {
                    rk = k;
                    w[k - 1] = w[k - 3] * (rk - 2.0) / (rk + 1.0);
                }

                break;
            }
        }

        return w;
    }
}