using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt.Sequence;

public static class Tribonacci
{
    public static int tribonacci_direct(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_direct() computes the N-th Tribonacci number directly.
        //
        //  Example:
        //
        //     N   T
        //    --  --
        //     1   0
        //     2   0
        //     3   1
        //     4   1
        //     5   2
        //     6   4
        //     7   7
        //     8  13
        //     9  24
        //    10  44
        //    11  81
        //    12 149
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int N: the index of the number to compute.
        //    N should be positive.
        //
        //  Output:
        //
        //    int TRIBONACCI_DIRECT: the value of the N-th number.
        //
    {
        double alpha = 0;
        Complex beta = new();
        Complex gamma = new();

        tribonacci_roots(ref alpha, ref beta, ref gamma);

        int t = n switch
        {
            <= 0 => 0,
            _ => (int) Math.Round((Complex.Pow(alpha, n) / (-Complex.Pow(alpha, 2) + 4.0 * alpha - 1.0) +
                                   Complex.Pow(beta, n) / (-Complex.Pow(beta, 2) + 4.0 * beta - 1.0) +
                                   Complex.Pow(gamma, n) / (-Complex.Pow(gamma, 2) + 4.0 * gamma - 1.0)).Real)
        };

        return t;
    }

    public static int[] tribonacci_recursive(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_recursive() computes the first N Tribonacci numbers.
        //
        //  Recursion:
        //
        //    F(1) = 0
        //    F(2) = 0
        //    F(3) = 1
        //
        //    F(N) = F(N-1) + F(N-2) + F(N-3)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int N, the highest number to compute.
        //
        //  Output:
        //
        //    int TRIBONACCI_RECURSIVE[N], the first N Tribonacci numbers.
        //
    {
        int i;

        int[] f = new int[n];

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case 0:
                case 1:
                    f[i] = 0;
                    break;
                case 2:
                    f[i] = 1;
                    break;
                default:
                    f[i] = f[i - 1] + f[i - 2] + f[i - 3];
                    break;
            }
        }

        return f;
    }

    public static void tribonacci_roots(ref double alpha, ref Complex beta,
            ref Complex gamma )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    tribonacci_roots() returns the Tribonacci roots.
        //
        //  Discussion:
        //
        //    The Nth Tribonacci number is defined by:
        //      T(N) = T(N-1) + T(N-2) + T(N-3)
        //    with
        //      T(1) = 0, T(2) = 0, T(3) = 1.
        //
        //    The related polynomial equation
        //      x^3 - x^2 - x - 1 = 0
        //
        //     ALPHA, BETA, and GAMMA are the roots of this equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2021
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    W R Spickerman,
        //    Binet's formula for the Tribonacci sequence,
        //    Fibonacci Quarterly, 
        //    Volume 20, Number 2, pages 118-120, May 1982.
        //
        //  Output:
        //
        //    double &ALPHA, complex <double> &BETA, complex <double> &GAMMA, the roots.
        //
    {
        double rho = typeMethods.r8_cube_root(19.0 + 3.0 * Math.Sqrt(33.0));
        double tau = typeMethods.r8_cube_root(19.0 - 3.0 * Math.Sqrt(33.0));

        double a = (2.0 - rho - tau) / 6.0;
        double b = Math.Sqrt(3.0) * (rho - tau) / 6.0;

        alpha = (1.0 + rho + tau) / 3.0;
        beta = new Complex(a, +b);
        gamma = new Complex(a, -b);
    }
}