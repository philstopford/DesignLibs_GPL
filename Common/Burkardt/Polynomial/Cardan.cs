using System;

namespace Burkardt.PolynomialNS;

public static class Cardan
{
    public static double[] cardan_poly(int n, double x, double s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDAN_POLY evaluates the Cardan polynomials.
        //
        //  First terms:
        //
        //     N  C(N,S,X)
        //
        //     0  2
        //     1  X
        //     2  X^2  -  2 S
        //     3  X^3  -  3 S X
        //     4  X^4  -  4 S X^2 +  2 S^2
        //     5  X^5  -  5 S X^3 +  5 S^2 X
        //     6  X^6  -  6 S X^4 +  9 S^2 X^2 -  2 S^3
        //     7  X^7  -  7 S X^5 + 14 S^2 X^3 -  7 S^3 X
        //     8  X^8  -  8 S X^6 + 20 S^2 X^4 - 16 S^3 X^2 +  2 S^4
        //     9  X^9  -  9 S X^7 + 27 S^2 X^5 - 30 S^3 X^3 +  9 S^4 X
        //    10  X^10 - 10 S X^8 + 35 S^2 X^6 - 50 S^3 X^4 + 25 S^4 X^2 -  2 S^5
        //    11  X^11 - 11 S X^9 + 44 S^2 X^7 - 77 S^3 X^5 + 55 S^4 X^3 - 11 S^5 X
        //
        //  Recursion:
        //
        //    Writing the N-th polynomial in terms of its coefficients:
        //
        //      C(N,S,X) = sum ( 0 <= I <= N ) D(N,I) * S^(N-I)/2 * X^I
        //
        //    then
        //
        //    D(0,0) = 1
        //
        //    D(1,1) = 1
        //    D(1,0) = 0
        //
        //    D(N,N) = 1
        //    D(N,K) = D(N-1,K-1) - D(N-2,K)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Osler,
        //    Cardan Polynomials and the Reduction of Radicals,
        //    Mathematics Magazine,
        //    Volume 74, Number 1, February 2001, pages 26-32.
        //
        //  Parameters:
        //
        //    Input, int N, the highest polynomial to compute.
        //
        //    Input, double X, the point at which the polynomials are to be computed.
        //
        //    Input, double S, the value of the parameter, which must be positive.
        //
        //    Output, double CARDAN_POLY[N+1], the values of the Cardan polynomials at X.
        //
    {
        double fact;
        int i;
        double s2;
        double[] v;
        double[] x2 = new double[1];

        s2 = Math.Sqrt(s);
        x2[0] = 0.5 * x / s2;

        v = Chebyshev.cheby_t_poly(1, n, x2);

        fact = 1.0;

        for (i = 0; i <= n; i++)
        {
            v[i] = 2.0 * fact * v[i];
            fact *= s2;
        }

        return v;
    }

    public static void cardan_poly_coef(int n, double s, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CARDAN_POLY_COEF computes the coefficients of the N-th Cardan polynomial.
        //
        //  First terms:
        //
        //    2
        //    0       1
        //   -2 S     0       1
        //    0      -3 S     0       1
        //    2 S^2   0      -4 S     0       1
        //    0       5 S^2   0      -5 S     0       1
        //   -2 S^3   0       9 S^2   0      -6 S     0       1
        //    0       7 S^3   0      14 S^2  0      -7 S     0       1
        //    2 S^4   0     -16 S^3   0      20 S^2   0      -8 S     0        1
        //    0       9 S^4   0     -30 S^3   0      27 S^2   0      -9 S      0     1
        //   -2 S^5   0      25 S^4   0     -50 S^3   0      35 S^2   0      -10 S   0   1
        //    0     -11 S^5   0      55 S^4   0     -77 S^3   0     +44 S^2    0   -11 S 0 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Osler,
        //    Cardan Polynomials and the Reduction of Radicals,
        //    Mathematics Magazine,
        //    Volume 74, Number 1, February 2001, pages 26-32.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial
        //
        //    Input, double S, the value of the parameter, which must be positive.
        //
        //    Output, double C[N+1], the coefficients.  C(0) is the constant term,
        //    and C(N) is the coefficient of X^N.
        //
    {
        double[] cm1;
        double[] cm2;
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return;
        }

        c[0] = 2.0;
        for (i = 1; i <= n; i++)
        {
            c[i] = 0.0;
        }

        switch (n)
        {
            case 0:
                return;
        }

        cm1 = new double[n + 1];
        cm2 = new double[n + 1];

        for (i = 0; i <= n; i++)
        {
            cm1[i] = c[i];
        }

        c[0] = 0.0;
        c[1] = 1.0;
        for (i = 2; i <= n; i++)
        {
            c[i] = 0.0;
        }

        for (i = 2; i <= n; i++)
        {

            for (j = 0; j <= i - 2; j++)
            {
                cm2[j] = cm1[j];
            }

            for (j = 0; j <= i - 1; j++)
            {
                cm1[j] = c[j];
            }

            c[0] = 0.0;
            for (j = 1; j <= i; j++)
            {
                c[j] = cm1[j - 1];
            }

            for (j = 0; j <= i - 2; j++)
            {
                c[j] -= s * cm2[j];
            }
        }
    }
}