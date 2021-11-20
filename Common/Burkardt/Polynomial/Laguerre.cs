using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

using QuadratureRule = Burkardt.Laguerre.QuadratureRule;

public static class Laguerre
{
    public static void gen_laguerre_poly(int n, double alpha, double x, ref double[] cx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEN_LAGUERRE_POLY evaluates generalized Laguerre polynomials.
        //
        //  Differential equation:
        //
        //    X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0
        //
        //  Recursion:
        //
        //    L(0,ALPHA,X) = 1
        //    L(1,ALPHA,X) = 1+ALPHA-X
        //
        //    L(N,ALPHA,X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA,X)
        //                   - (N-1+ALPHA) * L(N-2,ALPHA,X) ) / N
        //
        //  Restrictions:
        //
        //    -1 < ALPHA
        //
        //  Special values:
        //
        //    For ALPHA = 0, the generalized Laguerre polynomial L(N,ALPHA,X)
        //    is equal to the Laguerre polynomial L(N,X).
        //
        //    For ALPHA integral, the generalized Laguerre polynomial
        //    L(N,ALPHA,X) equals the associated Laguerre polynomial L(N,ALPHA,X).
        //
        //  Norm:
        //
        //    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,ALPHA,X)^2 dX
        //    = Gamma ( N + ALPHA + 1 ) / N!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order function to compute.
        //
        //    Input, double ALPHA, the parameter.  -1 < ALPHA is required.
        //
        //    Input, double X, the point at which the functions are to be
        //    evaluated.
        //
        //    Output, double CX[N+1], the polynomials of
        //    degrees 0 through N evaluated at the point X.
        //
    {
        int i;

        switch (alpha)
        {
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("GEN_LAGUERRE_POLY - Fatal error!");
                Console.WriteLine("  The input value of ALPHA is " + alpha + "");
                Console.WriteLine("  but ALPHA must be greater than -1.");
                return;
        }

        switch (n)
        {
            case < 0:
                return;
        }

        cx[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        cx[1] = 1.0 + alpha - x;

        for (i = 2; i <= n; i++)
        {
            cx[i] = ((2 * i - 1 + alpha - x) * cx[i - 1]
                     + (-i + 1 - alpha) * cx[i - 2])
                    / i;
        }
    }

    public static void laguerre_recur(ref double p2, ref double dp2, ref double p1, double x,
            int order, double alpha, double[] b, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Output, double *P2, the value of L(ORDER)(X).
        //
        //    Output, double *DP2, the value of L'(ORDER)(X).
        //
        //    Output, double *P1, the value of L(ORDER-1)(X).
        //
        //    Input, double X, the point at which polynomials are evaluated.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
        //    Input, double ALPHA, the exponent of the X factor in the
        //    integrand.
        //
        //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
        //
    {
        int i;

        p1 = 1.0;
        double dp1 = 0.0;

        p2 = x - alpha - 1.0;
        dp2 = 1.0;

        for (i = 1; i < order; i++)
        {
            double p0 = p1;
            double dp0 = dp1;

            p1 = p2;
            dp1 = dp2;

            p2 = (x - b[i]) * p1 - c[i] * p0;
            dp2 = (x - b[i]) * dp1 + p1 - c[i] * dp0;
        }
    }

    public static void laguerre_root(ref double x, int order, double alpha, ref double dp2,
            ref double p1, double[] b, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input/output, double *X, the approximate root, which
        //    should be improved on output.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //
        //    Output, double *DP2, the value of L'(ORDER)(X).
        //
        //    Output, double *P1, the value of L(ORDER-1)(X).
        //
        //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
        //
    {
        double p2 = 0;
        int step;
        const int step_max = 10;

        double eps = typeMethods.r8_epsilon();

        for (step = 1; step <= step_max; step++)
        {
            laguerre_recur(ref p2, ref dp2, ref p1, x, order, alpha, b, c);

            double d = p2 / dp2;
            x -= d;

            if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
            {
                break;
            }
        }
    }

    public static void laguerre_associated(int n, int m, double x, ref double[] cx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ASSOCIATED evaluates the associated Laguerre polynomials L(N,M,X) at X.
        //
        //  Differential equation:
        //
        //    X Y'' + (M+1-X) Y' + (N-M) Y = 0
        //
        //  First terms:
        //
        //    M = 0
        //
        //    L(0,0,X) =   1
        //    L(1,0,X) =  -X    +  1
        //    L(2,0,X) =   X^2 -  4 X     +  2
        //    L(3,0,X) =  -X^3 +  9 X^2 -  18 X    +    6
        //    L(4,0,X) =   X^4 - 16 X^3 +  72 X^2 -   96 X +      24
        //    L(5,0,X) =  -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120
        //    L(6,0,X) =   X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720
        //
        //    M = 1
        //
        //    L(0,1,X) =    0
        //    L(1,1,X) =   -1,
        //    L(2,1,X) =    2 X - 4,
        //    L(3,1,X) =   -3 X^2 + 18 X - 18,
        //    L(4,1,X) =    4 X^3 - 48 X^2 + 144 X - 96
        //
        //    M = 2
        //
        //    L(0,2,X) =    0
        //    L(1,2,X) =    0,
        //    L(2,2,X) =    2,
        //    L(3,2,X) =   -6 X + 18,
        //    L(4,2,X) =   12 X^2 - 96 X + 144
        //
        //    M = 3
        //
        //    L(0,3,X) =    0
        //    L(1,3,X) =    0,
        //    L(2,3,X) =    0,
        //    L(3,3,X) =   -6,
        //    L(4,3,X) =   24 X - 96
        //
        //    M = 4
        //
        //    L(0,4,X) =    0
        //    L(1,4,X) =    0
        //    L(2,4,X) =    0
        //    L(3,4,X) =    0
        //    L(4,4,X) =   24
        //
        //  Recursion:
        //
        //    if N = 0:
        //
        //      L(N,M,X)   = 0
        //
        //    if N = 1:
        //
        //      L(N,M,X)   = (M+1-X)
        //
        //    if N => 2:
        //
        //      L(N,M,X)   = ( (M+2*N-1-X) * L(N-1,M,X)
        //                  +   (1-M-N)     * L(N-2,M,X) ) / N
        //
        //  Special values:
        //
        //    For M = 0, the associated Laguerre polynomials L(N,M,X) are equal
        //    to the Laguerre polynomials L(N,X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, int M, the parameter.  M must be nonnegative.
        //
        //    Input, double X, the point at which the polynomials are to be evaluated.
        //
        //    Output, double CX[N+1], the associated Laguerre polynomials of
        //    degrees 0 through N evaluated at the point X.
        //
    {
        int i;

        switch (m)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("LAGUERRE_ASSOCIATED - Fatal error!");
                Console.WriteLine("  Input value of M = " + m + "");
                Console.WriteLine("  but M must be nonnegative.");
                return;
        }

        switch (n)
        {
            case < 0:
                return;
        }

        cx[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        cx[1] = m + 1 - x;

        for (i = 2; i <= n; i++)
        {
            cx[i] = ((2 * i + m - 1 - x) * cx[i - 1]
                     + (-i - m + 1) * cx[i - 2])
                    / i;
        }

    }

    public static void laguerre_poly(int n, double x, ref double[] cx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_POLY evaluates the Laguerre polynomials at X.
        //
        //  Differential equation:
        //
        //    X * Y'' + (1-X) * Y' + N * Y = 0
        //
        //  First terms:
        //
        //      1
        //     -X    +  1
        //   (  X^2 -  4 X     +  2 ) / 2
        //   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
        //   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
        //   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
        //   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
        //   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
        //      + 52920 X^2 - 35280 X + 5040 ) / 5040
        //
        //  Recursion:
        //
        //    L(0,X) = 1,
        //    L(1,X) = 1-X,
        //    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
        //
        //  Orthogonality:
        //
        //    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
        //    = 0 if N /= M
        //    = 1 if N == M
        //
        //  Special values:
        //
        //    L(N,0) = 1.
        //
        //  Relations:
        //
        //    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X, the point at which the polynomials are to be evaluated.
        //
        //    Output, double CX[N+1], the Laguerre polynomials of degree 0 through
        //    N evaluated at the point X.
        //
    {
        int i;

        switch (n)
        {
            case < 0:
                return;
        }

        cx[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        cx[1] = 1.0 - x;

        for (i = 2; i <= n; i++)
        {
            cx[i] = ((2 * i - 1 - x) * cx[i - 1]
                     + (-i + 1) * cx[i - 2])
                    / i;

        }
    }

    public static void laguerre_poly_coef(int n, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_POLY_COEF evaluates the Laguerre polynomial coefficients.
        //
        //  First terms:
        //
        //    0: 1
        //    1: 1  -1
        //    2: 1  -2  1/2
        //    3: 1  -3  3/2  1/6
        //    4: 1  -4  4   -2/3  1/24
        //    5: 1  -5  5   -5/3  5/24  -1/120
        //
        //  Recursion:
        //
        //    L(0) = ( 1,  0, 0, ..., 0 )
        //    L(1) = ( 1, -1, 0, ..., 0 )
        //    L(N) = (2*N-1-X) * L(N-1) - (N-1) * L(N-2) / N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Output, double C[(N+1)*(N+1)], the coefficients of the Laguerre 
        //    polynomials of degree 0 through N.  Each polynomial is stored as a row.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return;
        }

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        for (i = 0; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = 1.0;
        }

        switch (n)
        {
            case 0:
                return;
        }

        c[1 + 1 * (n + 1)] = -1.0;

        for (i = 2; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                c[i + j * (n + 1)] = (
                                         (2 * i - 1) * c[i - 1 + j * (n + 1)]
                                         + (-i + 1) * c[i - 2 + j * (n + 1)]
                                         - c[i - 1 + (j - 1) * (n + 1)])
                                     / i;
            }

        }
    }

    public static double[] l_exponential_product(int p, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_EXPONENTIAL_PRODUCT: exponential product table for L(n,x).
        //
        //  Discussion:
        //
        //    Let L(n,x) represent the Laguerre polynomial of degree n.  
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of exp(B*X) with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( 0 <= X < +oo ) exp(b*x) * L(i,x) * L(j,x) * exp (-x) dx
        //
        //    Because of the exponential factor, the quadrature will not be exact.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, double B, the coefficient of X in the exponential factor.
        //
        //    Output, double L_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of integrals.  
        //    TABLE(I,J) represents the weighted integral of exp(B*X) * L(i,x) * L(j,x).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];

        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = (3 * p + 4) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        QuadratureRule.l_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] l_table = l_polynomial(1, p, x_table, xIndex: +k);
            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] += w_table[k] * Math.Exp(b * x) * l_table[i] * l_table[j];
                }
            }
        }

        return table;
    }

    public static double l_integral(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_INTEGRAL evaluates a monomial integral associated with L(n,x).
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the exponent.
        //    0 <= N.
        //
        //    Output, double L_INTEGRAL, the value of the integral.
        //
    {
        double value = 0;

        value = typeMethods.r8_factorial(n);

        return value;
    }

    public static double[] l_polynomial(int m, int n, double[] x, int xIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_POLYNOMIAL evaluates the Laguerre polynomials L(n,x).
        //
        //  Differential equation:
        //
        //    X * Y'' + (1-X) * Y' + N * Y = 0
        //
        //  First terms:
        //
        //      1
        //     -X    +  1
        //   (  X^2 -  4 X     +  2 ) / 2
        //   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
        //   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
        //   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
        //   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
        //   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
        //      + 52920 X^2 - 35280 X + 5040 ) / 5040
        //
        //  Recursion:
        //
        //    L(0,X) = 1,
        //    L(1,X) = 1-X,
        //    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
        //
        //  Orthogonality:
        //
        //    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N,X) * L(M,X) dX
        //    = 0 if N /= M
        //    = 1 if N == M
        //
        //  Special values:
        //
        //    L(N,0) = 1.
        //
        //  Relations:
        //
        //    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double L_POLYNOMIAL[M*(N+1)], the function values.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = 1.0 - x[xIndex + i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = ((2 * j - 1 - x[xIndex + i]) * v[i + (j - 1) * m]
                                + (-j + 1) * v[i + (j - 2) * m])
                               / j;
            }
        }

        return v;
    }

    public static double[] l_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_POLYNOMIAL_COEFFICIENTS: coeffs for Laguerre polynomial L(n,x).
        //
        //  First terms:
        //
        //    0: 1
        //    1: 1  -1
        //    2: 1  -2  1/2
        //    3: 1  -3  3/2  1/6
        //    4: 1  -4  4   -2/3  1/24
        //    5: 1  -5  5   -5/3  5/24  -1/120
        //
        //  Recursion:
        //
        //    L(0,X) = ( 1,  0, 0, ..., 0 )
        //    L(1,X) = ( 1, -1, 0, ..., 0 )
        //    L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X) / N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Output, double L_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients
        //    of the Laguerre polynomials of degree 0 through N. 
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] c = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        for (i = 0; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = 1.0;
        }

        switch (n)
        {
            case 0:
                return c;
        }

        c[1 + 1 * (n + 1)] = -1.0;

        for (i = 2; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                c[i + j * (n + 1)] = (
                                         (2 * i - 1) * c[i - 1 + j * (n + 1)]
                                         + (-i + 1) * c[i - 2 + j * (n + 1)]
                                         - c[i - 1 + (j - 1) * (n + 1)])
                                     / i;
            }
        }

        return c;
    }

    public static void l_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial L(n,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      LaguerreL[n,x]
        //
        //  Differential equation:
        //
        //    X * Y'' + (1-X) * Y' + N * Y = 0;
        //
        //  First terms:
        //
        //      1
        //     -X    +  1
        //   (  X^2 -  4 X     +  2 ) / 2
        //   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
        //   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
        //   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
        //   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
        //   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
        //      + 52920 X^2 - 35280 X + 5040 ) / 5040
        //
        //  Recursion:
        //
        //    L(0,X) = 1,
        //    L(1,X) = 1-X,
        //    N * L(N,X) = (2*N-1-X) * L(N-1,X) - (N-1) * L(N-2,X)
        //
        //  Orthogonality:
        //
        //    Integral ( 0 <= X < oo ) exp ( - X ) * L(N,X) * L(M,X) dX
        //    = 0 if N /= M
        //    = 1 if N == M
        //
        //  Special values:
        //
        //    L(N,0) = 1.
        //
        //  Relations:
        //
        //    L(N,X) = (-1)^N / N! * exp ( x ) * (d/dx)^n ( exp ( - x ) * x^n )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the order of the polynomial.
        //
        //    Output, double &X, the point where the polynomial is evaluated.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 17;

        double[] fx_vec =
        {
            0.1000000000000000E+01,
            0.0000000000000000E+00,
            -0.5000000000000000E+00,
            -0.6666666666666667E+00,
            -0.6250000000000000E+00,
            -0.4666666666666667E+00,
            -0.2569444444444444E+00,
            -0.4047619047619048E-01,
            0.1539930555555556E+00,
            0.3097442680776014E+00,
            0.4189459325396825E+00,
            0.4801341790925124E+00,
            0.4962122235082305E+00,
            -0.4455729166666667E+00,
            0.8500000000000000E+00,
            -0.3166666666666667E+01,
            0.3433333333333333E+02
        };

        int[] n_vec =
        {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 11,
            12, 5, 5,
            5, 5
        };

        double[] x_vec =
        {
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            0.5E+00,
            3.0E+00,
            5.0E+00,
            1.0E+01
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            n = 0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

    public static double[] l_polynomial_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_POLYNOMIAL_ZEROS: zeros of the Laguerre polynomial L(n,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt.
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
        //    Input, int N, the order of the polynomial.
        //
        //    Output, double L_POLYNOMIAL_ZEROS[N], the zeros.
        //
    {
        int i;
        //
        //  Define the zero-th moment.
        //
        double zemu = 1.0;
        //
        //  Define the Jacobi matrix.
        //
        double[] bj = new double[n];
        for (i = 0; i < n; i++)
        {
            bj[i] = i + 1;
        }

        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 2 * i + 1;
        }

        double[] w = new double[n];
        w[0] = Math.Sqrt(zemu);
        for (i = 1; i < n; i++)
        {
            w[i] = 0.0;
        }

        //
        //  Diagonalize the Jacobi matrix.
        //
        IMTQLX.imtqlx(n, ref x, ref bj, ref w);


        return x;
    }

    public static double[] l_power_product(int p, int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L_POWER_PRODUCT: power product table for L(n,x).
        //
        //  Discussion:
        //
        //    Let L(n,x) represent the Laguerre polynomial of degree n.  
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of X^E with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( 0 <= X < +oo ) x^e * L(i,x) * L(j,x) * exp (-x) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, int E, the exponent of X.
        //    0 <= E.
        //
        //    Output, double L_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
        //    TABLE(I,J) represents the weighted integral of X^E * L(i,x) * L(j,x).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];

        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = p + 1 + (e + 1) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        QuadratureRule.l_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] l_table = l_polynomial(1, p, x_table, xIndex: +k);

            switch (e)
            {
                case 0:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * l_table[i] * l_table[j];
                        }
                    }

                    break;
                }
                default:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * Math.Pow(x, e) * l_table[i] * l_table[j];
                        }
                    }

                    break;
                }
            }
        }

        return table;
    }


}