using System;

namespace Burkardt.PolynomialNS;

public static class Chebyshev
{
    public static void cheb(int deg, double pt, ref double[] tcheb, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEB computes normalized Chebyshev polynomials.
        //
        //  Discussion:
        //
        //    This subroutine computes the array TCHEB of normalized Chebyshev 
        //    polynomials from degree 0 to DEG:
        //      T_0(x)=1, 
        //      T_j(x) = sqrt(2) * cos ( j * acos(x) ) 
        //    at the point x = PT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //  
        //  Modified:
        //
        //    14 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
        //    Marco Vianello.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Marco Caliari, Stefano de Marchi, Marco Vianello,
        //    Algorithm 886:
        //    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
        //    ACM Transactions on Mathematical Software,
        //    Volume 35, Number 3, October 2008, Article 21, 11 pages.
        //
        //  Parameters:
        //
        //    Input, int DEG, the degree.
        //    0 <= DEG.
        //
        //    Input, double PT, the evaluation point.
        //
        //    Output, double TCHEB[DEG+1], the value of the normalized
        //    Chebyshev polynomials of degrees 0 through DEG at the point PT.
        //
    {
        int j;
        const double sqrt2 = 1.4142135623730951;

        switch (deg)
        {
            case < 0:
                return;
        }

        tcheb[index + 0] = 1.0;

        switch (deg)
        {
            case < 1:
                return;
        }

        tcheb[index + 1] = sqrt2 * pt;

        switch (deg)
        {
            case < 2:
                return;
        }

        tcheb[index + 2] = 2.0 * pt * tcheb[index + 1] - sqrt2 * tcheb[index + 0];
        //
        //  Chebyshev recurrence.
        //
        for (j = 3; j <= deg; j++)
        {
            tcheb[index + j] = 2.0 * pt * tcheb[index + j - 1] - tcheb[index + j - 2];
        }
    }

    public static double[] cheby_t_poly(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY evaluates Chebyshev polynomials T(n,x).
        //
        //  Discussion:
        //
        //    Chebyshev polynomials are useful as a basis for representing the
        //    approximation of functions since they are well conditioned, in the sense
        //    that in the interval [-1,1] they each have maximum absolute value 1.
        //    Hence an error in the value of a coefficient of the approximation, of
        //    size epsilon, is exactly reflected in an error of size epsilon between
        //    the computed approximation and the theoretical approximation.
        //
        //    Typical usage is as follows, where we assume for the moment
        //    that the interval of approximation is [-1,1].  The value
        //    of N is chosen, the highest polynomial to be used in the
        //    approximation.  Then the function to be approximated is
        //    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
        //    Chebyshev polynomial.  Let these values be denoted by F(XJ).
        //
        //    The coefficients of the approximation are now defined by
        //
        //      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
        //
        //    except that C(0) is given a value which is half that assigned
        //    to it by the above formula,
        //
        //    and the representation is
        //
        //    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
        //
        //    Now note that, again because of the fact that the Chebyshev polynomials
        //    have maximum absolute value 1, if the higher order terms of the
        //    coefficients C are small, then we have the option of truncating
        //    the approximation by dropping these terms, and we will have an
        //    exact value for maximum perturbation to the approximation that
        //    this will cause.
        //
        //    It should be noted that typically the error in approximation
        //    is dominated by the first neglected basis function (some multiple of
        //    T(N+1,X) in the example above).  If this term were the exact error,
        //    then we would have found the minimax polynomial, the approximating
        //    polynomial of smallest maximum deviation from the original function.
        //    The minimax polynomial is hard to compute, and another important
        //    feature of the Chebyshev approximation is that it tends to behave
        //    like the minimax polynomial while being easy to compute.
        //
        //    To evaluate a sum like
        //
        //      sum ( 0 <= J <= N ) C(J) T(J,X),
        //
        //    Clenshaw's recurrence formula is recommended instead of computing the
        //    polynomial values, forming the products and summing.
        //
        //    Assuming that the coefficients C(J) have been computed
        //    for J = 0 to N, then the coefficients of the representation of the
        //    indefinite integral of the function may be computed by
        //
        //      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
        //
        //    with
        //
        //      C(N+1)=0
        //      B(0) arbitrary.
        //
        //    Also, the coefficients of the representation of the derivative of the
        //    function may be computed by:
        //
        //      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
        //
        //    with
        //
        //      D(N+1) = D(N)=0.
        //
        //    Some of the above may have to adjusted because of the irregularity of C(0).
        //
        //    The formula is:
        //
        //      T(N,X) = COS(N*ARCCOS(X))
        //
        //  Differential equation:
        //
        //    (1-X*X) Y'' - X Y' + N N Y = 0
        //
        //  First terms:
        //
        //    T(0,X) =  1
        //    T(1,X) =  1 X
        //    T(2,X) =  2 X^2 -   1
        //    T(3,X) =  4 X^3 -   3 X
        //    T(4,X) =  8 X^4 -   8 X^2 +  1
        //    T(5,X) = 16 X^5 -  20 X^3 +  5 X
        //    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
        //    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
        //
        //  Inequality:
        //
        //    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
        //
        //  Orthogonality:
        //
        //    For integration over [-1,1] with weight
        //
        //      W(X) = 1 / sqrt(1-X*X),
        //
        //    if we write the inner product of T(I,X) and T(J,X) as
        //
        //      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
        //
        //    then the result is:
        //
        //      0 if I /= J
        //      PI/2 if I == J /= 0
        //      PI if I == J == 0
        //
        //    A discrete orthogonality relation is also satisfied at each of
        //    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
        //                              = 0 if I /= J
        //                              = N/2 if I == J /= 0
        //                              = N if I == J == 0
        //
        //  Recursion:
        //
        //    T(0,X) = 1,
        //    T(1,X) = X,
        //    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
        //
        //    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
        //
        //  Special values:
        //
        //    T(N,1) = 1
        //    T(N,-1) = (-1)^N
        //    T(2N,0) = (-1)^N
        //    T(2N+1,0) = 0
        //    T(N,X) = (-1)^N * T(N,-X)
        //
        //  Zeroes:
        //
        //    M-th zero of T(N,X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
        //
        //  Extrema:
        //
        //    M-th extremum of T(N,X) is cos(PI*M/N), M = 0 to N
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
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest polynomial to compute.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double CHEBY_T_POLY[M*(N+1)], the values of the Chebyshev polynomials.
        //
    {
        int i;
        int j;
        double[] v;

        switch (n)
        {
            case < 0:
                return null;
        }

        v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i] = 1.0;
        }

        switch (n)
        {
            case < 1:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = 2.0 * x[i] * v[i + (j - 1) * m] - v[i + (j - 2) * m];
            }
        }

        return v;
    }

    public static double[] cheby_t_poly_coef(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY_COEF evaluates coefficients of Chebyshev polynomials T(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     1
        //     2     -1     0      2
        //     3      0    -3      0      4
        //     4      1     0     -8      0       8
        //     5      0     5      0    -20       0    16
        //     6     -1     0     18      0     -48     0     32
        //     7      0    -7      0     56       0  -112      0    64
        //
        //  Recursion:
        //
        //    T(0,X) = 1,
        //    T(1,X) = X,
        //    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2012
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
        //    Output, double CHEBY_T_POLY_COEF[(N+1)*(N+1)], the coefficients of 
        //    the Chebyshev T polynomials.
        //
    {
        double[] c;
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        c = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return c;
        }

        c[1 + 1 * (n + 1)] = 1.0;

        for (i = 2; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = -c[i - 2 + 0 * (n + 1)];
            for (j = 1; j <= i - 2; j++)
            {
                c[i + j * (n + 1)] = 2.0 * c[i - 1 + (j - 1) * (n + 1)] - c[i - 2 + j * (n + 1)];
            }

            c[i + (i - 1) * (n + 1)] = 2.0 * c[i - 1 + (i - 2) * (n + 1)];
            c[i + i * (n + 1)] = 2.0 * c[i - 1 + (i - 1) * (n + 1)];
        }

        return c;
    }
        
    public static double[] cheby_t_poly_zero ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_POLY_ZERO returns zeroes of Chebyshev polynomials T(n,x).
        //
        //  Discussion:
        //
        //    The I-th zero of T(N,X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Output, double CHEBY_T_POLY_ZERO[N], the zeroes of T(N,X).
        //
    {
        double angle;
        int i;
            
        double[] z;

        z = new double[n];

        for ( i = 0; i < n; i++ )
        {
            angle = (2 * i + 1) * Math.PI / (2 * n);
            z[i] = Math.Cos ( angle );
        }

        return z;
    }

    public static double[] cheby_u_poly(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY evaluates Chebyshev polynomials U(n,x).
        //
        //  Differential equation:
        //
        //    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
        //
        //  First terms:
        //
        //    U(0,X) =   1
        //    U(1,X) =   2 X
        //    U(2,X) =   4 X^2 -   1
        //    U(3,X) =   8 X^3 -   4 X
        //    U(4,X) =  16 X^4 -  12 X^2 +  1
        //    U(5,X) =  32 X^5 -  32 X^3 +  6 X
        //    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
        //    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
        //
        //  Recursion:
        //
        //    U(0,X) = 1,
        //    U(1,X) = 2 * X,
        //    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
        //
        //  Norm:
        //
        //    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest polynomial to compute.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double CHEBY_T_POLY[M*(N+1)], the values of the Chebyshev polynomials.
        //
    {
        int i;
        int j;
        double[] v;

        switch (n)
        {
            case < 0:
                return null;
        }

        v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i] = 1.0;
        }

        switch (n)
        {
            case < 1:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = 2.0 * x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = 2.0 * x[i] * v[i + (j - 1) * m] - v[i + (j - 2) * m];
            }
        }

        return v;
    }

    public static void cheby_u_poly_coef(int n, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY_COEF evaluates coefficients of Chebyshev polynomials U(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     2
        //     2     -1     0      4
        //     3      0    -4      0      8
        //     4      1     0    -12      0      16
        //     5      0     6      0    -32       0    32
        //     6     -1     0     24      0     -80     0     64
        //     7      0    -8      0     80       0  -192      0   128
        //
        //  Recursion:
        //
        //    U(0,X) = 1,
        //    U(1,X) = 2*X,
        //    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2003
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
        //    Output, double C[(N+1)*((N+1)], the coefficients of the Chebyshev U
        //    polynomials.
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

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        c[1 + 1 * (n + 1)] = 2.0;

        for (i = 2; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = -c[i - 2 + 0 * (n + 1)];
            for (j = 1; j <= i - 2; j++)
            {
                c[i + j * (n + 1)] = 2.0 * c[i - 1 + (j - 1) * (n + 1)] - c[i - 2 + j * (n + 1)];
            }

            c[i + (i - 1) * (n + 1)] = 2.0 * c[i - 1 + (i - 2) * (n + 1)];
            c[i + i * (n + 1)] = 2.0 * c[i - 1 + (i - 1) * (n + 1)];
        }
    }

    public static double[] cheby_u_poly_zero ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_POLY_ZERO returns zeroes of Chebyshev polynomials U(n,x).
        //
        //  Discussion:
        //
        //    The I-th zero of U(N,X) is cos((I-1)*PI/(N-1)), I = 1 to N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Output, double CHEBY_U_POLY_ZERO[N], the zeroes of U(N,X).
        //
    {
        double angle;
        int i;
            
        double[] z;

        z = new double[n];

        for ( i = 0; i < n; i++ )
        {
            angle = (i + 1) * Math.PI / (n + 1);
            z[i] = Math.Cos ( angle );
        }

        return z;
    }

    public static void chebyshev_discrete ( int n, int m, double x, ref double[] v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials at a point.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Walter Gautschi,
        //    Orthogonal Polynomials: Computation and Approximation,
        //    Oxford, 2004,
        //    ISBN: 0-19-850672-4,
        //    LC: QA404.5 G3555.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order of the polynomials to 
        //    be evaluated.  0 <= N <= M.
        //
        //    Input, int M, the maximum order of the polynomials.
        //    0 <= M.
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double V[N+1], the value of the polynomials at X.
        //
    {
        int i;

        switch (m)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV_DISCRETE - Fatal error!");
                Console.WriteLine("  Parameter M must be nonnegative.");
                return;
        }

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("CHEBYSHEV_DISCRETE - Fatal error!");
                Console.WriteLine("  Parameter N must be nonnegative.");
                return;
        }

        if ( m < n )
        {
            Console.WriteLine("");
            Console.WriteLine("CHEBYSHEV_DISCRETE - Fatal error!");
            Console.WriteLine("  Parameter N must be no greater than M.");
            return;
        }

        v[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        v[1] = 2.0 * x + (1 - m);

        switch (n)
        {
            case 1:
                return;
        }

        for ( i = 1; i < n; i++ )
        {
            v[i+1] = ( 
                (2 * i + 1) 
                * ( 2.0 * x + (1 - m) ) * v[i]
                - i * ( m + i ) * ( m - i ) * v[i-1]
            ) / (i + 1);
        }

    }

}