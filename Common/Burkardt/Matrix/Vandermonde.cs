using System;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class VandermondeMatrix
{
    public static double[] bivand1(int n, double[] alpha, double[] beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVAND1 returns a bidimensional Vandermonde1 matrix.
        //
        //  Discussion:
        //
        //    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
        //
        //    (x,y)   | (1,10)  (2,10)  (3,10)  (1,20)  (2,20)  (1,30)
        //    --------+-----------------------------------------------
        //    1       |     1       1       1       1       1       1  
        //    x       |     1       2       3       1       2       1
        //       y    |    10      10      10      20      20      30
        //    x^2     |     1       4       9       1       4       1
        //    x  y    |    10      20      30      20      40      30
        //    x^2y^2  |   100     100     100     400     400     900
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the data vectors.
        //
        //    Input, double ALPHA[N], BETA[N], the values that define A.
        //
        //    Output, double BIVAND1[((N+1)*N)/2*((N+1)*N)/2], the matrix.
        //
    {
        double[] a;
        int e;
        int e1;
        int e2;
        int ii;
        int j1;
        int j2;
        int jj;
        int n2;

        n2 = n * (n + 1) / 2;
        a = new double[n2 * n2];

        e1 = 0;
        e2 = 0;
        e = 0;

        for (ii = 0; ii < n2; ii++)
        {
            j1 = 0;
            j2 = 0;
            for (jj = 0; jj < n2; jj++)
            {
                a[ii + jj * n2] = ii switch
                {
                    0 => 1.0,
                    _ => Math.Pow(alpha[j1], e1) * Math.Pow(beta[j2], e2)
                };

                if (j1 + j2 < n - 1)
                {
                    j1 += 1;
                }
                else
                {
                    j1 = 0;
                    j2 += 1;
                }
            }

            if (e2 < e)
            {
                e1 -= 1;
                e2 += 1;
            }
            else
            {
                e += 1;
                e1 = e;
                e2 = 0;
            }
        }

        return a;
    }

    public static double[] bivand2(int n, double[] alpha, double[] beta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVAND2 returns a bidimensional Vandermonde1 matrix.
        //
        //  Discussion:
        //
        //    N = 3, ALPHA = ( 1, 2, 3 ), BETA = ( 10, 20, 30 )
        //
        //    (x,y)   | (1,10) (2,10) (3,10) (1,20) (2,20) (3,20) (1,30) (2,30) (3,30)
        //    --------+---------------------------------------------------------------
        //    1       |     1      1      1      1      1      1      1      1      1  
        //    x       |     1      2      3      1      2      1      1      2      3
        //    x^2     |     1      4      9      1      4      1      1      4      9
        //       y    |    10     10     10     20     20     20     30     30     30
        //    x  y    |    10     20     30     20     40     60     30     60     90
        //    x^2y    |    10     40     90     20     80    180     30    120    270
        //       y^2  |   100    100    100    400    400    400    900    900    900
        //    x  y^2  |   100    200    300    400    800   1200    900   1800   2700
        //    x^2y^2  |   100    400    900    400   1600   3600    900   3600   8100
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the data vectors.
        //
        //    Input, double ALPHA[N], BETA[N], the values that define A.
        //
        //    Output, double BIVAND1[(N*N)*(N*N)], the matrix.
        //
    {
        double[] a;
        int i;
        int ix;
        int iy;
        int j;
        int jx;
        int jy;

        a = new double[n * n * n * n];

        i = 0;
        for (iy = 0; iy < n; iy++)
        {
            for (ix = 0; ix < n; ix++)
            {
                j = 0;
                for (jy = 0; jy < n; jy++)
                {
                    for (jx = 0; jx < n; jx++)
                    {
                        a[i + j * n * n] = Math.Pow(alpha[jx], ix) * Math.Pow(beta[jy], iy);
                        j += 1;
                    }
                }

                i += 1;
            }
        }

        return a;
    }

    public static double[] dvand(int n, double[] alpha, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVAND solves a Vandermonde system A' * x = b.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ake Bjorck, Victor Pereyra,
        //    Solution of Vandermonde Systems of Equations,
        //    Mathematics of Computation,
        //    Volume 24, Number 112, October 1970, pages 893-903.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double ALPHA[N], the parameters that define the matrix.
        //    The values should be distinct.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double DVAND[N], the solution of the linear system.
        //
    {
        int j;
        int k;
        double[] x;

        x = typeMethods.r8vec_copy_new(n, b);

        for (k = 0; k < n - 1; k++)
        {
            for (j = n - 1; k < j; j--)
            {
                x[j] = (x[j] - x[j - 1]) / (alpha[j] - alpha[j - k - 1]);
            }
        }

        for (k = n - 2; 0 <= k; k--)
        {
            for (j = k; j < n - 1; j++)
            {
                x[j] -= alpha[k] * x[j + 1];
            }
        }

        return x;
    }

    public static void dvandprg(int n, double[] alpha, double[] b, ref double[] x, ref double[] c,
            ref double[] m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DVANDPRG solves a Vandermonde system A' * x = f progressively.
        //
        //  Discussion:
        //
        //    This function receives the solution to the system of equations A' * x = f
        //    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
        //    and new values alpha(n) and f(n).  It updates the solution.
        //
        //    To solve a system of Nbig equations, this function may be called 
        //    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
        //    current subsystem is returned.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ake Bjorck, Victor Pereyra,
        //    Solution of Vandermonde Systems of Equations,
        //    Mathematics of Computation,
        //    Volume 24, Number 112, October 1970, pages 893-903.
        //
        //  Parameters:
        //
        //    Input, int N, the new order of the matrix, which is 1 
        //    larger than on the previous call.  For the first call, N must be 1.
        //
        //    Input, double ALPHA[N], the parameters that define the matrix.
        //    The values should be distinct.  The value ALPHA(N) has just been
        //    added to the system.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N].  On input, the first N-1 entries 
        //    contain the solution of the N-1xN-1 linear system.  On output, the 
        //    solution to the NxN linear system.
        //
        //    Input/output, double C[N], M[N].  On input, the first N-1 
        //    entries contain factorization data for the N-1xN-1 linear system.  On 
        //    output, factorization data for the NxN linear system.
        //
    {
        double cn;
        int j;

        c[n - 1] = b[n - 1];
        for (j = n - 1; 1 <= j; j--)
        {
            c[j - 1] = (c[j] - c[j - 1]) / (alpha[n - 1] - alpha[j - 1]);
        }

        m[n - 1] = n switch
        {
            1 => 1.0,
            _ => 0.0
        };

        cn = c[0];
        x[n - 1] = c[0];

        for (j = n - 1; 1 <= j; j--)
        {
            m[j] -= alpha[n - 2] * m[j - 1];
            x[n - j - 1] += m[j] * cn;
        }
    }

    public static double[] pvand(int n, double[] alpha, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PVAND solves a Vandermonde system A * x = b.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ake Bjorck, Victor Pereyra,
        //    Solution of Vandermonde Systems of Equations,
        //    Mathematics of Computation,
        //    Volume 24, Number 112, October 1970, pages 893-903.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double ALPHA[N], the parameters that define the matrix.
        //    The values should be distinct.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double PVAND[N], the solution of the linear system.
        //
    {
        int j;
        int k;
        double[] x;

        x = typeMethods.r8vec_copy_new(n, b);

        for (k = 0; k < n - 1; k++)
        {
            for (j = n - 1; k < j; j--)
            {
                x[j] -= alpha[k] * x[j - 1];
            }
        }

        for (k = n - 2; 0 <= k; k--)
        {
            for (j = k + 1; j < n; j++)
            {
                x[j] /= (alpha[j] - alpha[j - k - 1]);
            }

            for (j = k; j < n - 1; j++)
            {
                x[j] -= x[j + 1];
            }
        }

        return x;
    }

    public static void pvandprg(int n, double[] alpha, double[] b, ref double[] x, ref double[] d,
            ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PVANDPRG solves a Vandermonde system A * x = f progressively.
        //
        //  Discussion:
        //
        //    This function receives the solution to the system of equations A * x = f
        //    where A is a Vandermonde matrix for alpha(0) through alpha(n-1),
        //    and new values alpha(n) and f(n).  It updates the solution.
        //
        //    To solve a system of Nbig equations, this function may be called 
        //    repeatedly, with N = 1, 2, ..., Nbig.  Each time, a solution to the 
        //    current subsystem is returned.
        //
        //    Note that the reference, which lists an Algol version of this algorithm, 
        //    omits a minus sign, writing
        //      u[j] := u[j] x delta;
        //    where
        //      u[j] := - u[j] x delta;
        //    is actually necessary.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Ake Bjorck, Victor Pereyra,
        //    Solution of Vandermonde Systems of Equations,
        //    Mathematics of Computation,
        //    Volume 24, Number 112, October 1970, pages 893-903.
        //
        //  Parameters:
        //
        //    Input, int N, the new order of the matrix, which is 1 
        //    larger than on the previous call.  For the first call, N must be 1.
        //
        //    Input, double ALPHA[N], the parameters that define the matrix.
        //    The values should be distinct.  The value ALPHA(N) has just been
        //    added to the system.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N]; on input, the solution of the 
        //    N-1xN-1 linear system.  On output, the solution of the NxN linear system.
        //
        //    Input/output, double D[N], U[N]; on input, factorization data 
        //    for the N-1xN-1 linear system.  On output, factorization data for the
        //    NxN linear system.
        //
    {
        double delta;
        double dn;
        int j;

        d[n - 1] = b[n - 1];
        for (j = n - 1; 1 <= j; j--)
        {
            d[j - 1] = d[j] - alpha[n - j - 1] * d[j - 1];
        }

        dn = d[0];
        u[n - 1] = 1.0;

        for (j = 1; j <= n - 1; j++)
        {
            delta = alpha[n - 1] - alpha[j - 1];
            u[j - 1] = -u[j - 1] * delta;
            u[n - 1] *= delta;
            x[j - 1] += dn / u[j - 1];
        }

        x[n - 1] = dn / u[n - 1];
    }

    public static double[] cheby_van1(int m, double a, double b, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_VAN1 returns the CHEBY_VAN1 matrix.
        //
        //  Discussion:
        //
        //    Normally, the Chebyshev polynomials are defined on -1 <= XI <= +1.
        //    Here, we assume the Chebyshev polynomials have been defined on the
        //    interval A <= X <= B, using the mapping
        //      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
        //    so that
        //      ChebyAB(A,B;X) = Cheby(XI).
        //
        //    if ( I == 1 ) then
        //      V(1,1:N) = 1;
        //    elseif ( I == 2 ) then
        //      V(2,1:N) = XI(1:N);
        //    else
        //      V(I,1:N) = 2.0 * XI(1:N) * V(I-1,1:N) - V(I-2,1:N);
        //
        //  Example:
        //
        //    N = 5, A = -1, B = +1, X = ( 1, 2, 3, 4, 5 )
        //
        //    1  1   1    1    1
        //    1  2   3    4    5
        //    1  7  17   31   49
        //    1 26  99  244  485
        //    1 97 577 1921 4801
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Nicholas Higham,
        //    Stability analysis of algorithms for solving confluent
        //    Vandermonde-like systems,
        //    SIAM Journal on Matrix Analysis and Applications,
        //    Volume 11, 1990, pages 23-41.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, double A, B, the interval.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, double X[N], the vector that defines the matrix.
        //
        //    Output, double CHEBY_VAN1[M*N], the matrix.
        //
    {
        int i;
        int j;
        double[] v;
        double xi;

        v = new double[m * n];

        for (j = 0; j < n; j++)
        {
            xi = (-(b - x[j]) + (x[j] - a)) / (b - a);
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = i switch
                {
                    0 => 1.0,
                    1 => xi,
                    _ => 2.0 * xi * v[i - 1 + j * m] - v[i - 2 + j * m]
                };
            }
        }

        return v;
    }

    public static double[] legendre_van(int m, double a, double b, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_VAN returns the LEGENDRE_VAN matrix.
        //
        //  Discussion:
        //
        //    The LEGENDRE_VAN matrix is the Legendre Vandermonde-like matrix.
        //
        //    Normally, the Legendre polynomials are defined on -1 <= XI <= +1.
        //    Here, we assume the Legendre polynomials have been defined on the
        //    interval A <= X <= B, using the mapping
        //      XI = ( - ( B - X ) + ( X - A ) ) / ( B - A )
        //    so that
        //      Lab(A,B;X) = L(XI).
        //
        //    if ( I = 1 ) then
        //      V(1,1:N) = 1
        //    else if ( I = 2 ) then
        //      V(2,1:N) = XI(1:N)
        //    else
        //      V(I,1:N) = ( (2*I-1) * XI(1:N) * V(I-1,1:N) - (I-1)*V(I-2,1:N) ) / I
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, double A, B, the interval.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, double X[N], the vector that defines the matrix.
        //
        //    Output, double LEGENDRE_VAN[M*N], the matrix.
        //
    {
        int i;
        int j;
        double[] v;
        double xi;

        v = new double[m * n];

        for (j = 0; j < n; j++)
        {
            xi = (-(b - x[j]) + (x[j] - a)) / (b - a);
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = i switch
                {
                    0 => 1.0,
                    1 => xi,
                    _ => ((2 * i - 1) * xi * v[i - 1 + j * m] + (-i + 1) * v[i - 2 + j * m]) / i
                };
            }
        }

        return v;
    }

    public static double[] vand1(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VAND1 returns the Vandermonde1 matrix A with 1's on the first row.
        //
        //  Formula:
        //
        //    A(I,J) = X(J)^(I-1)
        //
        //  Example:
        //
        //    N = 5, X = ( 2, 3, 4, 5, 6 )
        //
        //    1  1   1   1   1
        //    2  3   4   5   6
        //    4  9  16  25  36
        //    8 27  64 125  216
        //   16 81 256 625 1296
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //    A is nonsingular if, and only if, the X values are distinct.
        //
        //    det ( A ) = product ( 1 <= I <= N ) ( 1 <= J .lt. I ) ( X(I) - X(J) ).
        //             = product ( 1 <= J <= N ) X(J)
        //             * product ( 1 <= I .lt. J ) ( X(J) - X(I) ).
        //
        //    A is generally ill-conditioned.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969, page 27,
        //    LC: QA263.G68.
        //
        //    Nicholas Higham,
        //    Stability analysis of algorithms for solving confluent
        //    Vandermonde-like systems,
        //    SIAM Journal on Matrix Analysis and Applications,
        //    Volume 11, 1990, pages 23-41.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix desired.
        //
        //    Input, double X[N], the values that define A.
        //
        //    Output, double VAND1[N*N], the matrix.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i + j * n] = i switch
                {
                    0 when x[j] == 0.0 => 1.0,
                    _ => Math.Pow(x[j], i)
                };
            }
        }

        return a;
    }

    public static double[] vandermonde_approx_1d_matrix(int n, int m, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_APPROX_1D_MATRIX computes a Vandermonde 1D approximation matrix.
        //
        //  Discussion:
        //
        //    We assume the approximant has the form
        //
        //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m.
        //
        //    We have n data values (x(i),y(i)) which must be approximated:
        //
        //      p(x(i)) = c0 + c1 * x(i) + c2 * x(i)^2 + ... + cm * x(i)^m = y(i)
        //
        //    This can be cast as an Nx(M+1) linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 x1^2 ... x1^m ] [  c0 ] = [  y1 ]
        //      [ 1 x2 x2^2 ... x2^m ] [  c1 ] = [  y2 ]
        //      [ .................. ] [ ... ] = [ ... ]
        //      [ 1 xn xn^2 ... xn^m ] [  cm ] = [  yn ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double X(N), the data values.
        //
        //    Output, double VANDERMONDE_APPROX_1D_MATRIX[N*(M+1)], the Vandermonde matrix for X.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[n * (m + 1)];

        for (i = 0; i < n; i++)
        {
            a[i + 0 * n] = 1.0;
        }

        for (j = 1; j <= m; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = a[i + (j - 1) * n] * x[i];
            }
        }

        return a;
    }

    public static double[] vandermonde_approx_2d_matrix(int n, int m, int tm, double[] x,
            double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VANDERMONDE_APPROX_2D_MATRIX computes a Vandermonde 2D approximation matrix.
        //
        //  Discussion:
        //
        //    We assume the approximating function has the form of a polynomial
        //    in X and Y of total degree M.
        //
        //      p(x,y) = c00 
        //             + c10 * x                + c01 * y
        //             + c20 * x^2   + c11 * xy + c02 * y^2
        //             + ...
        //             + cm0 * x^(m) + ...      + c0m * y^m.
        //
        //    If we let T(K) = the K-th triangular number 
        //            = sum ( 1 <= I <= K ) I
        //    then the number of coefficients in the above polynomial is T(M+1).
        //
        //    We have n data locations (x(i),y(i)) and values z(i) to approximate:
        //
        //      p(x(i),y(i)) = z(i)
        //
        //    This can be cast as an NxT(M+1) linear system for the polynomial
        //    coefficients:
        //
        //      [ 1 x1 y1  x1^2 ... y1^m ] [ c00 ] = [  z1 ]
        //      [ 1 x2 y2  x2^2 ... y2^m ] [ c10 ] = [  z2 ]
        //      [ 1 x3 y3  x3^2 ... y3^m ] [ c01 ] = [  z3 ]
        //      [ ...................... ] [ ... ] = [ ... ]
        //      [ 1 xn yn  xn^2 ... yn^m ] [ c0m ] = [  zn ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, int TM, the M+1st triangular number.
        //
        //    Input, double X[N], Y[N], the data locations.
        //
        //    Output, double VANDERMONDE_APPROX_2D_MATRIX[N*TM], the Vandermonde matrix for X.
        //
    {
        double[] a;
        int ex;
        int ey;
        int i;
        int j;
        int s;

        a = new double[n * tm];
        j = 0;

        for (s = 0; s <= m; s++)
        {
            for (ex = s; 0 <= ex; ex--)
            {
                ey = s - ex;
                for (i = 0; i < n; i++)
                {
                    a[i + j * n] = Math.Pow(x[i], ex) * Math.Pow(y[i], ey);
                }

                j += 1;
            }
        }

        return a;
    }
}