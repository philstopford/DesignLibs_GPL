using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Vandermonde
    {
        public static double[] bivand1(int n, double[] alpha, double[] beta )

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

            n2 = (n * (n + 1)) / 2;
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
                    if (ii == 0)
                    {
                        a[ii + jj * n2] = 1.0;
                    }
                    else
                    {
                        a[ii + jj * n2] = Math.Pow(alpha[j1], e1) * Math.Pow(beta[j2], e2);
                    }

                    if (j1 + j2 < n - 1)
                    {
                        j1 = j1 + 1;
                    }
                    else
                    {
                        j1 = 0;
                        j2 = j2 + 1;
                    }
                }

                if (e2 < e)
                {
                    e1 = e1 - 1;
                    e2 = e2 + 1;
                }
                else
                {
                    e = e + 1;
                    e1 = e;
                    e2 = 0;
                }
            }

            return a;
        }

        public static double[] bivand2(int n, double[] alpha, double[] beta )

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
                            j = j + 1;
                        }
                    }

                    i = i + 1;
                }
            }

            return a;
        }

        public static double[] dvand(int n, double[] alpha, double[] b )

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
                    x[j] = x[j] - alpha[k] * x[j + 1];
                }
            }

            return x;
        }

        public static void dvandprg(int n, double[] alpha, double[] b, ref double[] x, ref double[] c,
        ref double[] m )

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

            if (n == 1)
            {
                m[n - 1] = 1.0;
            }
            else
            {
                m[n - 1] = 0.0;
            }

            cn = c[0];
            x[n - 1] = c[0];

            for (j = n - 1; 1 <= j; j--)
            {
                m[j] = m[j] - alpha[n - 2] * m[j - 1];
                x[n - j - 1] = x[n - j - 1] + m[j] * cn;
            }
        }

        public static double[] pvand(int n, double[] alpha, double[] b )

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
                    x[j] = x[j] - alpha[k] * x[j - 1];
                }
            }

            for (k = n - 2; 0 <= k; k--)
            {
                for (j = k + 1; j < n; j++)
                {
                    x[j] = x[j] / (alpha[j] - alpha[j - k - 1]);
                }

                for (j = k; j < n - 1; j++)
                {
                    x[j] = x[j] - x[j + 1];
                }
            }

            return x;
        }

        public static void pvandprg(int n, double[] alpha, double[] b, ref double[] x, ref double[] d,
        ref double[] u )

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
                u[n - 1] = u[n - 1] * delta;
                x[j - 1] = x[j - 1] + dn / u[j - 1];
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
                    if (i == 0)
                    {
                        v[i + j * m] = 1.0;
                    }
                    else if (i == 1)
                    {
                        v[i + j * m] = xi;
                    }
                    else
                    {
                        v[i + j * m] = 2.0 * xi * v[i - 1 + j * m] - v[i - 2 + j * m];
                    }
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
                    if (i == 0)
                    {
                        v[i + j * m] = 1.0;
                    }
                    else if (i == 1)
                    {
                        v[i + j * m] = xi;
                    }
                    else
                    {
                        v[i + j * m] = ((double) (2 * i - 1) * xi * v[i - 1 + j * m] +
                                        (double) (-i + 1) * v[i - 2 + j * m])
                                       / (double) (i);
                    }
                }
            }

            return v;
        }
        
        public static double[] vand1 ( int n, double[] x )

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

            a = new double[n*n];

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    if ( i == 0 && x[j] == 0.0 )
                    {
                        a[i+j*n] = 1.0;
                    }
                    else
                    {
                        a[i+j*n] = Math.Pow ( x[j], i );
                    }
                }
            }

            return a;
        }
    }
}