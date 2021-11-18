using System;
using System.Globalization;
using Burkardt.SubsetNS;

namespace Burkardt.Types;

public static partial class typeMethods
{

    public static double[] r8mat_vand2 ( int n, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_VAND2 returns the N by N row Vandermonde matrix A.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The row Vandermonde matrix returned by this routine reads "across"
        //    rather than down.  In particular, each row begins with a 1, followed by
        //    some value X, followed by successive powers of X.
        //
        //    The formula for the matrix entries is:
        //
        //      A(I,J) = X(I)^(J-1)
        //
        //  Properties:
        //
        //    A is nonsingular if, and only if, the X values are distinct.
        //
        //    The determinant of A is
        //
        //      det(A) = product ( 2 <= I <= N ) (
        //        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).
        //
        //    The matrix A is generally ill-conditioned.
        //
        //  Example:
        //
        //    N = 5, X = (2, 3, 4, 5, 6)
        //
        //    1 2  4   8   16
        //    1 3  9  27   81
        //    1 4 16  64  256
        //    1 5 25 125  625
        //    1 6 36 216 1296
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix desired.
        //
        //    Input, double X[N], the values that define A.
        //
        //    Output, double R8MAT_VAND2[N*N], the N by N row Vandermonde matrix.
        //
    {
        int i;

        double[] a = new double[n*n];

        for ( i = 0; i < n; i++ )
        {
            int j;
            for ( j = 0; j < n; j++ )
            {
                a[i + j * n] = j switch
                {
                    0 when x[i] == 0.0 => 1.0,
                    _ => Math.Pow(x[i], j)
                };
            }
        }

        return a;
    }
        
    public static void r8mat_shortest_path ( int n, ref double[] m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SHORTEST_PATH computes the shortest distance between all pairs of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Floyd,
        //    Algorithm 97, Shortest Path,
        //    Communications of the ACM,
        //    Volume 5, Number 6, June 1962, page 345.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, double M[N*N].
        //    On input, M(I,J) contains the length of the direct link between 
        //    nodes I and J, or HUGE if there is no direct link.
        //    On output, M(I,J) contains the distance between nodes I and J,
        //    that is, the length of the shortest path between them.  If there
        //    is no such path, then M(I,J) will remain HUGE.
        //
    {
        int i;
        const double r8_inf = 1.0E+30;

        for ( i = 0; i < n; i++ )
        {
            int j;
            for ( j = 0; j < n; j++ )
            {
                switch (m[j+i*n])
                {
                    case < r8_inf:
                    {
                        int k;
                        for ( k = 0; k < n; k++ )
                        {
                            switch (m[i+k*n])
                            {
                                case < r8_inf:
                                {
                                    double s = m[j+i*n] + m[i+k*n];
                                    if ( s < m[j+k*n] )
                                    {
                                        m[j+k*n] = s;
                                    }

                                    break;
                                }
                            }
                        }

                        break;
                    }
                }
            }
        }
    }
    public static double[] r8mat_sub_new ( int m, int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SUB_NEW computes C = A - B.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrices.
        //
        //    Input, double A[M*N], double B[M*N], the matrices.
        //
        //    Output, double R8MAT_SUB_NEW[M*N], the value of A-B.
        //
    {
        int j;

        double[] c = new double[m*n];

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < n; i++ )
            {
                c[i+j*m] = a[i+j*m] - b[i+j*m];
            }
        }

        return c;
    }
    public static void r8mat_mm ( int n1, int n2, int n3, double[] a, double[] b, ref double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MM multiplies two matrices.
        //
        //  Discussion: 							    
        //
        //    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double C[N1*N3], the product matrix C = A * B.
        //
    {
        int i;

        for ( i = 0; i < n1; i ++ )
        {
            int j;
            for ( j = 0; j < n3; j++ )
            {
                c[i+j*n1] = 0.0;
                int k;
                for ( k = 0; k < n2; k++ )
                {
                    c[i+j*n1] += a[i+k*n1] * b[k+j*n2];
                }
            }
        }
    }

    public static double[] r8mat_mm_new(int n1, int n2, int n3, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MM_NEW multiplies two matrices.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
        //
    {
        int i;

        double[] c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[i + k * n1] * b[k + j * n2];
                }
            }
        }

        return c;
    }
        
    public static double r8mat_diff_frobenius(int m, int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIFF_FROBENIUS returns the Frobenius norm of the difference of R8MAT's.
        //
        //  Discussion: 							    
        //
        //    An R8MAT is a doubly dimensioned array of double precision values, which
        //    may be stored as a vector in column-major order.
        //
        //    The Frobenius norm is defined as
        //
        //      R8MAT_NORM_FRO = sqrt (
        //        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
        //
        //    The matrix Frobenius norm is not derived from a vector norm, but
        //    is compatible with the vector L2 norm, so that:
        //
        //      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], double B[M*N], the matrices for which we
        //    want the Frobenius norm of the difference.
        //
        //    Output, double R8MAT_DIFF_FROBENIUS, the Frobenius norm of ( A - B ).
        //
    {
        int j;

        double value = 0.0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                value += Math.Pow(a[i + j * m] - b[i + j * m], 2);
            }
        }

        value = Math.Sqrt(value);

        return value;
    }

    public static double[] r8mat_hess(Func<int, double[], double> fx, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_HESS approximates a Hessian matrix via finite differences.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    H(I,J) = d2 F / d X(I) d X(J)
        //
        //    The values returned by this routine will be only approximate.
        //    In some cases, they will be so poor that they are useless.
        //    However, one of the best applications of this routine is for
        //    checking your own Hessian calculations, since as Heraclitus
        //    said, you'll never get the same result twice when you differentiate
        //    a complicated expression by hand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double *FX ( int N, double X[] ), the name of the user
        //    function routine.
        //
        //    Input, int N, the number of variables.
        //
        //    Input, double X[N], the values of the variables.
        //
        //    Output, double H[N*N], the approximated N by N Hessian matrix.
        //
    {
        double fmm;
        double fpp;
        int i;
        double xi;
        //
        //  Choose the stepsizes.
        //
        double[] s = new double[n];

        double eps = Math.Pow(r8_epsilon(), 0.33);

        for (i = 0; i < n; i++)
        {
            s[i] = eps * Math.Max(Math.Abs(x[i]), 1.0);
        }

        //
        //  Calculate the diagonal elements.
        //
        double[] h = new double[n * n];

        for (i = 0; i < n; i++)
        {
            xi = x[i];

            double f00 = fx(n, x);

            x[i] = xi + s[i];
            fpp = fx(n, x);

            x[i] = xi - s[i];
            fmm = fx(n, x);

            h[i + i * n] = (fpp - f00 + (fmm - f00)) / s[i] / s[i];

            x[i] = xi;
        }

        //
        //  Calculate the off diagonal elements.
        //
        for (i = 0; i < n; i++)
        {
            xi = x[i];

            int j;
            for (j = i + 1; j < n; j++)
            {
                double xj = x[j];

                x[i] = xi + s[i];
                x[j] = xj + s[j];
                fpp = fx(n, x);

                x[i] = xi + s[i];
                x[j] = xj - s[j];
                double fpm = fx(n, x);

                x[i] = xi - s[i];
                x[j] = xj + s[j];
                double fmp = fx(n, x);

                x[i] = xi - s[i];
                x[j] = xj - s[j];
                fmm = fx(n, x);

                h[j + i * n] = (fpp - fpm + (fmm - fmp)) / (4.0 * s[i] * s[j]);

                h[i + j * n] = h[j + i * n];

                x[j] = xj;
            }

            x[i] = xi;
        }

        return x;
    }
        
    public static bool r8mat_in_01 ( int m, int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Output, bool R8MAT_IN_01, is TRUE if every entry of A is
        //    between 0 and 1.
        //
    {
        int j;

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                switch (a[i+j*m])
                {
                    case < 0.0:
                    case > 1.0:
                        return false;
                }
            }
        }

        return true;
    }

    public static double[] r8mat_indicator_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_INDICATOR_NEW sets up an "indicator" R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The value of each entry suggests its location, as in:
        //
        //      11  12  13  14
        //      21  22  23  24
        //      31  32  33  34
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8MAT_INDICATOR_NEW[M*N], the table.
        //
    {
        int i;

        double[] a = new double[m * n];

        int fac = (int) Math.Pow(10, Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return a;
    }


    public static double[] r8mat_jac(int m, int n, double eps,
            Func<int, int, double[], double[]> fx, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_JAC estimates a dense jacobian matrix of the function FX.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    FPRIME(I,J) = d F(I) / d X(J).
        //
        //    The jacobian is assumed to be dense, and the LINPACK/LAPACK
        //    double precision general matrix storage mode ("DGE") is used.
        //
        //    Forward differences are used, requiring N+1 function evaluations.
        //
        //    Values of EPS have typically been chosen between
        //    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
        //    machine tolerance.
        //
        //    If EPS is too small, then F(X+EPS) will be the same as
        //    F(X), and the jacobian will be full of zero entries.
        //
        //    If EPS is too large, the finite difference estimate will
        //    be inaccurate.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of functions.
        //
        //    Input, int N, the number of variables.
        //
        //    Input, double EPS, a tolerance to be used for shifting the
        //    X values during the finite differencing.  No single value
        //    of EPS will be reliable for all vectors X and functions FX.
        //
        //    Input, double *(*FX) ( int m, int n, double x[] ), the name of
        //    the user written routine which evaluates the M-dimensional
        //    function at a given N-dimensional point X.
        //
        //    Input, double X[N], the point where the jacobian
        //    is to be estimated.
        //
        //    Output, double R8MAT_JAC[M*N], the estimated jacobian matrix.
        //
    {
        int j;

        double[] fprime = new double[m * n];
        //
        //  Evaluate the function at the base point, X.
        //
        double[] work2 = fx(m, n, x);
        //
        //  Now, one by one, vary each component J of the base point X, and
        //  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
        //
        for (j = 0; j < n; j++)
        {
            double xsave = x[j];
            double del = eps * (1.0 + Math.Abs(x[j]));
            x[j] += del;
            double[] work1 = fx(m, n, x);
            x[j] = xsave;
            int i;
            for (i = 0; i < m; i++)
            {
                fprime[i + j * m] = (work1[i] - work2[i]) / del;
            }

        }

        return fprime;
    }

    public static double[] r8mat_kronecker(int m1, int n1, double[] a, int m2, int n2,
            double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_KRONECKER computes the Kronecker product of two R8MAT's.
        //
        //  Discussion:
        //
        //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
        //
        //    If A is an M1 by N1 array, and B is an M2 by N2 array, then
        //    the Kronecker product of A and B is an M1*M2 by N1*N2 array
        //      C(I,J) = A(I1,J1) * B(I2,J2)
        //    where
        //      I1 =       I   / M2
        //      I2 = mod ( I,    M2 )
        //      J1 =       J   / N2
        //      J2 = mod ( J,    N2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M1, N1, the order of the first matrix.
        //
        //    Input, double A[M1*N1], the first matrix.
        //
        //    Input, int M2, N2, the order of the second matrix.
        //
        //    Input, double B[M2*N2], the second matrix.
        //
        //    Output, double R8MAT_KRONECKER[(M1*M2)*(N1*N2)], the Kronecker product.
        //
    {
        int j1;

        int m = m1 * m2;
        int n = n1 * n2;
        double[] c = new double[m * n];

        for (j1 = 0; j1 < n1; j1++)
        {
            int i1;
            for (i1 = 0; i1 < m1; i1++)
            {
                int i0 = i1 * m2;
                int j0 = j1 * n2;
                int j = j0;
                int j2;
                for (j2 = 0; j2 < n2; j2++)
                {
                    int i = i0;
                    int i2;
                    for (i2 = 0; i2 < m2; i2++)
                    {
                        c[i + j * m] = a[i1 + j1 * m1] * b[i2 + j2 * m2];
                        i += 1;
                    }

                    j += 1;
                }
            }
        }

        return c;
    }


    public static void r8mat_lu(int m, int n, double[] a, ref double[] l, ref double[] p, ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_LU computes the LU factorization of a rectangular R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The routine is given an M by N matrix A, and produces
        //
        //      L, an M by M unit lower triangular matrix,
        //      U, an M by N upper triangular matrix, and
        //      P, an M by M permutation matrix P,
        //
        //    so that
        //
        //      A = P' * L * U.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the M by N matrix to be factored.
        //
        //    Output, double L[M*M], the M by M unit lower triangular factor.
        //
        //    Output, double P[M*M], the M by M permutation matrix.
        //
        //    Output, double U[M*N], the M by N upper triangular factor.
        //
    {
        int i;
        int j;
        //
        //  Initialize:
        //
        //    U:=A
        //    L:=Identity
        //    P:=Identity
        //
        r8mat_copy(m, n, a, ref u);

        r8mat_zeros(m, m, ref l);
        r8mat_zeros(m, m, ref p);
        for (i = 0; i < m; i++)
        {
            l[i + i * m] = 1.0;
            p[i + i * m] = 1.0;
        }

        //
        //  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
        //
        for (j = 0; j < Math.Min(m - 1, n); j++)
        {
            double pivot = 0.0;
            int ipiv = -1;

            for (i = j; i < m; i++)
            {
                if (!(pivot < Math.Abs(u[i + j * m])))
                {
                    continue;
                }

                pivot = Math.Abs(u[i + j * m]);
                ipiv = i;
            }

            //
            //  Unless IPIV is zero, swap rows J and IPIV.
            //
            if (ipiv == -1)
            {
                continue;
            }

            int k;
            for (k = 0; k < n; k++)
            {
                double t = u[j + k * m];
                u[j + k * m] = u[ipiv + j * m];
                u[ipiv + k * m] = t;

                t = l[j + k * m];
                l[j + k * m] = l[ipiv + j * m];
                l[ipiv + k * m] = t;

                t = p[j + k * m];
                p[j + k * m] = p[ipiv + j * m];
                p[ipiv + k * m] = t;
            }

            //
            //  Zero out the entries in column J, from row J+1 to M.
            //
            for (i = j + 1; i < m; i++)
            {
                if (u[i + j * m] == 0.0)
                {
                    continue;
                }

                l[i + j * m] = u[i + j * m] / u[j + j * m];

                u[i + j * m] = 0.0;

                for (k = j + 1; k < n; k++)
                {
                    u[i + k * m] -= l[i + j * m] * u[j + k * m];
                }
            }
        }
    }


    public static double[] r8mat_mmt_new(int n1, int n2, int n3, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MMT_NEW computes C = A * B'.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
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
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.
        //
        //    Output, double R8MAT_MMT_NEW[N1*N3], the product matrix C = A * B'.
        //
    {
        int i;

        double[] c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[i + k * n1] * b[j + k * n3];
                }
            }
        }

        return c;
    }

    public static double[] r8mat_mtm_new(int n1, int n2, int n3, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_MTM_NEW computes C = A' * B.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8MAT_MTM_NEW[N1*N3], the product matrix C = A' * B.
        //
    {
        int i;

        double[] c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[k + i * n2] * b[k + j * n2];
                }
            }
        }

        return c;
    }
        
    public static void r8mat_nint(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_NINT rounds the entries of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the matrix.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                int s = a[i + j * m] switch
                {
                    < 0.0 => -1,
                    _ => 1
                };

                a[i + j * m] = s * (int) (Math.Abs(a[i + j * m]) + 0.5);
            }
        }
    }

    public static int r8mat_nonzeros(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_NONZEROS returns the number of nonzeros in an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Output, int R8MAT_NONZEROS, the number of nonzeros.
        //
    {
        int j;

        int value = 0;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (a[i + j * m] != 0.0)
                {
                    value += 1;
                }
            }
        }

        return value;
    }

    public static double[] r8mat_nullspace(int m, int n, double[] a, int nullspace_size)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_NULLSPACE computes the nullspace of a matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    Let A be an MxN matrix.
        //
        //    If X is an N-vector, and A*X = 0, then X is a null vector of A.
        //
        //    The set of all null vectors of A is called the nullspace of A.
        //
        //    The 0 vector is always in the null space.
        //
        //    If the 0 vector is the only vector in the nullspace of A, then A
        //    is said to have maximum column rank.  (Because A*X=0 can be regarded
        //    as a linear combination of the columns of A).  In particular, if A
        //    is square, and has maximum column rank, it is nonsingular.
        //
        //    The dimension of the nullspace is the number of linearly independent
        //    vectors that span the nullspace.  If A has maximum column rank,
        //    its nullspace has dimension 0.
        //
        //    This routine uses the reduced row echelon form of A to determine
        //    a set of NULLSPACE_SIZE independent null vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input, double A[M*N], the matrix to be analyzed.
        //
        //    Input, int NULLSPACE_SIZE, the size of the nullspace.
        //
        //    Output, double R8MAT_NULLSPACE[N*NULLSPACE_SIZE], vectors that
        //    span the nullspace.
        //
    {
        int i;
        int j;
        //
        //  Make a copy of A.
        //
        double[] rref = r8mat_copy_new(m, n, a);
        //
        //  Get the reduced row echelon form of A.
        //
        r8mat_rref(m, n, ref rref);
        //
        //  Note in ROW the columns of the leading nonzeros.
        //  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
        //
        int[] row = new int[m];
        for (i = 0; i < m; i++)
        {
            row[i] = 0;
        }

        int[] col = new int[n];
        for (j = 0; j < n; j++)
        {
            col[j] = -(j + 1);
        }

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (!(Math.Abs(rref[i + j * m] - 1.0) <= double.Epsilon))
                {
                    continue;
                }

                row[i] = j + 1;
                col[j] = j + 1;
                break;
            }
        }

        double[] nullspace = r8mat_zeros_new(n, nullspace_size);

        int j2 = 0;
        //
        //  If column J does not contain a leading 1, then it contains
        //  information about a null vector.
        //
        for (j = 0; j < n; j++)
        {
            switch (col[j])
            {
                case < 0:
                {
                    for (i = 0; i < m; i++)
                    {
                        if (rref[i + j * m] == 0.0)
                        {
                            continue;
                        }

                        int i2 = row[i] - 1;
                        nullspace[i2 + j2 * n] = -rref[i + j * m];
                    }

                    nullspace[j + j2 * n] = 1.0;
                    j2 += 1;
                    break;
                }
            }
        }

        return nullspace;
    }

    public static int r8mat_nullspace_size(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    Let A be an MxN matrix.
        //
        //    If X is an N-vector, and A*X = 0, then X is a null vector of A.
        //
        //    The set of all null vectors of A is called the nullspace of A.
        //
        //    The 0 vector is always in the null space.
        //
        //    If the 0 vector is the only vector in the nullspace of A, then A
        //    is said to have maximum column rank.  (Because A*X=0 can be regarded
        //    as a linear combination of the columns of A).  In particular, if A
        //    is square, and has maximum column rank, it is nonsingular.
        //
        //    The dimension of the nullspace is the number of linearly independent
        //    vectors that span the nullspace.  If A has maximum column rank,
        //    its nullspace has dimension 0.
        //
        //    This routine ESTIMATES the dimension of the nullspace.  Cases of
        //    singularity that depend on exact arithmetic will probably be missed.
        //
        //    The nullspace will be estimated by counting the leading 1's in the
        //    reduced row echelon form of A, and subtracting this from N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input, double A[M*N], the matrix to be analyzed.
        //
        //    Output, int R8MAT_NULLSPACE_SIZE, the estimated size
        //    of the nullspace.
        //
    {
        int i;
        //
        //  Make a copy of A.
        //
        double[] rref = r8mat_copy_new(m, n, a);
        //
        //  Get the reduced row echelon form of A.
        //
        r8mat_rref(m, n, ref rref);
        //
        //  Count the leading 1's in A.
        //
        int leading = 0;
        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (!(Math.Abs(rref[i + j * m] - 1.0) <= double.Epsilon))
                {
                    continue;
                }

                leading += 1;
                break;
            }
        }

        int nullspace_size = n - leading;

        return nullspace_size;
    }

    public static double[] r8mat_orth_uniform_new(int n, ref r8NormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ORTH_UNIFORM_NEW returns a random orthogonal matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The inverse of A is equal to A'.
        //
        //    A * A'  = A' * A = I.
        //
        //    Columns and rows of A have unit Euclidean norm.
        //
        //    Distinct pairs of columns of A are orthogonal.
        //
        //    Distinct pairs of rows of A are orthogonal.
        //
        //    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
        //
        //    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
        //
        //    The determinant of A is +1 or -1.
        //
        //    All the eigenvalues of A have modulus 1.
        //
        //    All singular values of A are 1.
        //
        //    All entries of A are between -1 and 1.
        //
        //  Discussion:
        //
        //    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
        //    National Academy of Sciences of Belarus, for convincingly
        //    pointing out the severe deficiencies of an earlier version of
        //    this routine.
        //
        //    Essentially, the computation involves saving the Q factor of the
        //    QR factorization of a matrix whose entries are normally distributed.
        //    However, it is only necessary to generate this matrix a column at
        //    a time, since it can be shown that when it comes time to annihilate
        //    the subdiagonal elements of column K, these (transformed) elements of
        //    column K are still normally distributed random values.  Hence, there
        //    is no need to generate them at the beginning of the process and
        //    transform them K-1 times.
        //
        //    For computational efficiency, the individual Householder transformations
        //    could be saved, as recommended in the reference, instead of being
        //    accumulated into an explicit matrix format.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Pete Stewart,
        //    Efficient Generation of Random Orthogonal Matrices With an Application
        //    to Condition Estimators,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 3, June 1980, pages 403-409.
        //
        //  Parameters:
        //
        //    Input, int N, the order of A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8MAT_ORTH_UNIFORM_NEW[N*N], the orthogonal matrix.
        //
    {
        int j;
        //
        //  Start with Q = the identity matrix.
        //
        double[] q = r8mat_identity_new(n);
        //
        //  Now behave as though we were computing the QR factorization of
        //  some other random matrix.  Generate the N elements of the first column,
        //  compute the Householder matrix H1 that annihilates the subdiagonal elements,
        //  and set Q := Q * H1' = Q * H.
        //
        //  On the second step, generate the lower N-1 elements of the second column,
        //  compute the Householder matrix H2 that annihilates them,
        //  and set Q := Q * H2' = Q * H2 = H1 * H2.
        //
        //  On the N-1 step, generate the lower 2 elements of column N-1,
        //  compute the Householder matrix HN-1 that annihilates them, and
        //  and set Q := Q * H(N-1)' = Q * H(N-1) = H1 * H2 * ... * H(N-1).
        //  This is our random orthogonal matrix.
        //
        double[] a_col = new double[n];

        for (j = 1; j < n; j++)
        {
            //
            //  Set the vector that represents the J-th column to be annihilated.
            //
            int i;
            for (i = 1; i < j; i++)
            {
                a_col[i - 1] = 0.0;
            }

            for (i = j; i <= n; i++)
            {
                a_col[i - 1] = r8_normal_01(ref data, ref seed);
            }

            //
            //  Compute the vector V that defines a Householder transformation matrix
            //  H(V) that annihilates the subdiagonal elements of A.
            //
            double[] v = r8vec_house_column(n, a_col, j);
            //
            //  Postmultiply the matrix Q by H'(V) = H(V).
            //
            double[] q2 = r8mat_house_axh_new(n, q, v);
                
            r8mat_copy(n, n, q2, ref q);

        }

        return q;
    }

    public static void r8mat_perm0 ( int n, ref double[] a, int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PERM0 permutes the rows and columns of a square R8MAT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 May 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double A[N*N].
        //    On input, the matrix to be permuted.
        //    On output, the permuted matrix.
        //
        //    Input, int P[N], a permutation to be applied to the rows
        //    and columns.  P[I] is the new number of row and column I.
        //
    {
        int i;
        const int iopt = 1;
        int isgn = 0;
        int ncycle = 0;

        Permutation.perm0_cycle ( n, p, ref isgn, ref ncycle, iopt );
        //
        //  Temporarily increment P by 1.
        //
        for ( i = 0; i < n; i++ )
        {
            p[i] += 1;
        }

        for ( i = 1; i <= n; i++ )
        {
            int i1 = - p[i-1];

            switch (i1)
            {
                case > 0:
                {
                    int lc = 0;

                    for ( ; ; )
                    {
                        i1 = p[i1-1];
                        lc += 1;

                        if ( i1 <= 0 )
                        {
                            break;
                        }

                    }

                    i1 = i;

                    int j;
                    for ( j = 1; j <= n; j++ )
                    {
                        switch (p[j-1])
                        {
                            case <= 0:
                            {
                                int j2 = j;
                                int k = lc;

                                for ( ; ; )
                                {
                                    int j1 = j2;
                                    double it = a[i1-1+(j1-1)*n];

                                    for ( ; ; )
                                    {
                                        i1 = Math.Abs ( p[i1-1] );
                                        j1 = Math.Abs ( p[j1-1] );

                                        (a[i1-1+(j1-1)*n], it) = (it, a[i1-1+(j1-1)*n]);

                                        if ( j1 != j2 )
                                        {
                                            continue;
                                        }

                                        k -= 1;

                                        if ( i1 == i )
                                        {
                                            break;
                                        }
                                    }

                                    j2 = Math.Abs ( p[j2-1] );

                                    if ( k == 0 ) 
                                    {
                                        break;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
        //
        //  Restore the positive signs of the data.
        //
        for ( i = 0; i < n; i++ )
        {
            p[i] = Math.Abs ( p[i] );
        }

        for ( i = 0; i < n; i++ )
        {
            p[i] -= 1;
        }
    }

    public static void r8mat_2perm0(int m, int n, ref double[] a, int[] p, int[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_2PERM0 permutes rows and columns of a rectangular R8MAT, in place.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int M, number of rows in the matrix.
        //
        //    Input, int N, number of columns in the matrix.
        //
        //    Input/output, double A[M*N].
        //    On input, the matrix to be permuted.
        //    On output, the permuted matrix.
        //
        //    Input, int P[M], the row permutation.  P(I) is the new number of row I.
        //
        //    Input, int Q[N], the column permutation.  Q(I) is the new number of
        //    column I.  
        //
    {
        int i;
        int is_ = 0;
        int j;
        int nc = 0;
        /*
        Wretched maneuvers to deal with necessity of 1-based values,
        and to handle case where P and Q are same vector.
        */
        int[] p1 = i4vec_copy_new(m, p);
        Permutation.perm0_cycle(m, p1, ref is_, ref nc, 1);
        for (i = 0; i < m; i++)
        {
            p1[i] += 1;
        }

        int[] q1 = i4vec_copy_new(n, q);
        Permutation.perm0_cycle(n, q1, ref is_, ref nc, 1);
        for (j = 0; j < n; j++)
        {
            q1[j] += 1;
        }

        for (i = 1; i <= m; i++)
        {
            int i1 = -p[i - 1];

            switch (i1)
            {
                case > 0:
                {
                    int lc = 0;

                    for (;;)
                    {
                        i1 = p[i1 - 1];
                        lc += 1;

                        if (i1 <= 0)
                        {
                            break;
                        }

                    }

                    i1 = i;

                    for (j = 1; j <= n; j++)
                    {
                        switch (q[j - 1])
                        {
                            case <= 0:
                            {
                                int j2 = j;
                                int k = lc;

                                for (;;)
                                {
                                    int j1 = j2;
                                    double t = a[i1 - 1 + (j1 - 1) * n];

                                    for (;;)
                                    {
                                        i1 = Math.Abs(p[i1 - 1]);
                                        j1 = Math.Abs(q[j1 - 1]);

                                        (a[i1 - 1 + (j1 - 1) * n], t) = (t, a[i1 - 1 + (j1 - 1) * n]);


                                        if (j1 != j2)
                                        {
                                            continue;
                                        }

                                        k -= 1;

                                        if (i1 == i)
                                        {
                                            break;
                                        }
                                    }

                                    j2 = Math.Abs(q[j2 - 1]);

                                    if (k == 0)
                                    {
                                        break;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static double r8mat_permanent(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PERMANENT computes the permanent of an R8MAT.
        //
        //  Discussion:
        //
        //    The permanent function is similar to the determinant.  Recall that
        //    the determinant of a matrix may be defined as the sum of all the
        //    products:
        //
        //      S * A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
        //
        //    where I is any permutation of the columns of the matrix, and S is the
        //    sign of the permutation.  By contrast, the permanent function is
        //    the (unsigned) sum of all the products
        //
        //      A(1,I(1)) * A(2,I(2)) * ... * A(N,I(N))
        //
        //    where I is any permutation of the columns of the matrix.  The only
        //    difference is that there is no permutation sign multiplying each summand.
        //
        //    Symbolically, then, the determinant of a 2 by 2 matrix
        //
        //      a b
        //      c d
        //
        //    is a*d-b*c, whereas the permanent of the same matrix is a*d+b*c.
        //
        //
        //    The permanent is invariant under row and column permutations.
        //    If a row or column of the matrix is multiplied by S, then the
        //      permanent is likewise multiplied by S.
        //    If the matrix is square, then the permanent is unaffected by
        //      transposing the matrix.
        //    Unlike the determinant, however, the permanent does change if
        //      one row is added to another, and it is not true that the
        //      permanent of the product is the product of the permanents.
        //
        //
        //    Note that if A is a matrix of all 1's and 0's, then the permanent
        //    of A counts exactly which permutations hit exactly 1's in the matrix.
        //    This fact can be exploited for various combinatorial purposes.
        //
        //    For instance, setting the diagonal of A to 0, and the offdiagonals
        //    to 1, the permanent of A counts the number of derangements of N
        //    objects.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, number of rows and columns in matrix.
        //
        //    Input, double A[N*N], the matrix whose permanent is desired.
        //
        //    Output, double R8MAT_PERMANENT, the value of the permanent of A.
        //
    {
        int i;
        int iadd = 0;
        int ncard = 0;

        bool more = false;

        int[] iwork = new int[n];
        double[] work = new double[n];

        for (i = 1; i <= n; i++)
        {
            work[i - 1] = a[i - 1 + (n - 1) * n];
            int j;
            for (j = 1; j <= n; j++)
            {
                work[i - 1] -= 0.5 * a[i - 1 + (j - 1) * n];
            }
        }

        double p = 0.0;
        double sgn = -1.0;

        for (;;)
        {
            sgn = -sgn;

            Subset.subset_gray_next(n - 1, ref iwork, ref more, ref ncard, ref iadd);

            if (ncard != 0)
            {
                double z = 2 * iwork[iadd - 1] - 1;
                for (i = 1; i <= n; i++)
                {
                    work[i - 1] += z * a[i - 1 + (iadd - 1) * n];
                }
            }

            double prod = 1.0;
            for (i = 0; i < n; i++)
            {
                prod *= work[i];
            }

            p += sgn * prod;

            if (!more)
            {
                break;
            }
        }

        double perm = p * (4 * (n % 2) - 2);

        return perm;
    }

    public static void r8mat_plot(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PLOT "plots" an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        int jlo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        for (jlo = 1; jlo <= n; jlo += 70)
        {
            int jhi = Math.Min(jlo + 70 - 1, n);
            Console.WriteLine("");
            string cout = "          ";
            int j;
            for (j = jlo; j <= jhi; j++)
            {
                cout += j % 10;
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            int i;
            for (i = 1; i <= m; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "    ";
                for (j = jlo; j <= jhi; j++)
                {
                    cout += r8mat_plot_symbol(a[i - 1 + (j - 1) * m]);
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static char r8mat_plot_symbol(double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PLOT_SYMBOL returns a symbol for entries of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, a value whose symbol is desired.
        //
        //    Output, char R8MAT_PLOT_SYMBOL, is
        //    '-' if R is negative,
        //    '0' if R is zero,
        //    '+' if R is positive.
        //
    {
        char c = r switch
        {
            < 0.0 => '-',
            0.0 => '0',
            _ => '+'
        };

        return c;
    }

    public static double[] r8mat_poly_char(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, double A[N*N], the N by N matrix.
        //
        //    Output, double R8MAT_POLY_CHAR[N+1], the coefficients of the characteristic
        //    polynomial of A.  P(N) contains the coefficient of X^N
        //    (which will be 1), P(I) contains the coefficient of X^I,
        //    and P(0) contains the constant term.
        //
    {
        int order;

        double[] p = new double[n + 1];
        //
        //  Initialize WORK1 to the identity matrix.
        //
        double[] work1 = r8mat_identity_new(n);

        p[n] = 1.0;

        for (order = n - 1; 0 <= order; order--)
        {
            //
            //  Work2 = A * WORK1.
            //
            double[] work2 = r8mat_mm_new(n, n, n, a, work1);
            //
            //  Take the trace.
            //
            double trace = r8mat_trace(n, work2);
            //
            //  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
            //
            p[order] = -trace / (n - order);
            //
            //  WORK1 := WORK2 + P(IORDER) * Identity.
            //

            r8mat_copy(n, n, work2, ref work1);

            int i;
            for (i = 0; i < n; i++)
            {
                work1[i + i * n] += p[order];
            }
        }

        return p;
    }

    public static double[] r8mat_power(int n, double[] a, int npow)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_POWER computes a nonnegative power of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The algorithm is:
        //
        //      B = I
        //      do NPOW times:
        //        B = A * B
        //      end
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of A.
        //
        //    Input, double A[N*N], the matrix to be raised to a power.
        //
        //    Input, int NPOW, the power to which A is to be raised.
        //    NPOW must be nonnegative.
        //
        //    Output, double B[N*N], the value of A^NPOW.
        //
    {
        int ipow;

        switch (npow)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("R8MAT_POWER - Fatal error!");
                Console.WriteLine("  Input value of NPOW < 0.");
                Console.WriteLine("  NPOW = " + npow + "");
                return null;
        }

        double[] b = r8mat_identity_new(n);

        for (ipow = 1; ipow <= npow; ipow++)
        {
            double[] c = r8mat_mm_new(n, n, n, a, b);
            r8mat_copy(n, n, c, ref b);
        }

        return b;
    }

    public static void r8mat_power_method(int n, double[] a, ref double r, ref double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_POWER_METHOD applies the power method to a matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    If the power method has not converged, then calling the routine
        //    again immediately with the output from the previous call will
        //    continue the iteration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of A.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Output, double *R, the estimated eigenvalue.
        //
        //    Input/output, double V[N], on input, an estimate
        //    for the eigenvector.  On output, an improved estimate for the
        //    eigenvector.
        //
    {
        int i;
        int it;
        const double it_eps = 0.0001;
        const int it_max = 100;
        const int it_min = 10;

        double eps = Math.Sqrt(r8_epsilon());

        r = r8vec_norm(n, v);

        switch (r)
        {
            case 0.0:
            {
                for (i = 0; i < n; i++)
                {
                    v[i] = 1.0;
                }

                r = Math.Sqrt(n);
                break;
            }
        }

        for (i = 0; i < n; i++)
        {
            v[i] /= r;
        }

        for (it = 1; it <= it_max; it++)
        {
            double[] av = r8mat_mv_new(n, n, a, v);

            double r_old = r;
            r = r8vec_norm(n, av);

            if (it_min < it)
            {
                if (Math.Abs(r - r_old) <= it_eps * (1.0 + Math.Abs(r)))
                {
                    break;
                }
            }

            r8vec_copy(n, av, ref v);

            if (r != 0.0)
            {
                for (i = 0; i < n; i++)
                {
                    v[i] /= r;
                }
            }

            //
            //  Perturb V a bit, to avoid cases where the initial guess is exactly
            //  the eigenvector of a smaller eigenvalue.
            //
            if (it >= it_max / 2)
            {
                continue;
            }

            int j = (it - 1) % n;
            v[j] += eps * (1.0 + Math.Abs(v[j]));
            double r2 = r8vec_norm(n, v);
            for (i = 0; i < n; i++)
            {
                v[i] /= r2;
            }
        }
    }

    public static double r8mat_product_elementwise(int m, int n, ref double[] a, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRODUCT_ELEMENTWISE returns the elementwise produce to two R8MAT's.
        //
        //  Example:
        //
        //    A = [ 1, 2, 3;    B = [ 1, 3, 5;    product = 86
        //         4, 5, 6 ]         2, 4, 6 ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows.
        //
        //    Input, int N, the number of columns.
        //
        //    Input, double A[M*N], B[M*N], the two matrices.
        //
        //    Output, double I4MAT_PRODUCT_ELEMENTWISE, the elementwise 
        //    product of A and B.
        //
    {
        int j;

        double value = 0.0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                value += a[i + j * m] * b[i + j * m];
            }
        }

        return value;
    }

    public static double r8mat_ref(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_REF computes the row echelon form of a matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    A matrix is in row echelon form if:
        //
        //    * The first nonzero entry in each row is 1.
        //
        //    * The leading 1 in a given row occurs in a column to
        //      the right of the leading 1 in the previous row.
        //
        //    * Rows which are entirely zero must occur last.
        //
        //  Example:
        //
        //    Input matrix:
        //
        //     1.0  3.0  0.0  2.0  6.0  3.0  1.0
        //    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
        //     3.0  9.0  0.0  0.0  6.0  6.0  2.0
        //    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
        //
        //    Output matrix:
        //
        //     1.0  3.0  0.0  2.0  6.0  3.0  1.0
        //     0.0  0.0  0.0  1.0  2.0  4.5  1.5
        //     0.0  0.0  0.0  0.0  0.0  1.0  0.3
        //     0.0  0.0  0.0  0.0  0.0  0.0  0.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Cullen,
        //    An Introduction to Numerical Linear Algebra,
        //    PWS Publishing Company, 1994,
        //    ISBN: 978-0534936903,
        //    LC: QA185.D37.C85.
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input/output, double A[M*N].  On input, the matrix to be
        //    analyzed.  On output, the REF form of the matrix.
        //
        //    Output, double R8MAT_REF, the pseudo-determinant.
        //
    {
        int i;
        int j;
        int r;

        double det = 1.0;
        double asum = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                asum += Math.Abs(a[i + j * m]);
            }
        }

        double tol = r8_epsilon() * asum;
        int lead = 0;

        for (r = 0; r < m; r++)
        {
            if (n - 1 < lead)
            {
                break;
            }

            i = r;

            while (Math.Abs(a[i + lead * m]) <= tol)
            {
                i += 1;

                if (m - 1 >= i)
                {
                    continue;
                }

                i = r;
                lead += 1;
                if (n - 1 >= lead)
                {
                    continue;
                }

                lead = -1;
                break;
            }

            if (lead < 0)
            {
                break;
            }

            double temp;
            for (j = 0; j < n; j++)
            {
                temp = a[i + j * m];
                a[i + j * m] = a[r + j * m];
                a[r + j * m] = temp;
            }

            det *= a[r + lead * m];
            temp = a[r + lead * m];

            for (j = 0; j < n; j++)
            {
                a[r + j * m] /= temp;
            }

            for (i = r + 1; i < m; i++)
            {
                temp = a[i + lead * m];
                for (j = 0; j < n; j++)
                {
                    a[i + j * m] -= temp * a[r + j * m];
                }
            }

            lead += 1;
        }

        return det;
    }

    public static double r8mat_rms(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_RMS returns the RMS norm of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8's.
        //
        //    The matrix RMS norm is defined as:
        //
        //      R8MAT_RMS =
        //        sqrt ( sum ( 0 <= J < N ) sum ( 0 <= I < M ) A[I,J]^2 / M / N ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the dimensions of the array.
        //
        //    Input, double A[M*N], the array.
        //
        //    Output, double R8MAT_RMS, the RMS norm of A.
        //
    {
        int j;

        double value = 0.0;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                value += a[i + j * m] * a[i + j * m];
            }

            value = Math.Sqrt(value / m / n);
        }

        return value;
    }

    public static double r8mat_rref(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_RREF computes the reduced row echelon form of a matrix.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    A matrix is in row echelon form if:
        //
        //    * The first nonzero entry in each row is 1.
        //
        //    * The leading 1 in a given row occurs in a column to
        //      the right of the leading 1 in the previous row.
        //
        //    * Rows which are entirely zero must occur last.
        //
        //    The matrix is in reduced row echelon form if, in addition to
        //    the first three conditions, it also satisfies:
        //
        //    * Each column containing a leading 1 has no other nonzero entries.
        //
        //  Example:
        //
        //    Input matrix:
        //
        //     1.0  3.0  0.0  2.0  6.0  3.0  1.0
        //    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
        //     3.0  9.0  0.0  0.0  6.0  6.0  2.0
        //    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0
        //
        //    Output matrix:
        //
        //     1.0  3.0  0.0  0.0  2.0  0.0  0.0
        //     0.0  0.0  0.0  1.0  2.0  0.0  0.0
        //     0.0  0.0  0.0  0.0  0.0  1.0  0.3
        //     0.0  0.0  0.0  0.0  0.0  0.0  0.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Cullen,
        //    An Introduction to Numerical Linear Algebra,
        //    PWS Publishing Company, 1994,
        //    ISBN: 978-0534936903,
        //    LC: QA185.D37.C85.
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input/output, double A[M*N].  On input, the matrix to be
        //    analyzed.  On output, the RREF form of the matrix.
        //
        //    Output, double R8MAT_RREF, the pseudo-determinant.
        //
    {
        int i;
        int j;
        int r;

        double det = 1.0;
        double asum = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                asum += Math.Abs(a[i + j * m]);
            }
        }

        double tol = r8_epsilon() * asum;
        int lead = 0;

        for (r = 0; r < m; r++)
        {
            if (n - 1 < lead)
            {
                break;
            }

            i = r;

            while (Math.Abs(a[i + lead * m]) <= tol)
            {
                i += 1;

                if (m - 1 >= i)
                {
                    continue;
                }

                i = r;
                lead += 1;
                if (n - 1 >= lead)
                {
                    continue;
                }

                lead = -1;
                break;
            }

            if (lead < 0)
            {
                break;
            }

            double temp;
            for (j = 0; j < n; j++)
            {
                temp = a[i + j * m];
                a[i + j * m] = a[r + j * m];
                a[r + j * m] = temp;
            }

            det *= a[r + lead * m];
            temp = a[r + lead * m];

            for (j = 0; j < n; j++)
            {
                a[r + j * m] /= temp;
            }

            for (i = 0; i < m; i++)
            {
                if (i == r)
                {
                    continue;
                }

                temp = a[i + lead * m];
                for (j = 0; j < n; j++)
                {
                    a[i + j * m] -= temp * a[r + j * m];
                }
            }

            lead += 1;

        }

        return det;
    }
    public static double r8mat_trace(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_TRACE computes the trace of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The trace of a square matrix is the sum of the diagonal elements.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, double A[N*N], the matrix whose trace is desired.
        //
        //    Output, double R8MAT_TRACE, the trace of the matrix.
        //
    {
        int i;

        double value = 0.0;
        for (i = 0; i < n; i++)
        {
            value += a[i + i * n];
        }

        return value;
    }
        
    public static double r8mat_vtmv ( int m, int n, double[] x, double[] a, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_VTMV multiplies computes the scalar x' * A * y.
        //
        //  Discussion:
        //
        //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix.
        //
        //    Input, double X[N], the first vector factor.
        //
        //    Input, double A[M*N], the M by N matrix.
        //
        //    Input, double Y[M], the second vector factor.
        //
        //    Output, double R8MAT_VTMV, the value of X' * A * Y.
        //
    {
        int j;

        double vtmv = 0.0;
        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                vtmv += x[i] * a[i+j*m] * y[j];
            }
        }
        return vtmv;
    }

    public static double r8mat_residual_norm ( int m, int n, double[] a, double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_RESIDUAL_NORM returns the norm of A*x-b.
        //
        //  Discussion:
        //
        //    A is an MxN R8MAT, a matrix of R8's.
        //
        //    X is an N R8VEC, and B is an M R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[M,N], the M by N matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the right hand side vector.
        //
        //    Output, double R8MAT_RESIDUAL_NORM, the norm of A*x-b.
        //
    {
        int i;

        double[] r = new double[m];

        for ( i = 0; i < m; i++ )
        {
            r[i] = - b[i];
            int j;
            for ( j = 0; j < n; j++ )
            {
                r[i] += a[i+j*m] * x[j];
            }
        }

        double r_norm = 0.0;
        for ( i = 0; i < m; i++ )
        {
            r_norm += r[i] * r[i];
        }
        r_norm = Math.Sqrt ( r_norm );
            
        return r_norm;
    }

        
    public static void r8mat_add ( int m, int n, double[] a, ref double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ADD adds one R8MAT to another.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the matrix to add.
        //
        //    Input/output, double B[M*N], the matrix to be incremented.
        //
    {
        int j;

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                b[i+j*m] += a[i+j*m];
            }
        }
    }

    public static void r8mat_add(int m, int n, double alpha, double[] a, double beta,
            double[] b, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 November 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double ALPHA, the multiplier for A.
        //
        //    Input, double A[M*N], the first matrix.
        //
        //    Input, double BETA, the multiplier for A.
        //
        //    Input, double B[M*N], the second matrix.
        //
        //    Output, double C[M*N], the sum of alpha*A+beta*B.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                c[i + j * m] = alpha * a[i + j * m] + beta * b[i + j * m];
            }
        }
    }

    public static double[] r8mat_add_new(int m, int n, double alpha, double[] a, double beta,
            double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_ADD_NEW computes C = alpha * A + beta * B for R8MAT's.
        //
        //  Discussion:
        //
        //    An R8MAT is an array of R8 values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double ALPHA, the multiplier for A.
        //
        //    Input, double A[M*N], the first matrix.
        //
        //    Input, double BETA, the multiplier for A.
        //
        //    Input, double B[M*N], the second matrix.
        //
        //    Output, double R8MAT_ADD_NEW[M*N], the sum of alpha*A+beta*B.
        //
    {
        int j;

        double[] c = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                c[i + j * m] = alpha * a[i + j * m] + beta * b[i + j * m];
            }
        }

        return c;
    }

    public static void r8mat_divide ( int m, int n, double s, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIVIDE divides an R8MAT by a scalar.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double S, the divisor
        //
        //    Input/output, double A[M*N], the matrix to be scaled.
        //
    {
        int j;

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                a[i+j*m] /= s;
            }
        }
    }
        
    public static double r8mat_dif_fro ( int m, int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DIF_FRO returns the Frobenius norm of the difference of R8MAT's.
        //
        //  Discussion: 							    
        //
        //    An R8MAT is a doubly dimensioned array of double precision values, which
        //    may be stored as a vector in column-major order.
        //
        //    The Frobenius norm is defined as
        //
        //      R8MAT_NORM_FRO = sqrt (
        //        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
        //
        //    The matrix Frobenius norm is not derived from a vector norm, but
        //    is compatible with the vector L2 norm, so that:
        //
        //      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], double B[M*N], the matrices for which we
        //    want the Frobenius norm of the difference.
        //
        //    Output, double R8MAT_DIF_FRO, the Frobenius norm of ( A - B ).
        //
    {
        int j;

        double value = 0.0;
        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                value += Math.Pow ( a[i+j*m] - b[i+j*m], 2 );
            }
        }
        value = Math.Sqrt ( value );

        return value;
    }
        
    public static void r8mat_latinize ( int m, int n, ref double[] table )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_LATINIZE "Latinizes" an R8MAT.
        //
        //  Discussion:
        //
        //    It is assumed, though not necessary, that the input dataset
        //    has points that lie in the unit hypercube.
        //
        //    In any case, the output dataset will have this property.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 December 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of cells.
        //
        //    Input/output, double TABLE[M*N].  On input, the dataset to
        //    be "Latinized".  On output, the Latinized dataset.
        //
    {
        int i;

        double[] column = new double[n];

        for ( i = 0; i < m; i++ )
        {
            int j;
            for ( j = 0; j < n; j++ )
            {
                column[j] = table[i+j*m];
            }
            int[] indx = r8vec_sort_heap_index_a ( n, column );

            for ( j = 0; j < n; j++ )
            {
                table[i+indx[j]*m] = (2 * j + 1) / ( double ) ( 2 * n );
            }
        }
    }
        
    public static double[] r8mat_gen ( int lda, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_GEN generates a random R8MAT..
        //
        //  Modified:
        //
        //    06 June 2005
        //
        //  Parameters:
        //
        //    Input, integer LDA, the leading dimension of the matrix.
        //
        //    Input, integer N, the order of the matrix.
        //
        //    Output, double R8MAT_GEN[LDA*N], the N by N matrix.
        //
    {
        int[] init = { 1, 2, 3, 1325 };
        int j;

        double[] a = new double[lda*n];

        for ( j = 1; j <= n; j++ )
        {
            int i;
            for ( i = 1; i <= n; i++ )
            {
                a[i-1+(j-1)*lda] = r8_random ( init ) - 0.5;
            }
        }

        return a;
    }
}