using System;
using System.IO;
using System.Linq;
using Burkardt.Table;

namespace Burkardt.Types
{
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
            double[] a;
            int i;
            int j;

            a = new double[n*n];

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    if ( j == 0 && x[i] == 0.0 )
                    {
                        a[i+j*n] = 1.0;
                    }
                    else
                    {
                        a[i+j*n] = Math.Pow ( x[i], j );
                    }
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
            int j;
            int k;
            const double r8_inf = 1.0E+30;
            double s;

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    if ( m[j+i*n] < r8_inf )
                    {
                        for ( k = 0; k < n; k++ )
                        {
                            if ( m[i+k*n] < r8_inf )
                            {
                                s = m[j+i*n] + m[i+k*n];
                                if ( s < m[j+k*n] )
                                {
                                    m[j+k*n] = s;
                                }
                            }
                        }
                    }
                }
            }
            return;
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
            double[] c;
            int i;
            int j;

            c = new double[m*n];

            for ( j = 0; j < n; j++ )
            {
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
            int j;
            int k;

            for ( i = 0; i < n1; i ++ )
            {
                for ( j = 0; j < n3; j++ )
                {
                    c[i+j*n1] = 0.0;
                    for ( k = 0; k < n2; k++ )
                    {
                        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
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
            double[] c;
            int i;
            int j;
            int k;

            c = new double[n1 * n3];

            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n3; j++)
                {
                    c[i + j * n1] = 0.0;
                    for (k = 0; k < n2; k++)
                    {
                        c[i + j * n1] = c[i + j * n1] + a[i + k * n1] * b[k + j * n2];
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
            int i;
            int j;
            double value;

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + Math.Pow(a[i + j * m] - b[i + j * m], 2);
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
            double eps;
            double f00;
            double fmm;
            double fmp;
            double fpm;
            double fpp;
            double[] h;
            int i;
            int j;
            double[] s;
            double xi;
            double xj;
            //
            //  Choose the stepsizes.
            //
            s = new double[n];

            eps = Math.Pow(double.Epsilon, 0.33);

            for (i = 0; i < n; i++)
            {
                s[i] = eps * Math.Max(Math.Abs(x[i]), 1.0);
            }

            //
            //  Calculate the diagonal elements.
            //
            h = new double[n * n];

            for (i = 0; i < n; i++)
            {
                xi = x[i];

                f00 = fx(n, x);

                x[i] = xi + s[i];
                fpp = fx(n, x);

                x[i] = xi - s[i];
                fmm = fx(n, x);

                h[i + i * n] = ((fpp - f00) + (fmm - f00)) / s[i] / s[i];

                x[i] = xi;
            }

            //
            //  Calculate the off diagonal elements.
            //
            for (i = 0; i < n; i++)
            {
                xi = x[i];

                for (j = i + 1; j < n; j++)
                {
                    xj = x[j];

                    x[i] = xi + s[i];
                    x[j] = xj + s[j];
                    fpp = fx(n, x);

                    x[i] = xi + s[i];
                    x[j] = xj - s[j];
                    fpm = fx(n, x);

                    x[i] = xi - s[i];
                    x[j] = xj + s[j];
                    fmp = fx(n, x);

                    x[i] = xi - s[i];
                    x[j] = xj - s[j];
                    fmm = fx(n, x);

                    h[j + i * n] = ((fpp - fpm) + (fmm - fmp)) / (4.0 * s[i] * s[j]);

                    h[i + j * n] = h[j + i * n];

                    x[j] = xj;
                }

                x[i] = xi;
            }

            return x;
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
            double[] a;
            int fac;
            int i;
            int j;

            a = new double[m * n];

            fac = (int) Math.Pow(10, Math.Log10(n) + 1);

            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    a[i - 1 + (j - 1) * m] = (double) (fac * i + j);
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
            double del;
            double[] fprime;
            int i;
            int j;
            double xsave;
            double[] work1;
            double[] work2;

            fprime = new double[m * n];
            //
            //  Evaluate the function at the base point, X.
            //
            work2 = fx(m, n, x);
            //
            //  Now, one by one, vary each component J of the base point X, and
            //  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
            //
            for (j = 0; j < n; j++)
            {
                xsave = x[j];
                del = eps * (1.0 + Math.Abs(x[j]));
                x[j] = x[j] + del;
                work1 = fx(m, n, x);
                x[j] = xsave;
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
            double[] c;
            int i;
            int i0;
            int i1;
            int i2;
            int j;
            int j0;
            int j1;
            int j2;
            int m;
            int n;

            m = m1 * m2;
            n = n1 * n2;
            c = new double[m * n];

            for (j1 = 0; j1 < n1; j1++)
            {
                for (i1 = 0; i1 < m1; i1++)
                {
                    i0 = i1 * m2;
                    j0 = j1 * n2;
                    j = j0;
                    for (j2 = 0; j2 < n2; j2++)
                    {
                        i = i0;
                        for (i2 = 0; i2 < m2; i2++)
                        {
                            c[i + j * m] = a[i1 + j1 * m1] * b[i2 + j2 * m2];
                            i = i + 1;
                        }

                        j = j + 1;
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
            int ipiv;
            int j;
            int k;
            double pivot;
            double t;
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
                pivot = 0.0;
                ipiv = -1;

                for (i = j; i < m; i++)
                {
                    if (pivot < Math.Abs(u[i + j * m]))
                    {
                        pivot = Math.Abs(u[i + j * m]);
                        ipiv = i;
                    }
                }

                //
                //  Unless IPIV is zero, swap rows J and IPIV.
                //
                if (ipiv != -1)
                {
                    for (k = 0; k < n; k++)
                    {
                        t = u[j + k * m];
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
                        if (u[i + j * m] != 0.0)
                        {
                            l[i + j * m] = u[i + j * m] / u[j + j * m];

                            u[i + j * m] = 0.0;

                            for (k = j + 1; k < n; k++)
                            {
                                u[i + k * m] = u[i + k * m] - l[i + j * m] * u[j + k * m];
                            }
                        }
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
            double[] c;
            int i;
            int j;
            int k;

            c = new double[n1 * n3];

            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n3; j++)
                {
                    c[i + j * n1] = 0.0;
                    for (k = 0; k < n2; k++)
                    {
                        c[i + j * n1] = c[i + j * n1] + a[i + k * n1] * b[j + k * n3];
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
            double[] c;
            int i;
            int j;
            int k;

            c = new double[n1 * n3];

            for (i = 0; i < n1; i++)
            {
                for (j = 0; j < n3; j++)
                {
                    c[i + j * n1] = 0.0;
                    for (k = 0; k < n2; k++)
                    {
                        c[i + j * n1] = c[i + j * n1] + a[k + i * n2] * b[k + j * n2];
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
            int i;
            int j;
            int s;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (a[i + j * m] < 0.0)
                    {
                        s = -1;
                    }
                    else
                    {
                        s = 1;
                    }

                    a[i + j * m] = s * (int) (Math.Abs(a[i + j * m]) + 0.5);
                }
            }

            return;
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
            int i;
            int j;
            int value;

            value = 0;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (a[i + j * m] != 0.0)
                    {
                        value = value + 1;
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
            int[] col;
            int i;
            int i2;
            int j;
            int j2;
            double[] nullspace;
            int[] row;
            double[] rref;
            //
            //  Make a copy of A.
            //
            rref = r8mat_copy_new(m, n, a);
            //
            //  Get the reduced row echelon form of A.
            //
            r8mat_rref(m, n, ref rref);
            //
            //  Note in ROW the columns of the leading nonzeros.
            //  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
            //
            row = new int[m];
            for (i = 0; i < m; i++)
            {
                row[i] = 0;
            }

            col = new int[n];
            for (j = 0; j < n; j++)
            {
                col[j] = -(j + 1);
            }

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (rref[i + j * m] == 1.0)
                    {
                        row[i] = (j + 1);
                        col[j] = (j + 1);
                        break;
                    }
                }
            }

            nullspace = r8mat_zeros_new(n, nullspace_size);

            j2 = 0;
            //
            //  If column J does not contain a leading 1, then it contains
            //  information about a null vector.
            //
            for (j = 0; j < n; j++)
            {
                if (col[j] < 0)
                {
                    for (i = 0; i < m; i++)
                    {
                        if (rref[i + j * m] != 0.0)
                        {
                            i2 = row[i] - 1;
                            nullspace[i2 + j2 * n] = -rref[i + j * m];
                        }
                    }

                    nullspace[j + j2 * n] = 1.0;
                    j2 = j2 + 1;
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
            int j;
            int leading;
            int nullspace_size;
            double[] rref;
            //
            //  Make a copy of A.
            //
            rref = r8mat_copy_new(m, n, a);
            //
            //  Get the reduced row echelon form of A.
            //
            r8mat_rref(m, n, ref rref);
            //
            //  Count the leading 1's in A.
            //
            leading = 0;
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (rref[i + j * m] == 1.0)
                    {
                        leading = leading + 1;
                        break;
                    }
                }
            }

            nullspace_size = n - leading;

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
            double[] a_col;
            double[] q;
            double[] q2;
            int i;
            int j;
            double[] v;
            //
            //  Start with Q = the identity matrix.
            //
            q = r8mat_identity_new(n);
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
            a_col = new double[n];

            for (j = 1; j < n; j++)
            {
                //
                //  Set the vector that represents the J-th column to be annihilated.
                //
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
                v = r8vec_house_column(n, a_col, j);
                //
                //  Postmultiply the matrix Q by H'(V) = H(V).
                //
                q2 = r8mat_house_axh_new(n, q, v);
                
                r8mat_copy(n, n, q2, ref q);

            }

            return q;
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
            int i;
            int j;
            int jhi;
            int jlo;

            Console.WriteLine("");
            Console.WriteLine(title + "");

            for (jlo = 1; jlo <= n; jlo = jlo + 70)
            {
                jhi = Math.Min(jlo + 70 - 1, n);
                Console.WriteLine("");
                string cout = "          ";
                for (j = jlo; j <= jhi; j++)
                {
                    cout += (j % 10);
                }

                Console.WriteLine(cout);
                Console.WriteLine("");

                for (i = 1; i <= m; i++)
                {
                    cout = i.ToString().PadLeft(6) + "    ";
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
            char c;

            if (r < 0.0)
            {
                c = '-';
            }
            else if (r == 0.0)
            {
                c = '0';
            }
            else
            {
                c = '+';
            }

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
            int i;
            int order;
            double[] p;
            double trace;
            double[] work1;
            double[] work2;

            p = new double[n + 1];
            //
            //  Initialize WORK1 to the identity matrix.
            //
            work1 = r8mat_identity_new(n);

            p[n] = 1.0;

            for (order = n - 1; 0 <= order; order--)
            {
                //
                //  Work2 = A * WORK1.
                //
                work2 = r8mat_mm_new(n, n, n, a, work1);
                //
                //  Take the trace.
                //
                trace = r8mat_trace(n, work2);
                //
                //  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
                //
                p[order] = -trace / (double) (n - order);
                //
                //  WORK1 := WORK2 + P(IORDER) * Identity.
                //

                r8mat_copy(n, n, work2, ref work1);

                for (i = 0; i < n; i++)
                {
                    work1[i + i * n] = work1[i + i * n] + p[order];
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
            double[] b;
            double[] c;
            int ipow;

            if (npow < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8MAT_POWER - Fatal error!");
                Console.WriteLine("  Input value of NPOW < 0.");
                Console.WriteLine("  NPOW = " + npow + "");
                return null;
            }

            b = r8mat_identity_new(n);

            for (ipow = 1; ipow <= npow; ipow++)
            {
                c = r8mat_mm_new(n, n, n, a, b);
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
            double[] av;
            double eps;
            int i;
            int it;
            double it_eps = 0.0001;
            int it_max = 100;
            int it_min = 10;
            int j;
            double r2;
            double r_old;

            eps = Math.Sqrt(double.Epsilon);

            r = r8vec_norm(n, v);

            if (r == 0.0)
            {
                for (i = 0; i < n; i++)
                {
                    v[i] = 1.0;
                }

                r = Math.Sqrt((double) n);
            }

            for (i = 0; i < n; i++)
            {
                v[i] = v[i] / r;
            }

            for (it = 1; it <= it_max; it++)
            {
                av = r8mat_mv_new(n, n, a, v);

                r_old = r;
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
                        v[i] = v[i] / r;
                    }
                }

                //
                //  Perturb V a bit, to avoid cases where the initial guess is exactly
                //  the eigenvector of a smaller eigenvalue.
                //
                if (it < it_max / 2)
                {
                    j = ((it - 1) % n);
                    v[j] = v[j] + eps * (1.0 + Math.Abs(v[j]));
                    r2 = r8vec_norm(n, v);
                    for (i = 0; i < n; i++)
                    {
                        v[i] = v[i] / r2;
                    }
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
            int i;
            int j;
            double value;

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + a[i + j * m] * b[i + j * m];
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
            double asum;
            double det;
            int i;
            int j;
            int lead;
            int r;
            double temp;
            double tol;

            det = 1.0;
            asum = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    asum = asum + Math.Abs(a[i + j * m]);
                }
            }

            tol = double.Epsilon * asum;
            lead = 0;

            for (r = 0; r < m; r++)
            {
                if (n - 1 < lead)
                {
                    break;
                }

                i = r;

                while (Math.Abs(a[i + lead * m]) <= tol)
                {
                    i = i + 1;

                    if (m - 1 < i)
                    {
                        i = r;
                        lead = lead + 1;
                        if (n - 1 < lead)
                        {
                            lead = -1;
                            break;
                        }
                    }
                }

                if (lead < 0)
                {
                    break;
                }

                for (j = 0; j < n; j++)
                {
                    temp = a[i + j * m];
                    a[i + j * m] = a[r + j * m];
                    a[r + j * m] = temp;
                }

                det = det * a[r + lead * m];
                temp = a[r + lead * m];

                for (j = 0; j < n; j++)
                {
                    a[r + j * m] = a[r + j * m] / temp;
                }

                for (i = r + 1; i < m; i++)
                {
                    temp = a[i + lead * m];
                    for (j = 0; j < n; j++)
                    {
                        a[i + j * m] = a[i + j * m] - temp * a[r + j * m];
                    }
                }

                lead = lead + 1;
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
            int i;
            int j;
            double value;

            value = 0.0;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + a[i + j * m] * a[i + j * m];
                }

                value = Math.Sqrt(value / (double) (m) / (double) (n));
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
            double asum;
            double det;
            int i;
            int j;
            int lead;
            int r;
            double temp;
            double tol;

            det = 1.0;
            asum = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    asum = asum + Math.Abs(a[i + j * m]);
                }
            }

            tol = double.Epsilon * asum;
            lead = 0;

            for (r = 0; r < m; r++)
            {
                if (n - 1 < lead)
                {
                    break;
                }

                i = r;

                while (Math.Abs(a[i + lead * m]) <= tol)
                {
                    i = i + 1;

                    if (m - 1 < i)
                    {
                        i = r;
                        lead = lead + 1;
                        if (n - 1 < lead)
                        {
                            lead = -1;
                            break;
                        }
                    }
                }

                if (lead < 0)
                {
                    break;
                }

                for (j = 0; j < n; j++)
                {
                    temp = a[i + j * m];
                    a[i + j * m] = a[r + j * m];
                    a[r + j * m] = temp;
                }

                det = det * a[r + lead * m];
                temp = a[r + lead * m];

                for (j = 0; j < n; j++)
                {
                    a[r + j * m] = a[r + j * m] / temp;
                }

                for (i = 0; i < m; i++)
                {
                    if (i != r)
                    {
                        temp = a[i + lead * m];
                        for (j = 0; j < n; j++)
                        {
                            a[i + j * m] = a[i + j * m] - temp * a[r + j * m];
                        }
                    }
                }

                lead = lead + 1;

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
            double value;

            value = 0.0;
            for (i = 0; i < n; i++)
            {
                value = value + a[i + i * n];
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
            int i;
            int j;
            double vtmv;

            vtmv = 0.0;
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    vtmv = vtmv + x[i] * a[i+j*m] * y[j];
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
            int j;
            double[] r;
            double r_norm;

            r = new double[m];

            for ( i = 0; i < m; i++ )
            {
                r[i] = - b[i];
                for ( j = 0; j < n; j++ )
                {
                    r[i] = r[i] + a[i+j*m] * x[j];
                }
            }

            r_norm = 0.0;
            for ( i = 0; i < m; i++ )
            {
                r_norm = r_norm + r[i] * r[i];
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
            int i;
            int j;

            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    b[i+j*m] = b[i+j*m] + a[i+j*m];
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
            int i;
            int j;

            for (j = 0; j < n; j++)
            {
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
            double[] c;
            int i;
            int j;

            c = new double[m * n];

            for (j = 0; j < n; j++)
            {
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
            int i;
            int j;

            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    a[i+j*m] = a[i+j*m] / s;
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
            int i;
            int j;
            double value;

            value = 0.0;
            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    value = value + Math.Pow ( a[i+j*m] - b[i+j*m], 2 );
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
            double[] column;
            int i;
            int[] indx;
            int j;

            column = new double[n];

            for ( i = 0; i < m; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    column[j] = table[i+j*m];
                }
                indx = r8vec_sort_heap_index_a ( n, column );

                for ( j = 0; j < n; j++ )
                {
                    table[i+indx[j]*m] = ( double ) ( 2 * j + 1 ) / ( double ) ( 2 * n );
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
            double[] a;
            int i;
            int[] init = { 1, 2, 3, 1325 };
            int j;

            a = new double[lda*n];

            for ( j = 1; j <= n; j++ )
            {
                for ( i = 1; i <= n; i++ )
                {
                    a[i-1+(j-1)*lda] = r8_random ( init ) - 0.5;
                }
            }

            return a;
        }
    }
}