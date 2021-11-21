using System;
using Burkardt.Types;

namespace Burkardt.BLAS;

public static class BLAS2D
{
    public static void dgemv(char trans, int m, int n, double alpha, double[] a, int lda,
            double[] x, int incx, double beta, ref double[] y, int incy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEMV computes y := alpha * A * x + beta * y for general matrix A.
        //
        //  Discussion:
        //
        //    DGEMV performs one of the matrix-vector operations
        //      y := alpha*A *x + beta*y
        //    or
        //      y := alpha*A'*x + beta*y,
        //    where alpha and beta are scalars, x and y are vectors and A is an
        //    m by n matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 April 2014
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Parameters:
        //
        //    Input, char TRANS, specifies the operation to be performed:
        //    'n' or 'N'   y := alpha*A *x + beta*y.
        //    't' or 'T'   y := alpha*A'*x + beta*y.
        //    'c' or 'C'   y := alpha*A'*x + beta*y.
        //
        //    Input, int M, the number of rows of the matrix A.
        //    0 <= M.
        //
        //    Input, int N, the number of columns of the matrix A.
        //    0 <= N.
        //
        //    Input, double ALPHA, the scalar multiplier for A * x.
        //
        //    Input, double A[LDA*N].  The M x N subarray contains
        //    the matrix A.
        //
        //    Input, int LDA, the the first dimension of A as declared
        //    in the calling routine.  max ( 1, M ) <= LDA.
        //
        //    Input, double X[*], an array containing the vector to be 
        //    multiplied by the matrix A.  
        //    If TRANS = 'N' or 'n', then X must contain N entries, stored in INCX 
        //    increments in a space of at least ( 1 + ( N - 1 ) * abs ( INCX ) ) 
        //    locations.
        //    Otherwise, X must contain M entries, store in INCX increments
        //    in a space of at least ( 1 + ( M - 1 ) * abs ( INCX ) ) locations.
        //
        //    Input, int INCX, the increment for the elements of
        //    X.  INCX must not be zero.
        //
        //    Input, double BETA, the scalar multiplier for Y.
        //
        //    Input/output, double Y[*], an array containing the vector to
        //    be scaled and incremented by A*X.
        //    If TRANS = 'N' or 'n', then Y must contain M entries, stored in INCY
        //    increments in a space of at least ( 1 + ( M - 1 ) * abs ( INCY ) ) 
        //    locations.
        //    Otherwise, Y must contain N entries, store in INCY increments
        //    in a space of at least ( 1 + ( N - 1 ) * abs ( INCY ) ) locations.
        //
        //    Input, int INCY, the increment for the elements of
        //    Y.  INCY must not be zero.
        //
    {
        int i;
        int iy;
        int j;
        int lenx;
        int leny;
        double temp;
        //
        //  Test the input parameters.
        //
        int info = 0;
        if (trans != 'N' &&
            trans != 'T' &&
            trans != 'C')
        {
            info = 1;
        }
        else
        {
            switch (m)
            {
                case < 0:
                    info = 2;
                    break;
                default:
                {
                    switch (n)
                    {
                        case < 0:
                            info = 3;
                            break;
                        default:
                        {
                            if (lda < Math.Max(1, m))
                            {
                                info = 6;
                            }
                            else
                            {
                                info = incx switch
                                {
                                    0 => 8,
                                    _ => incy switch
                                    {
                                        0 => 11,
                                        _ => info
                                    }
                                };
                            }

                            break;
                        }
                    }

                    break;
                }
            }
        }

        if (info != 0)
        {
            typeMethods.xerbla("DGEMV", info);
            return;
        }

        //
        //  Quick return if possible.
        //
        if (m == 0 ||
            n == 0 ||
            alpha == 0.0 && Math.Abs(beta - 1.0) <= double.Epsilon)
        {
            return;
        }

        switch (trans)
        {
            //
            //  Set LENX and LENY, the lengths of the vectors x and y, and set
            //  up the start points in X and Y.
            //
            case 'N':
                lenx = n;
                leny = m;
                break;
            default:
                lenx = m;
                leny = n;
                break;
        }

        int kx = incx switch
        {
            > 0 => 0,
            _ => 0 - (lenx - 1) * incx
        };

        int ky = incy switch
        {
            > 0 => 0,
            _ => 0 - (leny - 1) * incy
        };

        //
        //  Start the operations. In this version the elements of A are
        //  accessed sequentially with one pass through A.
        //
        //  First form  y := beta*y.
        //
        if (Math.Abs(beta - 1.0) > double.Epsilon)
        {
            switch (incy)
            {
                case 1 when beta == 0.0:
                {
                    for (i = 0; i < leny; i++)
                    {
                        y[i] = 0.0;
                    }

                    break;
                }
                case 1:
                {
                    for (i = 0; i < leny; i++)
                    {
                        y[i] = beta * y[i];
                    }

                    break;
                }
                default:
                {
                    iy = ky;
                    switch (beta)
                    {
                        case 0.0:
                        {
                            for (i = 0; i < leny; i++)
                            {
                                y[iy] = 0.0;
                                iy += incy;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < leny; i++)
                            {
                                y[iy] = beta * y[iy];
                                iy += incy;
                            }

                            break;
                        }
                    }

                    break;
                }
            }
        }

        switch (alpha)
        {
            case 0.0:
                return;
        }

        switch (trans)
        {
            //
            //  Form y := alpha*A*x + y.
            //
            case 'N':
            {
                int jx = kx;
                switch (incy)
                {
                    case 1:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (x[jx] != 0.0)
                            {
                                temp = alpha * x[jx];
                                for (i = 0; i < m; i++)
                                {
                                    y[i] += temp * a[i + j * lda];
                                }
                            }

                            jx += incx;
                        }

                        break;
                    }
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (x[jx] != 0.0)
                            {
                                temp = alpha * x[jx];
                                iy = ky;
                                for (i = 0; i < m; i++)
                                {
                                    y[iy] += temp * a[i + j * lda];
                                    iy += incy;
                                }
                            }

                            jx += incx;
                        }

                        break;
                    }
                }

                break;
            }
            /*
        Form y := alpha*A'*x + y.
        */
            default:
            {
                int jy = ky;
                switch (incx)
                {
                    case 1:
                    {
                        for (j = 0; j < n; j++)
                        {
                            temp = 0.0;
                            for (i = 0; i < m; i++)
                            {
                                temp += a[i + j * lda] * x[i];
                            }

                            y[jy] += alpha * temp;
                            jy += incy;
                        }

                        break;
                    }
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            temp = 0.0;
                            int ix = kx;
                            for (i = 0; i < m; i++)
                            {
                                temp += a[i + j * lda] * x[ix];
                                ix += incx;
                            }

                            y[jy] += alpha * temp;
                            jy += incy;
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void dger(int m, int n, double alpha, double[] x, int incx, double[] y,
            int incy, ref double[] a, int lda )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGER computes A := alpha*x*y' + A.
        //
        //  Discussion:
        //
        //    DGER performs the rank 1 operation
        //
        //      A := alpha*x*y' + A,
        //
        //    where alpha is a scalar, x is an m element vector, y is an n element
        //    vector and A is an m by n matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 April 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra,  Jeremy Du Croz,  
        //    Sven Hammarling,  Richard Hanson.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix A.
        //    0 <= M.
        //
        //    Input, int N, the number of columns of the matrix A.
        //    0 <= N.
        //
        //    Input, double ALPHA, the scalar multiplier.
        //
        //    Input, double X[1+(M-1)*abs(INCX)], the first vector.
        //
        //    Input, int INCX, the increment for elements of X.
        //    INCX must not be zero.
        //
        //    Input, double Y[1+(N-1)*abs(INCY)], the second vector.
        //
        //    Input, int INCY, the increment for elements of Y.
        //    INCY must not be zero.
        //
        //    Input/output, double A[LDA*N].  On entry, the leading M by N 
        //    part of the array contains the matrix of coefficients. On exit, A is
        //    overwritten by the updated matrix.
        //
        //    Input, int LDA, the first dimension of A as declared
        //    in the calling program. max ( 1, M ) <= LDA.
        //
    {
        int i;
        int j;
        double temp;
        //
        //  Test the input parameters.
        //
        int info = 0;
        switch (m)
        {
            case < 0:
                info = 1;
                break;
            default:
            {
                switch (n)
                {
                    case < 0:
                        info = 2;
                        break;
                    default:
                    {
                        switch (incx)
                        {
                            case 0:
                                info = 5;
                                break;
                            default:
                            {
                                switch (incy)
                                {
                                    case 0:
                                        info = 7;
                                        break;
                                    default:
                                    {
                                        if (lda < Math.Max(1, m))
                                        {
                                            info = 9;
                                        }

                                        break;
                                    }
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

        if (info != 0)
        {
            typeMethods.xerbla("DGER", info);
            return;
        }

        //
        //  Quick return if possible.
        //
        if (m == 0 || n == 0 || alpha == 0.0)
        {
            return;
        }

        int jy = incy switch
        {
            //
            //  Start the operations. In this version the elements of A are
            //  accessed sequentially with one pass through A.
            //
            > 0 => 0,
            _ => 0 - (n - 1) * incy
        };

        switch (incx)
        {
            case 1:
            {
                for (j = 0; j < n; j++)
                {
                    if (y[jy] != 0.0)
                    {
                        temp = alpha * y[jy];
                        for (i = 0; i < m; i++)
                        {
                            a[i + j * lda] += x[i] * temp;
                        }
                    }

                    jy += incy;
                }

                break;
            }
            default:
            {
                int kx = incx switch
                {
                    > 0 => 0,
                    _ => 0 - (m - 1) * incx
                };

                for (j = 0; j < n; j++)
                {
                    if (y[jy] != 0.0)
                    {
                        temp = alpha * y[jy];
                        int ix = kx;
                        for (i = 0; i < m; i++)
                        {
                            a[i + j * lda] += x[ix] * temp;
                            ix += incx;
                        }
                    }

                    jy += incy;
                }

                break;
            }
        }
    }

    public static void dtrmv(char uplo, char trans, char diag, int n, double[] a, int lda,
            ref double[] x, int incx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRMV computes x: = A*x or x = A'*x for a triangular matrix A.
        //
        //  Discussion:
        //
        //    DTRMV performs one of the matrix-vector operations
        //
        //      x := A*x,   or   x := A'*x,
        //
        //    where x is an n element vector and  A is an n by n unit, or non-unit,
        //    upper or lower triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 April 2014
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, char UPLO, specifies whether the matrix is an upper or
        //    lower triangular matrix as follows:
        //    'u' or 'U': A is an upper triangular matrix.
        //    'l' or 'L': A is a lower triangular matrix.
        //
        //    Input, char TRANS, specifies the operation to be performed as
        //    follows:
        //    'n' or 'N': x := A*x.
        //    't' or 'T': x := A'*x.
        //    'c' or 'C': x := A'*x.
        //
        //    Input, char DIAG, specifies whether or not A is unit
        //    triangular as follows:
        //    'u' or 'U': A is assumed to be unit triangular.
        //    'n' or 'N': A is not assumed to be unit triangular.
        //
        //    Input, int N, the order of the matrix A.
        //    0 <= N.
        //
        //    Input, double A[LDA*N].
        //    Before entry with  UPLO = 'u' or 'U', the leading n by n
        //    upper triangular part of the array A must contain the upper
        //    triangular matrix and the strictly lower triangular part of
        //    A is not referenced.
        //    Before entry with UPLO = 'l' or 'L', the leading n by n
        //    lower triangular part of the array A must contain the lower
        //    triangular matrix and the strictly upper triangular part of
        //    A is not referenced.
        //    Note that when  DIAG = 'u' or 'U', the diagonal elements of
        //    A are not referenced either, but are assumed to be unity.
        //
        //    Input, int LDA, the first dimension of A as declared
        //    in the calling program. max ( 1, N ) <= LDA.
        //
        //    Input/output, double X[1+(N-1)*abs( INCX)].
        //    Before entry, the incremented array X must contain the n
        //    element vector x. On exit, X is overwritten with the
        //    tranformed vector x.
        //
        //    Input, int INCX, the increment for the elements of
        //    X.  INCX must not be zero.
        //
    {
        int ix;
        int j;
        int jx;
        int kx = 0;
        double temp;
        //
        //  Test the input parameters.
        //
        int info = 0;
        if (uplo != 'U' && uplo != 'L')
        {
            info = 1;
        }
        else if (trans != 'N' && trans != 'T' &&
                 trans != 'C')
        {
            info = 2;
        }
        else if (diag != 'U' && diag != 'N')
        {
            info = 3;
        }
        else
        {
            switch (n)
            {
                case < 0:
                    info = 4;
                    break;
                default:
                {
                    if (lda < Math.Max(1, n))
                    {
                        info = 6;
                    }
                    else
                    {
                        info = incx switch
                        {
                            0 => 8,
                            _ => info
                        };
                    }

                    break;
                }
            }
        }

        if (info != 0)
        {
            typeMethods.xerbla("DTRMV", info);
            return;
        }

        switch (n)
        {
            //
            //  Quick return if possible.
            //
            case 0:
                return;
        }

        bool nounit = diag == 'N';
        switch (incx)
        {
            //
            //  Set up the start point in X if the increment is not unity. This
            //  will be  ( N - 1 ) * INCX  too small for descending loops.
            //
            case <= 0:
                kx = 0 - (n - 1) * incx;
                break;
            default:
            {
                if (incx != 1)
                {
                    kx = 0;
                }

                break;
            }
        }

        switch (trans)
        {
            //
            //  Start the operations. In this version the elements of A are
            //  accessed sequentially with one pass through A.
            //
            //
            //  Form x := A*x.
            //
            case 'N' when uplo == 'U':
            {
                switch (incx)
                {
                    case 1:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (x[j] != 0.0)
                            {
                                temp = x[j];
                                for (int i = 0; i < j; i++)
                                {
                                    x[i] += temp * a[i + j * lda];
                                }

                                switch (nounit)
                                {
                                    case true:
                                        x[j] *= a[j + j * lda];
                                        break;
                                }
                            }
                        }

                        break;
                    }
                    default:
                    {
                        jx = kx;
                        for (j = 0; j < n; j++)
                        {
                            if (x[jx] != 0.0)
                            {
                                temp = x[jx];
                                ix = kx;
                                for (int i = 0; i < j; i++)
                                {
                                    x[ix] += temp * a[i + j * lda];
                                    ix += incx;
                                }

                                switch (nounit)
                                {
                                    case true:
                                        x[jx] *= a[j + j * lda];
                                        break;
                                }
                            }

                            jx += incx;
                        }

                        break;
                    }
                }

                break;
            }
            case 'N' when incx == 1:
            {
                for (j = n - 1; 0 <= j; j--)
                {
                    if (x[j] != 0.0)
                    {
                        temp = x[j];
                        for (int i = n - 1; j < i; i--)
                        {
                            x[i] += temp * a[i + j * lda];
                        }

                        switch (nounit)
                        {
                            case true:
                                x[j] *= a[j + j * lda];
                                break;
                        }
                    }
                }

                break;
            }
            case 'N':
            {
                kx += (n - 1) * incx;
                jx = kx;
                for (j = n - 1; 0 <= j; j--)
                {
                    if (x[jx] != 0.0)
                    {
                        temp = x[jx];
                        ix = kx;
                        for (int i = n - 1; j < i; i--)
                        {
                            x[ix] += temp * a[i + j * lda];
                            ix -= incx;
                        }

                        switch (nounit)
                        {
                            case true:
                                x[jx] *= a[j + j * lda];
                                break;
                        }
                    }

                    jx -= incx;
                }

                break;
            }
            //
            default:
            {
                switch (uplo)
                {
                    case 'U' when incx == 1:
                    {
                        for (j = n - 1; 0 <= j; j--)
                        {
                            temp = x[j];
                            switch (nounit)
                            {
                                case true:
                                    temp *= a[j + j * lda];
                                    break;
                            }

                            for (int i = j - 1; 0 <= i; i--)
                            {
                                temp += a[i + j * lda] * x[i];
                            }

                            x[j] = temp;
                        }

                        break;
                    }
                    case 'U':
                    {
                        jx = kx + (n - 1) * incx;
                        for (j = n - 1; 0 <= j; j--)
                        {
                            temp = x[jx];
                            ix = jx;
                            switch (nounit)
                            {
                                case true:
                                    temp *= a[j + j * lda];
                                    break;
                            }

                            for (int i = j - 1; 0 <= i; i--)
                            {
                                ix -= incx;
                                temp += a[i + j * lda] * x[ix];
                            }

                            x[jx] = temp;
                            jx -= incx;
                        }

                        break;
                    }
                    default:
                    {
                        switch (incx)
                        {
                            case 1:
                            {
                                for (j = 0; j < n; j++)
                                {
                                    temp = x[j];
                                    switch (nounit)
                                    {
                                        case true:
                                            temp *= a[j + j * lda];
                                            break;
                                    }

                                    for (int i = j + 1; i < n; i++)
                                    {
                                        temp += a[i + j * lda] * x[i];
                                    }

                                    x[j] = temp;
                                }

                                break;
                            }
                            default:
                            {
                                jx = kx;
                                for (j = 0; j < n; j++)
                                {
                                    temp = x[jx];
                                    ix = jx;
                                    switch (nounit)
                                    {
                                        case true:
                                            temp *= a[j + j * lda];
                                            break;
                                    }

                                    for (int i = j + 1; i < n; i++)
                                    {
                                        ix += incx;
                                        temp += a[i + j * lda] * x[ix];
                                    }

                                    x[jx] = temp;
                                    jx += incx;
                                }

                                break;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }
}