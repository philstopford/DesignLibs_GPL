using System;
using Burkardt.Types;

namespace Burkardt.BLAS;

public static class BLAS3D
{
    public static void dgemm(char transa, char transb, int m, int n, int k,
            double alpha, double[] a, int lda, double[] b, int ldb, double beta,
            ref double[] c, int ldc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEMM computes C = alpha * A * B and related operations.
        //
        //  Discussion:
        //
        //    DGEMM performs one of the matrix-matrix operations
        //
        //     C := alpha * op ( A ) * op ( B ) + beta * C,
        //
        //    where op ( X ) is one of
        //
        //      op ( X ) = X   or   op ( X ) = X',
        //
        //    ALPHA and BETA are scalars, and A, B and C are matrices, with op ( A )
        //    an M by K matrix, op ( B ) a K by N matrix and C an N by N matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 February 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, char TRANSA, specifies the form of op( A ) to be used in
        //    the matrix multiplication as follows:
        //    'N' or 'n', op ( A ) = A.
        //    'T' or 't', op ( A ) = A'.
        //    'C' or 'c', op ( A ) = A'.
        //
        //    Input, char TRANSB, specifies the form of op ( B ) to be used in
        //    the matrix multiplication as follows:
        //    'N' or 'n', op ( B ) = B.
        //    'T' or 't', op ( B ) = B'.
        //    'C' or 'c', op ( B ) = B'.
        //
        //    Input, int M, the number of rows of the  matrix op ( A ) and of the  
        //    matrix C.  0 <= M.
        //
        //    Input, int N, the number  of columns of the matrix op ( B ) and the 
        //    number of columns of the matrix C.  0 <= N.
        //
        //    Input, int K, the number of columns of the matrix op ( A ) and the 
        //    number of rows of the matrix op ( B ).  0 <= K.
        //
        //    Input, double ALPHA, the scalar multiplier 
        //    for op ( A ) * op ( B ).
        //
        //    Input, double A(LDA,KA), where:
        //    if TRANSA is 'N' or 'n', KA is equal to K, and the leading M by K
        //    part of the array contains A;
        //    if TRANSA is not 'N' or 'n', then KA is equal to M, and the leading
        //    K by M part of the array must contain the matrix A.
        //
        //    Input, int LDA, the first dimension of A as declared in the calling 
        //    routine.  When TRANSA = 'N' or 'n' then LDA must be at least max ( 1, M ), 
        //    otherwise LDA must be at least max ( 1, K ).
        //
        //    Input, double B(LDB,KB), where:
        //    if TRANSB is 'N' or 'n', kB is N, and the leading K by N 
        //    part of the array contains B;
        //    if TRANSB is not 'N' or 'n', then KB is equal to K, and the leading
        //    N by K part of the array must contain the matrix B.
        //
        //    Input, int LDB, the first dimension of B as declared in the calling 
        //    routine.  When TRANSB = 'N' or 'n' then LDB must be at least max ( 1, K ), 
        //    otherwise LDB must be at least max ( 1, N ).
        //
        //    Input, double BETA, the scalar multiplier for C.
        //
        //    Input/output, double C[LDC*N].
        //    On input, the leading M by N part of this array must contain the 
        //    matrix C, except when BETA is 0.0, in which case C need not be set.
        //    On output, the array C is overwritten by the M by N matrix
        //    alpha * op ( A ) * op ( B ) + beta * C.
        //
        //    Input, int LDC, the first dimension of C as declared in the calling 
        //    routine.  max ( 1, M ) <= LDC.
        //
    {
        int i;
        //int info;
        int j;
        int l;
        //int ncola;
        int nrowa;
        int nrowb;
        bool nota;
        bool notb;
        double temp;
        //
        //  Set NOTA and NOTB as true if A and B respectively are not
        //  transposed and set NROWA, NCOLA and NROWB as the number of rows
        //  and columns of A and the number of rows of B respectively.
        //
        nota = transa == 'N' || transa == 'n';

        nrowa = nota switch
        {
            true => m,
            _ => k
        };

        notb = transb == 'N' || transb == 'n';

        nrowb = notb switch
        {
            true => k,
            _ => n
        };
        //
        //  Test the input parameters.
        //
        //info = 0;

        switch ((transa == 'N' || transa == 'n' ||
                 transa == 'C' || transa == 'c' ||
                 transa == 'T' || transa == 't'))
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("DGEMM - Fatal error!");
                Console.WriteLine("  Input TRANSA has an illegal value.");
                return;
        }

        switch ((transb == 'N' || transb == 'n' ||
                 transb == 'C' || transb == 'c' ||
                 transb == 'T' || transb == 't'))
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("DGEMM - Fatal error!");
                Console.WriteLine("  Input TRANSB has an illegal value.");
                return;
        }

        switch (m)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("DGEMM - Fatal error!");
                Console.WriteLine("  Input M has an illegal value.");
                return;
        }

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("DGEMM - Fatal error!");
                Console.WriteLine("  Input N has an illegal value.");
                return;
        }

        switch (k)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("DGEMM - Fatal error!");
                Console.WriteLine("  Input K has an illegal value.");
                return;
        }

        if (lda < Math.Max(1, nrowa))
        {
            Console.WriteLine("");
            Console.WriteLine("DGEMM - Fatal error!");
            Console.WriteLine("  Input LDA has an illegal value.");
            return;
        }

        if (ldb < Math.Max(1, nrowb))
        {
            Console.WriteLine("");
            Console.WriteLine("DGEMM - Fatal error!");
            Console.WriteLine("  Input LDB has an illegal value.");
            return;
        }

        if (ldc < Math.Max(1, m))
        {
            Console.WriteLine("");
            Console.WriteLine("DGEMM - Fatal error!");
            Console.WriteLine("  Input LDC has an illegal value.");
            return;
        }

        switch (m)
        {
            //
            //  Quick return if possible.
            //
            case 0:
                return;
        }

        switch (n)
        {
            case 0:
                return;
        }

        switch (alpha)
        {
            case 0.0 or 0 when Math.Abs(beta - 1.0) <= double.Epsilon:
                return;
            //
            //  And if alpha is 0.0.
            //
            case 0.0:
            {
                switch (beta)
                {
                    case 0.0:
                    {
                        for (j = 0; j < n; j++)
                        {
                            for (i = 0; i < m; i++)
                            {
                                c[i + j * ldc] = 0.0;
                            }
                        }

                        break;
                    }
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            for (i = 0; i < m; i++)
                            {
                                c[i + j * ldc] = beta * c[i + j * ldc];
                            }
                        }

                        break;
                    }
                }

                return;
            }
        }

        switch (notb)
        {
            //
            //  Start the operations.
            //
            //
            //  Form  C := alpha*A*B + beta*C.
            //
            case true when nota:
            {
                for (j = 0; j < n; j++)
                {
                    switch (beta)
                    {
                        case 0.0:
                        {
                            for (i = 0; i < m; i++)
                            {
                                c[i + j * ldc] = 0.0;
                            }

                            break;
                        }
                        default:
                        {
                            if (Math.Abs(beta - 1.0) > double.Epsilon)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    c[i + j * ldc] = beta * c[i + j * ldc];
                                }
                            }

                            break;
                        }
                    }

                    for (l = 0; l < k; l++)
                    {
                        if (b[l + j * ldb] != 0.0)
                        {
                            temp = alpha * b[l + j * ldb];
                            for (i = 0; i < m; i++)
                            {
                                c[i + j * ldc] += temp * a[i + l * lda];
                            }
                        }
                    }
                }

                break;
            }
            //
            //  Form  C := alpha*A'*B + beta*C
            //
            case true:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        temp = 0.0;
                        for (l = 0; l < k; l++)
                        {
                            temp += a[l + i * lda] * b[l + j * ldb];
                        }

                        c[i + j * ldc] = beta switch
                        {
                            0.0 => alpha * temp,
                            _ => alpha * temp + beta * c[i + j * ldc]
                        };
                    }
                }

                break;
            }
            //
            default:
            {
                switch (nota)
                {
                    case true:
                    {
                        for (j = 0; j < n; j++)
                        {
                            switch (beta)
                            {
                                case 0.0:
                                {
                                    for (i = 0; i < m; i++)
                                    {
                                        c[i + j * ldc] = 0.0;
                                    }

                                    break;
                                }
                                default:
                                {
                                    if (Math.Abs(beta - 1.0) > double.Epsilon)
                                    {
                                        for (i = 0; i < m; i++)
                                        {
                                            c[i + j * ldc] = beta * c[i + j * ldc];
                                        }
                                    }

                                    break;
                                }
                            }

                            for (l = 0; l < k; l++)
                            {
                                if (b[j + l * ldb] != 0.0)
                                {
                                    temp = alpha * b[j + l * ldb];
                                    for (i = 0; i < m; i++)
                                    {
                                        c[i + j * ldc] += temp * a[i + l * lda];
                                    }
                                }
                            }
                        }

                        break;
                    }
                    //
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            for (i = 0; i < m; i++)
                            {
                                temp = 0.0;
                                for (l = 0; l < k; l++)
                                {
                                    temp += a[l + i * lda] * b[j + l * ldb];
                                }

                                c[i + j * ldc] = beta switch
                                {
                                    0.0 => alpha * temp,
                                    _ => alpha * temp + beta * c[i + j * ldc]
                                };
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }

    public static void dtrmm(char side, char uplo, char transa, char diag, int m, int n,
            double alpha, double[] a, int lda, ref double[] b, int ldb )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRMM performs B:=A*B or B:=B*A, A triangular, B rectangular.
        //
        //  Discussion:
        //
        //    This routine performs one of the matrix-matrix operations
        //      B := alpha*op( A )*B,
        //    or
        //      B := alpha*B*op( A ),
        //    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
        //    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
        //      op( A ) = A   or   op( A ) = A'.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 April 2014
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, char SIDE, specifies whether op(A) multiplies B from
        //    the left or right as follows:
        //    'L' or 'l': B := alpha*op( A )*B.
        //    'R' or 'r': B := alpha*B*op( A ).
        //
        //    Input, char UPLO, specifies whether the matrix A is an upper or
        //    lower triangular matrix as follows:
        //    'U' or 'u': A is an upper triangular matrix.
        //    'L' or 'l': A is a lower triangular matrix.
        //
        //    Input, char TRANS, specifies the form of op( A ) to be used in
        //    the matrix multiplication as follows:
        //    'N' or 'n': op( A ) = A.
        //    'T' or 't': op( A ) = A'.
        //    'C' or 'c': op( A ) = A'.
        //
        //    Input, char DIAG, specifies whether or not A is unit triangular
        //    as follows:
        //    'U' or 'u': A is assumed to be unit triangular.
        //    'N' or 'n': A is not assumed to be unit triangular.
        //
        //    Input, int M, the number of rows of B.  0 <= M.
        //
        //    Input, int N, the number of columns of B.  
        //    0 <= N.
        //
        //    Input, double ALPHA, the scalar  alpha.  When alpha is
        //    0.0, A is not referenced and B need not be set before entry.
        //
        //    Input, double A[LDA*K], where k is m when  SIDE = 'L' or 'l'  
        //    and is  n  when  SIDE = 'R' or 'r'.
        //    Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
        //    upper triangular part of the array  A must contain the upper
        //    triangular matrix  and the strictly lower triangular part of
        //    A is not referenced.
        //    Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
        //    lower triangular part of the array  A must contain the lower
        //    triangular matrix  and the strictly upper triangular part of
        //    A is not referenced.
        //    Note that when  DIAG = 'U' or 'u',  the diagonal elements of
        //    A  are not referenced either,  but are assumed to be  unity.
        //
        //    Input, integer LDA, the first dimension of A as declared
        //    in the calling program.  When SIDE = 'L' or 'l' then LDA must be at 
        //    least max ( 1, M ); when SIDE = 'R' or 'r', LDA must be at least 
        //    max ( 1, N ).
        //
        //    Input/output, double B[LDB*N].
        //    Before entry, the leading m by n part of the array  B must contain 
        //    the matrix  B, and on exit is overwritten by the transformed matrix.
        //
        //    Input, integer LDB, the first dimension of B as declared
        //    in the calling program.   max ( 1, M ) <= LDB.
        //
    {
        int i;
        int info;
        int j;
        int k;
        bool lside;
        bool nounit;
        int nrowa;
        double temp;
        bool upper;
        //
        //  Test the input parameters.
        //
        lside = side == 'L';

        nrowa = lside switch
        {
            true => m,
            _ => n
        };

        nounit = diag == 'N';
        upper = uplo == 'U';

        info = 0;
        switch (lside)
        {
            case false when side != 'R':
                info = 1;
                break;
            default:
            {
                switch (upper)
                {
                    case false when uplo != 'L':
                        info = 2;
                        break;
                    default:
                    {
                        if (transa != 'N' && transa != 'T' &&
                            transa != 'C')
                        {
                            info = 3;
                        }
                        else if (diag != 'U' && diag != 'N')
                        {
                            info = 4;
                        }
                        else
                        {
                            switch (m)
                            {
                                case < 0:
                                    info = 5;
                                    break;
                                default:
                                {
                                    switch (n)
                                    {
                                        case < 0:
                                            info = 6;
                                            break;
                                        default:
                                        {
                                            if (lda < Math.Max(1, nrowa))
                                            {
                                                info = 9;
                                            }
                                            else if (ldb < Math.Max(1, m))
                                            {
                                                info = 11;
                                            }

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

                break;
            }
        }

        if (info != 0)
        {
            typeMethods.xerbla("DTRMM", info);
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

        switch (alpha)
        {
            //
            //  And when alpha is 0.0.
            //
            case 0.0:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        b[i + j * ldb] = 0.0;
                    }
                }

                return;
            }
        }

        switch (lside)
        {
            //
            //  Start the operations.
            //
            //
            //  Form  B := alpha*A*B.
            //
            case true when transa == 'N':
            {
                switch (upper)
                {
                    case true:
                    {
                        for (j = 0; j < n; j++)
                        {
                            for (k = 0; k < m; k++)
                            {
                                if (b[k + j * ldb] != 0.0)
                                {
                                    temp = alpha * b[k + j * ldb];
                                    for (i = 0; i < k; i++)
                                    {
                                        b[i + j * ldb] += temp * a[i + k * lda];
                                    }

                                    switch (nounit)
                                    {
                                        case true:
                                            temp *= a[k + k * lda];
                                            break;
                                    }

                                    b[k + j * ldb] = temp;
                                }
                            }
                        }

                        break;
                    }
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            for (k = m - 1; 0 <= k; k--)
                            {
                                if (b[k + j * ldb] != 0.0)
                                {
                                    temp = alpha * b[k + j * ldb];
                                    b[k + j * ldb] = temp;
                                    switch (nounit)
                                    {
                                        case true:
                                            b[k + j * ldb] *= a[k + k * lda];
                                            break;
                                    }

                                    for (i = k + 1; i < m; i++)
                                    {
                                        b[i + j * ldb] += temp * a[i + k * lda];
                                    }
                                }
                            }
                        }

                        break;
                    }
                }

                break;
            }
            //
            //  Form  B := alpha*A'*B.
            //
            case true when upper:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = m - 1; 0 <= i; i--)
                    {
                        temp = b[i + j * ldb];
                        switch (nounit)
                        {
                            case true:
                                temp *= a[i + i * lda];
                                break;
                        }

                        for (k = 0; k < i; k++)
                        {
                            temp += a[k + i * lda] * b[k + j * ldb];
                        }

                        b[i + j * ldb] = alpha * temp;
                    }
                }

                break;
            }
            case true:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        temp = b[i + j * ldb];
                        switch (nounit)
                        {
                            case true:
                                temp *= a[i + i * lda];
                                break;
                        }

                        for (k = i + 1; k < m; k++)
                        {
                            temp += a[k + i * lda] * b[k + j * ldb];
                        }

                        b[i + j * ldb] = alpha * temp;
                    }
                }

                break;
            }
            //
            default:
            {
                switch (transa)
                {
                    case 'n' when upper:
                    {
                        for (j = n - 1; 0 <= j; j--)
                        {
                            temp = alpha;
                            switch (nounit)
                            {
                                case true:
                                    temp *= a[j + j * lda];
                                    break;
                            }

                            for (i = 0; i < m; i++)
                            {
                                b[i + j * ldb] = temp * b[i + j * ldb];
                            }

                            for (k = 0; k < j; k++)
                            {
                                if (a[k + j * lda] != 0.0)
                                {
                                    temp = alpha * a[k + j * lda];
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] += temp * b[i + k * ldb];
                                    }
                                }
                            }
                        }

                        break;
                    }
                    case 'n':
                    {
                        for (j = 0; j < n; j++)
                        {
                            temp = alpha;
                            switch (nounit)
                            {
                                case true:
                                    temp *= a[j + j * lda];
                                    break;
                            }

                            for (i = 0; i < m; i++)
                            {
                                b[i + j * ldb] = temp * b[i + j * ldb];
                            }

                            for (k = j + 1; k < n; k++)
                            {
                                if (a[k + j * lda] != 0.0)
                                {
                                    temp = alpha * a[k + j * lda];
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] += temp * b[i + k * ldb];
                                    }
                                }
                            }
                        }

                        break;
                    }
                    //
                    default:
                    {
                        switch (upper)
                        {
                            case true:
                            {
                                for (k = 0; k < n; k++)
                                {
                                    for (j = 0; j < k; j++)
                                    {
                                        if (a[j + k * lda] != 0.0)
                                        {
                                            temp = alpha * a[j + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + j * ldb] += temp * b[i + k * ldb];
                                            }
                                        }
                                    }

                                    temp = alpha;
                                    switch (nounit)
                                    {
                                        case true:
                                            temp *= a[k + k * lda];
                                            break;
                                    }

                                    if (Math.Abs(temp - 1.0) > double.Epsilon)
                                    {
                                        for (i = 0; i < m; i++)
                                        {
                                            b[i + k * ldb] = temp * b[i + k * ldb];
                                        }
                                    }
                                }

                                break;
                            }
                            default:
                            {
                                for (k = n - 1; 0 <= k; k--)
                                {
                                    for (j = k + 1; j < n; j++)
                                    {
                                        if (a[j + k * lda] != 0.0)
                                        {
                                            temp = alpha * a[j + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + j * ldb] += temp * b[i + k * ldb];
                                            }
                                        }
                                    }

                                    temp = alpha;
                                    switch (nounit)
                                    {
                                        case true:
                                            temp *= a[k + k * lda];
                                            break;
                                    }

                                    if (Math.Abs(temp - 1.0) > double.Epsilon)
                                    {
                                        for (i = 0; i < m; i++)
                                        {
                                            b[i + k * ldb] = temp * b[i + k * ldb];
                                        }
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
    }

    public static void dtrsm(char side, char uplo, char transa, char diag, int m, int n,
            double alpha, double[] a, int lda, ref double[] b, int ldb )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRSM solves A*X=alpha*B or X*A=alpha*B, for triangular A, rectangular B.
        //
        //  Discussion:
        //
        //    DTRSM solves one of the matrix equations
        //      op( A )*X = alpha*B,   
        //    or
        //      X*op( A ) = alpha*B,
        //    where alpha is a scalar, X and B are m by n matrices, A is a unit, or
        //    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
        //      op( A ) = A   or   op( A ) = A'.
        //    The matrix X is overwritten on B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, char SIDE, specifies whether op( A ) appears on the left
        //    or right of X as follows:
        //    'L' or 'l': op( A )*X = alpha*B.
        //    'R' or 'r': X*op( A ) = alpha*B.
        //
        //    Input, char UPLO, specifies whether the matrix A is an upper or
        //    lower triangular matrix as follows:
        //    'U' or 'u': A is an upper triangular matrix.
        //    'L' or 'l': A is a lower triangular matrix.
        //
        //    Input, char TRANSA, specifies the form of op( A ) to be used in
        //    the matrix multiplication as follows:
        //    'N' or 'n': op( A ) = A.
        //    'T' or 't': op( A ) = A'.
        //    'C' or 'c': op( A ) = A'.
        //
        //    Input, char DIAG, specifies whether or not A is unit triangular
        //    as follows:
        //    'U' or 'u': A is assumed to be unit triangular.
        //    'N' or 'n': A is not assumed to be unit triangular.
        //
        //    Input, int M, the number of rows of B.  0 <= M.
        //
        //    Input, int N, the number of columns of B.  0 <= N.
        //
        //    Input, double ALPHA, the scalar alpha.  When alpha is
        //    0.0 then A is not referenced and B need not be set before entry.
        //
        //    Input, double A[LDA*K] where K is M when SIDE = 'L' or 'l'  
        //    and K is N when SIDE = 'R' or 'r'.
        //    Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
        //    upper triangular part of the array  A must contain the upper
        //    triangular matrix  and the strictly lower triangular part of
        //    A is not referenced.
        //    Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
        //    lower triangular part of the array  A must contain the lower
        //    triangular matrix  and the strictly upper triangular part of
        //    A is not referenced.
        //    Note that when  DIAG = 'U' or 'u',  the diagonal elements of
        //    A  are not referenced either,  but are assumed to be  unity.
        //
        //    Input, int LDA, the first dimension of A as declared
        //    in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
        //    LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
        //    then LDA must be at least max( 1, n ).
        //
        //    Input/output, double B[LDB*N].
        //    Before entry, the leading m by n part of the array B must
        //    contain the right-hand side matrix B, and on exit is
        //    overwritten by the solution matrix X.
        //
        //    Input, int LDB, the first dimension of B as declared
        //    in the calling program.  LDB must be at least max ( 1, M ).
        //
    {
        int i;
        int info;
        int j;
        int k;
        bool lside;
        bool nounit;
        int nrowa;
        double temp;
        bool upper;
        //
        //  Test the input parameters.
        //
        lside = side == 'L';

        nrowa = lside switch
        {
            true => m,
            _ => n
        };

        nounit = diag == 'N';
        upper = uplo == 'U';

        info = 0;

        switch (lside)
        {
            case false when side != 'R':
                info = 1;
                break;
            default:
            {
                switch (upper)
                {
                    case false when uplo != 'L':
                        info = 2;
                        break;
                    default:
                    {
                        if (transa != 'N' &&
                            transa != 'T' &&
                            transa != 'C')
                        {
                            info = 3;
                        }
                        else if (diag != 'U' && diag != 'N')
                        {
                            info = 4;
                        }
                        else
                        {
                            switch (m)
                            {
                                case < 0:
                                    info = 5;
                                    break;
                                default:
                                {
                                    switch (n)
                                    {
                                        case < 0:
                                            info = 6;
                                            break;
                                        default:
                                        {
                                            if (lda < Math.Max(1, nrowa))
                                            {
                                                info = 9;
                                            }
                                            else if (ldb < Math.Max(1, m))
                                            {
                                                info = 11;
                                            }

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

                break;
            }
        }

        if (info != 0)
        {
            typeMethods.xerbla("DTRSM", info);
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

        switch (alpha)
        {
            //
            //  and when alpha is 0.0.
            //
            case 0.0:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        b[i + j * ldb] = 0.0;
                    }
                }

                return;
            }
        }

        switch (lside)
        {
            //
            //  Start the operations.
            //
            //
            //  Form  B := alpha*inv( a )*B.
            //
            case true when transa == 'N':
            {
                switch (upper)
                {
                    case true:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (Math.Abs(alpha - 1.0) > double.Epsilon)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    b[i + j * ldb] = alpha * b[i + j * ldb];
                                }
                            }

                            for (k = m - 1; 0 <= k; k--)
                            {
                                if (b[k + j * ldb] != 0.0)
                                {
                                    switch (nounit)
                                    {
                                        case true:
                                            b[k + j * ldb] /= a[k + k * lda];
                                            break;
                                    }

                                    for (i = 0; i < k; i++)
                                    {
                                        b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda];
                                    }
                                }
                            }
                        }

                        break;
                    }
                    default:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (Math.Abs(alpha - 1.0) > double.Epsilon)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    b[i + j * ldb] = alpha * b[i + j * ldb];
                                }
                            }

                            for (k = 0; k < m; k++)
                            {
                                if (b[k + j * ldb] != 0.0)
                                {
                                    switch (nounit)
                                    {
                                        case true:
                                            b[k + j * ldb] /= a[k + k * lda];
                                            break;
                                    }

                                    for (i = k + 1; i < m; i++)
                                    {
                                        b[i + j * ldb] -= b[k + j * ldb] * a[i + k * lda];
                                    }
                                }
                            }
                        }

                        break;
                    }
                }

                break;
            }
            //
            //  Form  B := alpha*inv( A' )*B.
            //
            case true when upper:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        temp = alpha * b[i + j * ldb];
                        for (k = 0; k < i; k++)
                        {
                            temp -= a[k + i * lda] * b[k + j * ldb];
                        }

                        switch (nounit)
                        {
                            case true:
                                temp /= a[i + i * lda];
                                break;
                        }

                        b[i + j * ldb] = temp;
                    }
                }

                break;
            }
            case true:
            {
                for (j = 0; j < n; j++)
                {
                    for (i = m - 1; 0 <= i; i--)
                    {
                        temp = alpha * b[i + j * ldb];
                        for (k = i + 1; k < m; k++)
                        {
                            temp -= a[k + i * lda] * b[k + j * ldb];
                        }

                        switch (nounit)
                        {
                            case true:
                                temp /= a[i + i * lda];
                                break;
                        }

                        b[i + j * ldb] = temp;
                    }
                }

                break;
            }
            //
            default:
            {
                switch (transa)
                {
                    case 'N' when upper:
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (Math.Abs(alpha - 1.0) > double.Epsilon)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    b[i + j * ldb] = alpha * b[i + j * ldb];
                                }
                            }

                            for (k = 0; k < j; k++)
                            {
                                if (a[k + j * lda] != 0.0)
                                {
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb];
                                    }
                                }
                            }

                            switch (nounit)
                            {
                                case true:
                                {
                                    temp = 1.0 / a[j + j * lda];
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] = temp * b[i + j * ldb];
                                    }

                                    break;
                                }
                            }
                        }

                        break;
                    }
                    case 'N':
                    {
                        for (j = n - 1; 0 <= j; j--)
                        {
                            if (Math.Abs(alpha - 1.0) > double.Epsilon)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    b[i + j * ldb] = alpha * b[i + j * ldb];
                                }
                            }

                            for (k = j + 1; k < n; k++)
                            {
                                if (a[k + j * lda] != 0.0)
                                {
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] -= a[k + j * lda] * b[i + k * ldb];
                                    }
                                }
                            }

                            switch (nounit)
                            {
                                case true:
                                {
                                    temp = 1.0 / a[j + j * lda];
                                    for (i = 0; i < m; i++)
                                    {
                                        b[i + j * ldb] = temp * b[i + j * ldb];
                                    }

                                    break;
                                }
                            }
                        }

                        break;
                    }
                    //
                    default:
                    {
                        switch (upper)
                        {
                            case true:
                            {
                                for (k = n - 1; 0 <= k; k--)
                                {
                                    switch (nounit)
                                    {
                                        case true:
                                        {
                                            temp = 1.0 / a[k + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + k * ldb] = temp * b[i + k * ldb];
                                            }

                                            break;
                                        }
                                    }

                                    for (j = 0; j < k; j++)
                                    {
                                        if (a[j + k * lda] != 0.0)
                                        {
                                            temp = a[j + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + j * ldb] -= temp * b[i + k * ldb];
                                            }
                                        }
                                    }

                                    if (Math.Abs(alpha - 1.0) > double.Epsilon)
                                    {
                                        for (i = 0; i < m; i++)
                                        {
                                            b[i + k * ldb] = alpha * b[i + k * ldb];
                                        }
                                    }
                                }

                                break;
                            }
                            default:
                            {
                                for (k = 0; k < n; k++)
                                {
                                    switch (nounit)
                                    {
                                        case true:
                                        {
                                            temp = 1.0 / a[k + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + k * ldb] = temp * b[i + k * ldb];
                                            }

                                            break;
                                        }
                                    }

                                    for (j = k + 1; j < n; j++)
                                    {
                                        if (a[j + k * lda] != 0.0)
                                        {
                                            temp = a[j + k * lda];
                                            for (i = 0; i < m; i++)
                                            {
                                                b[i + j * ldb] -= temp * b[i + k * ldb];
                                            }
                                        }
                                    }

                                    if (Math.Abs(alpha - 1.0) > double.Epsilon)
                                    {
                                        for (i = 0; i < m; i++)
                                        {
                                            b[i + k * ldb] = alpha * b[i + k * ldb];
                                        }
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
    }

}