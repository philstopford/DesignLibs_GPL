using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DPOFA
{
    public static int dpofa(ref double[] a, int lda, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPOFA factors a real symmetric positive definite matrix.
        //
        //  Discussion:
        //
        //    DPOFA is usually called by DPOCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
        //    Pete Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input/output, double A[LDA*N].  On input, the symmetric matrix
        //    to be  factored.  Only the diagonal and upper triangle are used.
        //    On output, an upper triangular matrix R so that A = R'*R
        //    where R' is the transpose.  The strict lower triangle is unaltered.
        //    If INFO /= 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int DPOFA, error flag.
        //    0, for normal return.
        //    K, signals an error condition.  The leading minor of order K is not
        //    positive definite.
        //
    {
        int info;
        int j;

        for (j = 1; j <= n; j++)
        {
            double s = 0.0;

            int k;
            for (k = 1; k <= j - 1; k++)
            {
                double t = a[k - 1 + (j - 1) * lda] - BLAS1D.ddot(k - 1, a, 1, a, 1, xIndex: + 0 + (k - 1) * lda, yIndex: + 0 + (j - 1) * lda);
                t /= a[k - 1 + (k - 1) * lda];
                a[k - 1 + (j - 1) * lda] = t;
                s += t * t;
            }

            s = a[j - 1 + (j - 1) * lda] - s;

            switch (s)
            {
                case <= 0.0:
                    info = j;
                    return info;
                default:
                    a[j - 1 + (j - 1) * lda] = Math.Sqrt(s);
                    break;
            }
        }

        info = 0;

        return info;
    }
}