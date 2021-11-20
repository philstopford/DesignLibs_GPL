using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPOFA
{
    public static int zpofa(ref Complex[] a, int lda, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPOFA factors a complex hermitian positive definite matrix.
        //
        //  Discussion:
        //
        //    ZPOFA is usually called by ZPOCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //    (time for ZPOCO) = (1 + 18/N) * (time for ZPOFA).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //
        //  Parameters:
        //
        //    Input/output, Complex A[LDA*N]; On input, the hermitian matrix to be
        //    factored.  On output, an upper triangular matrix R so that
        //      A = hermitian(R)*R
        //    where hermitian(R) is the conjugate transpose.  The strict lower
        //    triangle is unaltered.  If INFO /= 0, the factorization is not
        //    complete.  Only the diagonal and upper triangle are used.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int ZPOFA.
        //    0, for normal return.
        //    K, signals an error condition.  The leading minor of order K is
        //    not positive definite.
        //
    {
        int j;

        int info = 0;

        for (j = 1; j <= n; j++)
        {
            double s = 0.0;
            int k;
            for (k = 1; k <= j - 1; k++)
            {
                Complex t = a[k - 1 + (j - 1) * lda] - BLAS1Z.zdotc(k - 1, a, 1, a, 1, xIndex: +0 + (k - 1) * lda,
                    yIndex: +0 + (j - 1) * lda);
                t /= a[k - 1 + (k - 1) * lda];
                a[k - 1 + (j - 1) * lda] = t;
                s += (t * Complex.Conjugate(t)).Real;
            }

            s = a[j - 1 + (j - 1) * lda].Real - s;

            if (s <= 0.0 || a[j - 1 + (j - 1) * lda].Imaginary != 0.0)
            {
                info = j;
                break;
            }

            a[j - 1 + (j - 1) * lda] = new Complex(Math.Sqrt(s), 0.0);
        }

        return info;
    }

}