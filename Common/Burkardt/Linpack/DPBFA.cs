using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DPBFA
{
    public static int dpbfa(ref double[] abd, int lda, int n, int m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPBFA factors a symmetric positive definite matrix stored in band form.
        //
        //  Discussion:
        //
        //    DPBFA is usually called by DPBCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //
        //    If A is a symmetric positive definite band matrix,
        //    the following program segment will set up the input.
        //
        //      m = (band width above diagonal)
        //      do j = 1, n
        //        i1 = max ( 1, j-m )
        //        do i = i1, j
        //          k = i-j+m+1
        //          abd(k,j) = a(i,j)
        //        end do
        //      end do
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
        //    Input/output, double ABD[LDA*N].  On input, the matrix to be
        //    factored.  The columns of the upper triangle are stored in the columns
        //    of ABD and the diagonals of the upper triangle are stored in the
        //    rows of ABD.  On output, an upper triangular matrix R, stored in band
        //    form, so that A = R' * R.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //    M+1 <= LDA is required.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //
        //    Output, int DPBFA, error indicator.
        //    0, for normal return.
        //    K, if the leading minor of order K is not positive definite.
        //
    {
        int ik;
        int info;
        int j;
        int jk;
        int k;
        int mu;
        double s;
        double t;

        for (j = 1; j <= n; j++)
        {
            s = 0.0;
            ik = m + 1;
            jk = Math.Max(j - m, 1);
            mu = Math.Max(m + 2 - j, 1);

            for (k = mu; k <= m; k++)
            {
                t = abd[k - 1 + (j - 1) * lda]
                    - BLAS1D.ddot(k - mu, abd, 1, abd, 1, xIndex:  + ik - 1 + (jk - 1) * lda, yIndex: + mu - 1 + (j - 1) * lda);
                t /= abd[m + (jk - 1) * lda];
                abd[k - 1 + (j - 1) * lda] = t;
                s += t * t;
                ik -= 1;
                jk += 1;
            }

            s = abd[m + (j - 1) * lda] - s;

            switch (s)
            {
                case <= 0.0:
                    info = j;
                    return info;
                default:
                    abd[m + (j - 1) * lda] = Math.Sqrt(s);
                    break;
            }
        }

        info = 0;

        return info;
    }
}