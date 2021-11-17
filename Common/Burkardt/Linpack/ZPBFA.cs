using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPBFA
{
    public static int zpbfa(ref Complex[] abd, int lda, int n, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPBFA factors a complex hermitian positive definite band matrix.
        //
        //  Discussion:
        //
        //    ZPBFA is usually called by ZPBCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
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
        //    Input/output, Complex ABD[LDA*N]; on input, the matrix to be factored.
        //    The columns of the upper triangle are stored in the columns of ABD
        //    and the diagonals of the upper triangle are stored in the rows of ABD.
        //    On output, an upper triangular matrix R, stored in band form, so that
        //    A = hermitian(R)*R.
        //
        //    Input, int LDA, the leading dimension of ABD.
        //    LDA must be at least M+1.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //    0 <= M < N.
        //
        //    Output, int ZSPFA.
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
        Complex t;

        info = 0;

        for (j = 1; j <= n; j++)
        {
            s = 0.0;
            ik = m + 1;
            jk = Math.Max(j - m, 1);
            mu = Math.Max(m + 2 - j, 1);

            for (k = mu; k <= m; k++)
            {
                t = abd[k - 1 + (j - 1) * lda]
                    - BLAS1Z.zdotc(k - mu, abd, 1, abd, 1, xIndex: +ik - 1 + (jk - 1) * lda,
                        yIndex: +mu - 1 + (j - 1) * lda);
                t /= abd[m + (jk - 1) * lda];
                abd[k - 1 + (j - 1) * lda] = t;
                s += (t * Complex.Conjugate(t)).Real;
                ik -= 1;
                jk += 1;
            }

            s = abd[m + (j - 1) * lda].Real - s;

            if (s <= 0.0 || abd[m + (j - 1) * lda].Imaginary != 0.0)
            {
                info = j;
                break;
            }

            abd[m + (j - 1) * lda] = new Complex(Math.Sqrt(s), 0.0);
        }

        return info;
    }

}