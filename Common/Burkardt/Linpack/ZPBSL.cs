using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPBSL
{
    public static void zpbsl(Complex[] abd, int lda, int n, int m,
            ref Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPBSL solves a complex hermitian positive definite band system.
        //
        //  Discussion:
        //
        //    The system matrix must have been factored by ZPBCO or ZPBFA.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
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
        //    Input, Complex ABD[LDA*N], the output from ZPBCO or ZPBFA.
        //
        //    Input, int LDA, the leading dimension of ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        int la;
        int lb;
        int lm;
        Complex t;
        //
        //  Solve hermitian(R) * Y = B.
        //
        for (k = 1; k <= n; k++)
        {
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            t = BLAS1Z.zdotc(lm, abd, 1, b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
            b[k - 1] = (b[k - 1] - t) / abd[m + (k - 1) * lda];
        }

        //
        //  Solve R * X = Y.
        //
        for (k = n; 1 <= k; k--)
        {
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            b[k - 1] /= abd[m + (k - 1) * lda];
            t = -b[k - 1];
            BLAS1Z.zaxpy(lm, t, abd, 1, ref b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
        }

    }

}