using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZGBFA
{
    public static int zgbfa(ref Complex[] abd, int lda, int n, int ml, int mu, ref int[] ipvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGBFA factors a complex band matrix by elimination.
        //
        //  Discussion:
        //
        //    ZGBFA is usually called by ZGBCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //
        //  Band storage:
        //
        //    If A is a band matrix, the following program segment
        //    will set up the input.
        //
        //      ml = (band width below the diagonal)
        //      mu = (band width above the diagonal)
        //      m = ml + mu + 1
        //      do j = 1, n
        //        i1 = max ( 1, j - mu )
        //        i2 = min ( n, j + ml )
        //        do i = i1, i2
        //          k = i - j + m
        //          abd(k,j) = a(i,j)
        //        end do
        //      end do
        //
        //    This uses rows ML+1 through 2*ML+MU+1 of ABD.
        //    In addition, the first ML rows in ABD are used for
        //    elements generated during the triangularization.
        //    The total number of rows needed in ABD is 2*ML+MU+1.
        //    The ML+MU by ML+MU upper left triangle and the
        //    ML by ML lower right triangle are not referenced.
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
        //    Input/output, new Complex ABD[LDA*N], on input, contains the matrix in
        //    band storage.  The columns of the matrix are stored in the columns
        //    of ABD and the diagonals of the matrix are stored in rows ML+1
        //    through 2*ML+MU+1 of ABD.  On output, an upper triangular matrix
        //    in band storage and the multipliers which were used to obtain it.
        //    The factorization can be written A = L*U where L is a product of
        //    permutation and unit lower triangular matrices and U is upper triangular.
        //
        //    Input, int LDA, the leading dimension of ABD.
        //    LDA must be at least 2*ML+MU+1.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the number of diagonals below the main diagonal.
        //    0 <= ML < N.
        //
        //    Input, int MU, the number of diagonals above the main diagonal.
        //    0 <= MU < N.  More efficient if ML <= MU.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, int ZGBFA.
        //    0, normal value.
        //    K, if U(K,K) == 0.0.  This is not an error condition for this
        //    subroutine, but it does indicate that ZGBSL will divide by zero if
        //    called.  Use RCOND in ZGBCO for a reliable indication of singularity.
        //
    {
        int i;
        int i0;
        int info;
        int j;
        int j0;
        int j1;
        int ju;
        int jz;
        int k;
        int l;
        int lm;
        int m;
        int mm;
        Complex t;

        m = ml + mu + 1;
        info = 0;
        //
        //  Zero initial fill-in columns.
        //
        j0 = mu + 2;
        j1 = Math.Min(n, m) - 1;

        for (jz = j0; jz <= j1; jz++)
        {
            i0 = m + 1 - jz;
            for (i = i0; i <= ml; i++)
            {
                abd[i - 1 + (jz - 1) * lda] = new Complex(0.0, 0.0);
            }
        }

        jz = j1;
        ju = 0;
        //
        //  Gaussian elimination with partial pivoting.
        //
        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Zero next fill-in column
            //
            jz += 1;
            if (jz <= n)
            {
                for (i = 1; i <= ml; i++)
                {
                    abd[i - 1 + (jz - 1) * lda] = new Complex(0.0, 0.0);
                }
            }

            //
            //  Find L = pivot index.
            //
            lm = Math.Min(ml, n - k);
            l = BLAS1Z.izamax(lm + 1, abd, 1, index: +m - 1 + (k - 1) * lda) + m - 1;
            ipvt[k - 1] = l + k - m;
            //
            //  Zero pivot implies this column already triangularized.
            //
            if (typeMethods.zabs1(abd[l - 1 + (k - 1) * lda]) == 0.0)
            {
                info = k;
                continue;
            }

            //
            //  Interchange if necessary.
            //
            if (l != m)
            {
                t = abd[l - 1 + (k - 1) * lda];
                abd[l - 1 + (k - 1) * lda] = abd[m - 1 + (k - 1) * lda];
                abd[m - 1 + (k - 1) * lda] = t;
            }

            //
            //  Compute multipliers.
            //
            t = -new Complex(1.0, 0.0) / abd[m - 1 + (k - 1) * lda];
            BLAS1Z.zscal(lm, t, ref abd, 1, index: +m + (k - 1) * lda);
            //
            //  Row elimination with column indexing.
            //
            ju = Math.Min(Math.Max(ju, mu + ipvt[k - 1]), n);
            mm = m;

            for (j = k + 1; j <= ju; j++)
            {
                l -= 1;
                mm -= 1;
                t = abd[l - 1 + (j - 1) * lda];
                if (l != mm)
                {
                    abd[l - 1 + (j - 1) * lda] = abd[mm - 1 + (j - 1) * lda];
                    abd[mm - 1 + (j - 1) * lda] = t;
                }

                BLAS1Z.zaxpy(lm, t, abd, 1, ref abd, 1, xIndex: +m + (k - 1) * lda, yIndex: +mm + (j - 1) * lda);
            }
        }

        ipvt[n - 1] = n;

        if (typeMethods.zabs1(abd[m - 1 + (n - 1) * lda]) == 0.0)
        {
            info = n;
        }

        return info;
    }

}