using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class ZGEFA
    {
        public static int zgefa(ref Complex[] a, int lda, int n, ref int[] ipvt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZGEFA factors a complex matrix by Gaussian elimination.
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
            //    Input/output, Complex A[LDA*N]; on input, the matrix to be factored.
            //    On output, an upper triangular matrix and the multipliers which were
            //    used to obtain it.  The factorization can be written A = L*U where
            //    L is a product of permutation and unit lower triangular matrices and
            //    U is upper triangular.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int IPVT[N], the pivot indices.
            //
            //    Output, int ZGEFA,
            //    0, normal value.
            //    K, if U(K,K) == 0.0.  This is not an error condition for this
            //    subroutine, but it does indicate that ZGESL or ZGEDI will divide by zero
            //    if called.  Use RCOND in ZGECO for a reliable indication of singularity.
            //
        {
            int info;
            int j;
            int k;
            int l;
            Complex t;
            //
            //  Gaussian elimination with partial pivoting.
            //
            info = 0;

            for (k = 1; k <= n - 1; k++)
            {
                //
                //  Find L = pivot index.
                //
                l = BLAS1Z.izamax(n - k + 1, a, 1, index: +(k - 1) + (k - 1) * lda) + k - 1;
                ipvt[k - 1] = l;
                //
                //  Zero pivot implies this column already triangularized.
                //
                if (typeMethods.zabs1(a[l - 1 + (k - 1) * lda]) == 0.0)
                {
                    info = k;
                    continue;
                }

                //
                //  Interchange if necessary.
                //
                if (l != k)
                {
                    t = a[l - 1 + (k - 1) * lda];
                    a[l - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda];
                    a[k - 1 + (k - 1) * lda] = t;
                }

                //
                //  Compute multipliers
                //
                t = -new Complex(1.0, 0.0) / a[k - 1 + (k - 1) * lda];
                BLAS1Z.zscal(n - k, t, ref a, 1, index: +k + (k - 1) * lda);
                //
                //  Row elimination with column indexing
                //
                for (j = k + 1; j <= n; j++)
                {
                    t = a[l - 1 + (j - 1) * lda];
                    if (l != k)
                    {
                        a[l - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda];
                        a[k - 1 + (j - 1) * lda] = t;
                    }

                    BLAS1Z.zaxpy(n - k, t, a, 1, ref a, 1, xIndex: +k + (k - 1) * lda, yIndex: +k + (j - 1) * lda);
                }

            }

            ipvt[n - 1] = n;

            if (typeMethods.zabs1(a[n - 1 + (n - 1) * lda]) == 0.0)
            {
                info = n;
            }

            return info;
        }

    }
}