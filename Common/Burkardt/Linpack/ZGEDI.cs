using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZGEDI
{
    public static void zgedi(ref Complex[] a, int lda, int n, int[] ipvt,
            ref Complex[] det, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGEDI computes the determinant and inverse of a matrix.
        //
        //  Discussion:
        //
        //    The matrix must have been factored by ZGECO or ZGEFA.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal and the inverse is requested.
        //    It will not occur if the subroutines are called correctly
        //    and if ZGECO has set 0.0 < RCOND or ZGEFA has set
        //    INFO == 0.
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
        //    Input/output, Complex A[LDA*N]; on input, the factor information
        //    from ZGECO or ZGEFA.  On output, the inverse matrix, if it
        //    was requested,
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int IPVT[N], the pivot vector from ZGECO or ZGEFA.
        //
        //    Output, Complex DET[2], the determinant of the original matrix,
        //    if requested.  Otherwise not referenced.
        //    Determinant = DET(1) * 10.0**DET(2) with
        //    1.0 <= typeMethods.zabs1 ( DET(1) ) < 10.0 or DET(1) == 0.0.
        //    Also, DET(2) is strictly real.
        //
        //    Input, int JOB.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
    {
        int i;
        int j;
        int k;
        int l;
        Complex t;
        Complex[] work;
        //
        //  Compute the determinant.
        //
        if (job / 10 != 0)
        {
            det[0] = new Complex(1.0, 0.0);
            det[1] = new Complex(0.0, 0.0);

            for (i = 1; i <= n; i++)
            {
                if (ipvt[i - 1] != i)
                {
                    det[0] = -det[0];
                }

                det[0] = a[i - 1 + (i - 1) * lda] * det[0];

                if (typeMethods.zabs1(det[0]) == 0.0)
                {
                    break;
                }

                while (typeMethods.zabs1(det[0]) < 1.0)
                {
                    det[0] *= new Complex(10.0, 0.0);
                    det[1] -= new Complex(1.0, 0.0);
                }

                while (10.0 <= typeMethods.zabs1(det[0]))
                {
                    det[0] /= new Complex(10.0, 0.0);
                    det[1] += new Complex(1.0, 0.0);
                }
            }
        }

        //
        //  Compute inverse(U).
        //
        if (job % 10 != 0)
        {
            work = new Complex[n];

            for (k = 1; k <= n; k++)
            {
                a[k - 1 + (k - 1) * lda] = new Complex(1.0, 0.0) / a[k - 1 + (k - 1) * lda];
                t = -a[k - 1 + (k - 1) * lda];
                BLAS1Z.zscal(k - 1, t, ref a, 1, index: +0 + (k - 1) * lda);

                for (j = k + 1; j <= n; j++)
                {
                    t = a[k - 1 + (j - 1) * lda];
                    a[k - 1 + (j - 1) * lda] = new Complex(0.0, 0.0);
                    BLAS1Z.zaxpy(k, t, a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda, yIndex: +0 + (j - 1) * lda);
                }
            }

            //
            //  Form inverse(U) * inverse(L).
            //
            for (k = n - 1; 1 <= k; k--)
            {
                for (i = k + 1; i <= n; i++)
                {
                    work[i - 1] = a[i - 1 + (k - 1) * lda];
                    a[i - 1 + (k - 1) * lda] = new Complex(0.0, 0.0);
                }

                for (j = k + 1; j <= n; j++)
                {
                    t = work[j - 1];
                    BLAS1Z.zaxpy(n, t, a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda, yIndex: +0 + (k - 1) * lda);
                }

                l = ipvt[k - 1];

                if (l != k)
                {
                    BLAS1Z.zswap(n, ref a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda, yIndex: +0 + (l - 1) * lda);
                }
            }
        }
    }

}