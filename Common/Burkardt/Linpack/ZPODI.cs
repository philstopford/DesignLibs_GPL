using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPODI
{
    public static void zpodi(ref Complex[] a, int lda, int n, ref double[] det, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPODI: determinant, inverse of a complex hermitian positive definite matrix.
        //
        //  Discussion:
        //
        //    The matrix is assumed to have been factored by ZPOCO, ZPOFA or ZQRDC.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal and the inverse is requested.
        //    It will not occur if the subroutines are called correctly
        //    and if ZPOCO or ZPOFA has set INFO == 0.
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
        //    Input/output, Complex A[LDA*N]; on input, the output A from ZPOCO or
        //    ZPOFA, or the output X from ZQRDC.  On output, if ZPOCO or ZPOFA was
        //    used to factor A, then ZPODI produces the upper half of inverse(A).
        //    If ZQRDC was used to decompose X, then ZPODI produces the upper half
        //    of inverse(hermitian(X)*X) where hermitian(X) is the conjugate transpose.
        //    Elements of A below the diagonal are unchanged.
        //    If the units digit of JOB is zero, A is unchanged.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double DET[2], if requested, the determinant of A or of
        //    hermitian(X)*X.  Determinant = DET(1) * 10.0**DET(2) with
        //    1.0 <= abs ( DET(1) ) < 10.0 or DET(1) = 0.0.
        //
        //    Input, int JOB.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
    {
        //
        //  Compute determinant
        //
        if (job / 10 != 0)
        {
            det[0] = 1.0;
            det[1] = 0.0;

            int i;
            for (i = 0; i < n; i++)
            {
                det[0] = det[0] * a[i + i * lda].Real * a[i + i * lda].Real;

                if (det[0] == 0.0)
                {
                    break;
                }

                while (det[0] < 1.0)
                {
                    det[0] *= 10.0;
                    det[1] -= 1.0;
                }

                while (10.0 <= det[0])
                {
                    det[0] /= 10.0;
                    det[1] += 1.0;
                }
            }
        }

        //
        //  Compute inverse(R).
        //
        if (job % 10 != 0)
        {
            int j;
            int k;
            Complex t;
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
            //  Form inverse(R) * hermitian(inverse(R)).
            //
            for (j = 1; j <= n; j++)
            {
                for (k = 1; k <= j - 1; k++)
                {
                    t = Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                    BLAS1Z.zaxpy(k, t, a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda, yIndex: +0 + (k - 1) * lda);
                }

                t = Complex.Conjugate(a[j - 1 + (j - 1) * lda]);
                BLAS1Z.zscal(j, t, ref a, 1, index: +0 + (j - 1) * lda);
            }
        }
    }

}