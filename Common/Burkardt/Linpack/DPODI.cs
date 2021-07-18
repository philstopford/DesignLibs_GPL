using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DPODI
    {
        public static void dpodi(ref double[] a, int lda, int n, ref double[] det, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPODI computes the determinant and inverse of a certain matrix.
        //
        //  Discussion:
        //
        //    The matrix is real symmetric positive definite.
        //    DPODI uses the factors computed by DPOCO, DPOFA or DQRDC.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal and the inverse is requested.
        //    It will not occur if the subroutines are called correctly
        //    and if DPOCO or DPOFA has set INFO == 0.
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
        //    Input/output, double A[LDA*N].  On input, the output A from
        //    DPOCO or DPOFA, or the output X from DQRDC.  On output, if DPOCO or
        //    DPOFA was used to factor A then DPODI produces the upper half of
        //    inverse(A).  If DQRDC was used to decompose X then DPODI produces
        //    the upper half of inverse(X'*X) where X' is the transpose.
        //    Elements of A below the diagonal are unchanged.  If the units digit
        //    of JOB is zero, A is unchanged.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, int JOB, specifies the task.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
        //    Output, double DET[2], the determinant of A or of X'*X
        //    if requested.
        //      determinant = DET[0] * 10.0**DET[1]
        //    with 1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
        //
        {
            int i;
            int j;
            int k;
            double s;
            double t;
            //
            //  Compute the determinant.
            //
            if (job / 10 != 0)
            {
                det[0] = 1.0;
                det[1] = 0.0;
                s = 10.0;

                for (i = 1; i <= n; i++)
                {
                    det[0] = det[0] * a[i - 1 + (i - 1) * lda] * a[i - 1 + (i - 1) * lda];

                    if (det[0] == 0.0)
                    {
                        break;
                    }

                    while (det[0] < 1.0)
                    {
                        det[0] = det[0] * s;
                        det[1] = det[1] - 1.0;
                    }

                    while (s <= det[0])
                    {
                        det[0] = det[0] / s;
                        det[1] = det[1] + 1.0;
                    }
                }
            }

            //
            //  Compute inverse(R).
            //
            if ((job % 10) != 0)
            {
                for (k = 1; k <= n; k++)
                {
                    a[k - 1 + (k - 1) * lda] = 1.0 / a[k - 1 + (k - 1) * lda];
                    t = -a[k - 1 + (k - 1) * lda];
                    BLAS1D.dscal(k - 1, t, ref a, 1, index: + 0 + (k - 1) * lda);

                    for (j = k + 1; j <= n; j++)
                    {
                        t = a[k - 1 + (j - 1) * lda];
                        a[k - 1 + (j - 1) * lda] = 0.0;
                        BLAS1D.daxpy(k, t, a, 1, ref a, 1, xIndex: + 0 + (k - 1) * lda, yIndex: + 0 + (j - 1) * lda);
                    }
                }

                //
                //  Form inverse(R) * (inverse(R))'.
                //
                for (j = 1; j <= n; j++)
                {
                    for (k = 1; k <= j - 1; k++)
                    {
                        t = a[k - 1 + (j - 1) * lda];
                        BLAS1D.daxpy(k, t, a, 1, ref a, 1, xIndex:  + 0 + (j - 1) * lda, yIndex:  + 0 + (k - 1) * lda);
                    }

                    t = a[j - 1 + (j - 1) * lda];
                    BLAS1D.dscal(j, t, ref a, 1, index:  + 0 + (j - 1) * lda);
                }
            }
        }
    }
}