using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DGEDI
    {
        public static void dgedi(ref double[] a, int lda, int n, int[] ipvt, ref double[] det,
        double[] work, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGEDI computes the determinant and inverse of a matrix factored by DGECO or DGEFA.
        //
        //  Discussion:
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal and the inverse is requested.
        //    It will not occur if the subroutines are called correctly
        //    and if DGECO has set 0.0 < RCOND or DGEFA has set INFO == 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2005
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
        //    Input/output, double A[LDA*N], on input, the LU factor information,
        //    as output by DGECO or DGEFA.  On output, the inverse
        //    matrix if requested.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
        //
        //    Workspace, double WORK[N].
        //
        //    Output, double DET[2], the determinant of original matrix if
        //    requested.  The determinant = DET[0] * pow ( 10.0, DET[1] )
        //    with  1.0 <= abs ( DET[0] ) < 10.0 or DET[0] == 0.0.
        //
        //    Input, int JOB, specifies what is to be computed.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
        {
            int i;
            int j;
            int k;
            int l;
            double t;
            //
            //  Compute the determinant.
            //
            if (job / 10 != 0)
            {
                det[0] = 1.0;
                det[1] = 0.0;

                for (i = 1; i <= n; i++)
                {
                    if (ipvt[i - 1] != i)
                    {
                        det[0] = -det[0];
                    }

                    det[0] = det[0] * a[i - 1 + (i - 1) * lda];

                    if (det[0] == 0.0)
                    {
                        break;
                    }

                    while (Math.Abs(det[0]) < 1.0)
                    {
                        det[0] = det[0] * 10.0;
                        det[1] = det[1] - 1.0;
                    }

                    while (10.0 <= Math.Abs(det[0]))
                    {
                        det[0] = det[0] / 10.0;
                        det[1] = det[1] + 1.0;
                    }
                }
            }

            //
            //  Compute inverse(U).
            //
            if ((job % 10) != 0)
            {
                for (k = 1; k <= n; k++)
                {
                    a[k - 1 + (k - 1) * lda] = 1.0 / a[k - 1 + (k - 1) * lda];
                    t = -a[k - 1 + (k - 1) * lda];
                    BLAS1D.dscal(k - 1, t, ref a, 1, index:  + 0 + (k - 1) * lda);

                    for (j = k + 1; j <= n; j++)
                    {
                        t = a[k - 1 + (j - 1) * lda];
                        a[k - 1 + (j - 1) * lda] = 0.0;
                        BLAS1D.daxpy(k, t, a, 1, ref a, 1, xIndex:  + 0 + (k - 1) * lda, yIndex:  + 0 + (j - 1) * lda);
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
                        a[i - 1 + (k - 1) * lda] = 0.0;
                    }

                    for (j = k + 1; j <= n; j++)
                    {
                        t = work[j - 1];
                        BLAS1D.daxpy(n, t, a, 1, ref a, 1, xIndex:  + 0 + (j - 1) * lda, yIndex:  + 0 + (k - 1) * lda);
                    }

                    l = ipvt[k - 1];
                    if (l != k)
                    {
                        BLAS1D.dswap(n, ref a, 1, ref a, 1, xIndex:  + 0 + (k - 1) * lda, yIndex:  + 0 + (l - 1) * lda);
                    }
                }
            }
        }
    }
}