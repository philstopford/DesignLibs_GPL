using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DPPDI
{
    public static void dppdi(ref double[] ap, int n, ref double[] det, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPPDI computes the determinant and inverse of a matrix factored by DPPCO or DPPFA.
        //
        //  Discussion:
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
        //    24 May 2005
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
        //    Input/output, double AP[N*(N+1)/2].  On input, the output from
        //    DPPCO or DPPFA.  On output, the upper triangular half of the
        //    inverse, if requested.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double DET[2], the determinant of the original matrix
        //    if requested.
        //      determinant = DET[0] * 10.0**DET[1]
        //    with  1.0D+00 <= DET[0] < 10.0D+00 or DET[0] == 0.0D+00.
        //
        //    Input, int JOB, job request.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
    {
        //
        //  Compute the determinant.
        //
        if (job / 10 != 0)
        {
            det[0] = 1.0;
            det[1] = 0.0;
            double s = 10.0;
            int ii = 0;

            int i;
            for (i = 1; i <= n; i++)
            {
                ii += i;

                det[0] = det[0] * ap[ii - 1] * ap[ii - 1];

                if (det[0] == 0.0)
                {
                    break;
                }

                while (det[0] < 1.0)
                {
                    det[0] *= s;
                    det[1] -= 1.0;
                }

                while (s <= det[0])
                {
                    det[0] /= s;
                    det[1] += 1.0;
                }
            }
        }

        //
        //  Compute inverse(R).
        //
        if (job % 10 == 0)
        {
            return;
        }

        int kk = 0;

        int k;
        int k1;
        int kj;
        int j1;
        double t;
        int j;
        for (k = 1; k <= n; k++)
        {
            k1 = kk + 1;
            kk += k;
            ap[kk - 1] = 1.0 / ap[kk - 1];
            t = -ap[kk - 1];
            BLAS1D.dscal(k - 1, t, ref ap, 1, index: + k1 - 1);
            j1 = kk + 1;
            kj = kk + k;

            for (j = k + 1; j <= n; j++)
            {
                t = ap[kj - 1];
                ap[kj - 1] = 0.0;
                BLAS1D.daxpy(k, t, ap, 1, ref ap, 1, xIndex: + k1 - 1, yIndex: + j1 - 1);
                j1 += j;
                kj += j;
            }
        }

        //
        //  Form inverse(R) * (inverse(R))'.
        //
        int jj = 0;

        for (j = 1; j <= n; j++)
        {
            j1 = jj + 1;
            jj += j;
            k1 = 1;
            kj = j1;

            for (k = 1; k <= j - 1; k++)
            {
                t = ap[kj - 1];
                BLAS1D.daxpy(k, t, ap, 1, ref ap, 1, xIndex: + j1 - 1, yIndex: + k1 - 1);
                k1 += k;
                kj += 1;
            }

            t = ap[jj - 1];
            BLAS1D.dscal(j, t, ref ap, 1, index: + j1 - 1);
        }
    }

}