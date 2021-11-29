using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPPDI
{
    public static void zppdi(ref Complex[] ap, int n, ref double[] det, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPPDI: determinant, inverse of a complex hermitian positive definite matrix.
        //
        //  Discussion:
        //
        //    The matrix is assumed to have been factored by ZPPCO or ZPPFA.
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
        //    Input/output, complex <double> A[(N*(N+1))/2]; on input, the output from ZPPCO
        //    or ZPPFA.  On output, the upper triangular half of the inverse.
        //    The strict lower triangle is unaltered.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double DET[2], the determinant of original matrix if requested.
        //    Otherwise not referenced.  Determinant = DET(1) * 10.0**DET(2)
        //    with 1.0 <= DET(1) < 10.0 or DET(1) == 0.0.
        //
        //    Input, int JOB.
        //    11, both determinant and inverse.
        //    01, inverse only.
        //    10, determinant only.
        //
    {
        //
        //  Compute determinant.
        //
        if (job / 10 != 0)
        {
            det[0] = 1.0;
            det[1] = 0.0;
            int ii = 0;

            int i;
            for (i = 1; i <= n; i++)
            {
                ii += i;
                det[0] = det[0] * ap[ii - 1].Real * ap[ii - 1].Real;

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
        //  Compute inverse ( R ).
        //
        if (job % 10 == 0)
        {
            return;
        }

        int kk = 0;

        int k1;
        int kj;
        int k;
        Complex t;
        int j;
        int j1;
        for (k = 1; k <= n; k++)
        {
            k1 = kk + 1;
            kk += k;
            ap[kk - 1] = new Complex(1.0, 0.0) / ap[kk - 1];
            t = -ap[kk - 1];
            BLAS1Z.zscal(k - 1, t, ref ap, 1, index: + k1 - 1);
            int kp1 = k + 1;
            j1 = kk + 1;
            kj = kk + k;

            for (j = kp1; j <= n; j++)
            {
                t = ap[kj - 1];
                ap[kj - 1] = new Complex(0.0, 0.0);
                BLAS1Z.zaxpy(k, t, ap, 1, ref ap, 1, xIndex: + k1 - 1, yIndex: + j1 - 1);
                j1 += j;
                kj += j;
            }
        }

        //
        //  Form inverse ( R ) * hermitian ( inverse ( R ) ).
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
                t = Complex.Conjugate(ap[kj - 1]);
                BLAS1Z.zaxpy(k, t, ap, 1, ref ap, 1, xIndex: + j1 - 1, yIndex: k1 - 1);
                k1 += k;
                kj += 1;
            }

            t = Complex.Conjugate(ap[jj - 1]);
            BLAS1Z.zscal(j, t, ref ap, 1, index: + j1 - 1);
        }
    }

}