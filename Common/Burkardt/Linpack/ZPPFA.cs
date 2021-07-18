using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZPPFA
    {
        public static int zppfa(ref Complex[] ap, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZPPFA factors a complex hermitian positive definite packed matrix.
            //
            //  Discussion:
            //
            //    The following program segment will pack the upper triangle of a
            //    hermitian matrix.
            //
            //      k = 0
            //      do j = 1, n
            //        do i = 1, j
            //          k = k + 1
            //          ap(k) = a(i,j)
            //        end do
            //      end do
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
            //    Input/output, Complex AP[N*(N+1)/2]; on input, the packed form
            //    of a hermitian matrix A.  The columns of the upper triangle are
            //    stored sequentially in a one-dimensional array.  On output, an
            //    upper triangular matrix R, stored in packed form, so that
            //      A = hermitian(R) * R.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int ZPPFA.
            //    0, for normal return.
            //    K, if the leading minor of order K is not positive definite.
            //
        {
            int info;
            int j;
            int jj;
            int k;
            int kj;
            int kk;
            double s;
            Complex t;

            info = 0;
            jj = 0;

            for (j = 1; j <= n; j++)
            {
                s = 0.0;
                kj = jj;
                kk = 0;

                for (k = 1; k <= j - 1; k++)
                {
                    kj = kj + 1;
                    t = ap[kj - 1] - BLAS1Z.zdotc(k - 1, ap, 1, ap, 1, xIndex: +kk, yIndex: +jj);
                    kk = kk + k;
                    t = t / ap[kk - 1];
                    ap[kj - 1] = t;
                    s = s + (t * Complex.Conjugate(t)).Real;
                }

                jj = jj + j;
                s = (ap[jj - 1].Real) - s;

                if (s <= 0.0 || (ap[jj - 1].Imaginary) != 0.0)
                {
                    info = j;
                    break;
                }

                ap[jj - 1] = new Complex(Math.Sqrt(s), 0.0);
            }

            return info;
        }

    }
}