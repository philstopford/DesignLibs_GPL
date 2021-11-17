using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DPPFA
{
    public static int dppfa(ref double[] ap, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPPFA factors a real symmetric positive definite matrix in packed form.
        //
        //  Discussion:
        //
        //    DPPFA is usually called by DPPCO, but it can be called
        //    directly with a saving in time if RCOND is not needed.
        //
        //  Packed storage:
        //
        //    The following program segment will pack the upper
        //    triangle of a symmetric matrix.
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
        //    Input/output, double AP[N*(N+1)/2].  On input, the packed
        //    form of a symmetric matrix A.  The columns of the upper triangle are
        //    stored sequentially in a one-dimensional array.  On output, an upper
        //    triangular matrix R, stored in packed form, so that A = R'*R.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int DPPFA, error flag.
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
        double t;

        info = 0;
        jj = 0;

        for (j = 1; j <= n; j++)
        {
            s = 0.0;
            kj = jj;
            kk = 0;

            for (k = 1; k <= j - 1; k++)
            {
                kj += 1;
                t = ap[kj - 1] - BLAS1D.ddot(k - 1, ap, 1, ap, 1, xIndex: + kk, yIndex: + jj);
                kk += k;
                t /= ap[kk - 1];
                ap[kj - 1] = t;
                s += t * t;
            }

            jj += j;
            s = ap[jj - 1] - s;

            switch (s)
            {
                case <= 0.0:
                    info = j;
                    return info;
                default:
                    ap[jj - 1] = Math.Sqrt(s);
                    break;
            }
        }

        return info;
    }

}