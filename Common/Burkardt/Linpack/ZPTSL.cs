using System.Numerics;

namespace Burkardt.Linpack;

public static class ZPTSL
{
    public static void zptsl(int n, ref Complex[] d, ref Complex[] e,
            ref Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPTSL solves a Hermitian positive definite tridiagonal linear system.
        //
        //  Discussion;
        //
        //    The system does not have to be factored first.
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
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, Complex D[N].  On input, the diagonal of the
        //    matrix.  On output, this has been overwritten by other information.
        //
        //    Input/output, Complex E[N].  On input, the superdiagonal
        //    entries of the matrix in locations E(1:N-1).  On output, this has
        //    been overwritten by other information.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        int kbm1;
        int ke;
        int kf;
        int kp1;
        int nm1d2;
        Complex t1;
        Complex t2;
        switch (n)
        {
            //
            //  Check for 1 x 1 case.
            //
            case 1:
                b[0] /= d[0];
                return;
        }

        nm1d2 = (n - 1) / 2;

        if (n != 2)
        {
            kbm1 = n - 1;
            //
            //  Zero top half of subdiagonal and bottom half of superdiagonal.
            //
            for (k = 1; k <= nm1d2; k++)
            {
                t1 = Complex.Conjugate(e[k - 1]) / d[k - 1];
                d[k] -= t1 * e[k - 1];
                b[k] -= t1 * b[k - 1];
                t2 = e[kbm1 - 1] / d[kbm1];
                d[kbm1 - 1] -= t2 * Complex.Conjugate(e[kbm1 - 1]);
                b[kbm1 - 1] -= t2 * b[kbm1];
                kbm1 -= 1;
            }
        }

        kp1 = nm1d2 + 1;
        switch (n % 2)
        {
            //
            //  Clean up for possible 2 x 2 block at center.
            //
            case 0:
                t1 = Complex.Conjugate(e[kp1 - 1]) / d[kp1 - 1];
                d[kp1] -= t1 * e[kp1 - 1];
                b[kp1] -= t1 * b[kp1 - 1];
                kp1 += 1;
                break;
        }

        //
        //  Back solve starting at the center, going towards the top and bottom.
        //
        b[kp1 - 1] /= d[kp1 - 1];

        if (n != 2)
        {
            k = kp1 - 1;
            ke = kp1 + nm1d2 - 1;

            for (kf = kp1; kf <= ke; kf++)
            {
                b[k - 1] = (b[k - 1] - e[k - 1] * b[k]) / d[k - 1];
                b[kf] = (b[kf] - Complex.Conjugate(e[kf - 1]) * b[kf - 1]) / d[kf];
                k -= 1;
            }
        }

        b[0] = (n % 2) switch
        {
            0 => (b[0] - e[0] * b[1]) / d[0],
            _ => b[0]
        };
    }

}