namespace Burkardt.Linpack;

public static class DPTSL
{
    public static void dptsl(int n, ref double[] d, double[] e, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPTSL solves a positive definite tridiagonal linear system.
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
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double D[N], on input the diagonal of the
        //    tridiagonal matrix.  On output, D is destroyed.
        //
        //    Input, double E[N], the offdiagonal of the tridiagonal matrix in
        //    entries E(1:N-1).
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        int kbm1;
        int ke;
        int kf;
        int kp1;
        int nm1d2;
        double t1;
        double t2;
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

        switch (n)
        {
            case > 2:
            {
                kbm1 = n - 1;
                //
                //  Zero top half of subdiagonal and bottom half of superdiagonal.
                //
                for (k = 1; k <= nm1d2; k++)
                {
                    t1 = e[k - 1] / d[k - 1];
                    d[k] -= t1 * e[k - 1];
                    b[k] -= t1 * b[k - 1];
                    t2 = e[kbm1 - 1] / d[kbm1];
                    d[kbm1 - 1] -= t2 * e[kbm1 - 1];
                    b[kbm1 - 1] -= t2 * b[kbm1];
                    kbm1 -= 1;
                }

                break;
            }
        }

        kp1 = nm1d2 + 1;
        switch (n % 2)
        {
            //
            //  Clean up for possible 2 x 2 block at center.
            //
            case 0:
                t1 = e[kp1 - 1] / d[kp1 - 1];
                d[kp1] -= t1 * e[kp1 - 1];
                b[kp1] -= t1 * b[kp1 - 1];
                kp1 += 1;
                break;
        }

        //
        //  Back solve starting at the center, going towards the top and bottom.
        //
        b[kp1 - 1] /= d[kp1 - 1];

        switch (n)
        {
            case > 2:
            {
                k = kp1 - 1;
                ke = kp1 + nm1d2 - 1;

                for (kf = kp1; kf <= ke; kf++)
                {
                    b[k - 1] = (b[k - 1] - e[k - 1] * b[k]) / d[k - 1];
                    b[kf] = (b[kf] - e[kf - 1] * b[kf - 1]) / d[kf];
                    k -= 1;
                }

                break;
            }
        }

        b[0] = (n % 2) switch
        {
            0 => (b[0] - e[0] * b[1]) / d[0],
            _ => b[0]
        };
    }

}