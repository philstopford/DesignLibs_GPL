using System;

namespace Burkardt.Linpack;

public static class DGTSL
{
    public static int dgtsl(int n, ref double[] c, ref double[] d, ref double[] e, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGTSL solves a general tridiagonal linear system.
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
        //    Input, int N, the order of the tridiagonal matrix.
        //
        //    Input/output, double C[N], contains the subdiagonal of the
        //    tridiagonal matrix in entries C(2:N).  On output, C is destroyed.
        //
        //    Input/output, double D[N].  On input, the diagonal of the
        //    matrix.  On output, D is destroyed.
        //
        //    Input/output, double E[N], contains the superdiagonal of the
        //    tridiagonal matrix in entries E(1:N-1).  On output E is destroyed.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Output, int DGTSL, error flag.
        //    0, normal value.
        //    K, the K-th element of the diagonal becomes exactly zero.  The
        //       subroutine returns if this error condition is detected.
        //
    {
        int k;

        int info = 0;
        c[0] = d[0];

        switch (n)
        {
            case >= 2:
            {
                d[0] = e[0];
                e[0] = 0.0;
                e[n - 1] = 0.0;

                for (k = 1; k <= n - 1; k++)
                {
                    //
                    //  Find the larger of the two rows.
                    //
                    double t;
                    if (Math.Abs(c[k - 1]) <= Math.Abs(c[k]))
                    {
                        //
                        //  Interchange rows.
                        //
                        t = c[k];
                        c[k] = c[k - 1];
                        c[k - 1] = t;

                        t = d[k];
                        d[k] = d[k - 1];
                        d[k - 1] = t;

                        t = e[k];
                        e[k] = e[k - 1];
                        e[k - 1] = t;

                        t = b[k];
                        b[k] = b[k - 1];
                        b[k - 1] = t;
                    }

                    switch (c[k - 1])
                    {
                        //
                        //  Zero elements.
                        //
                        case 0.0:
                            info = k;
                            return info;
                    }

                    t = -c[k] / c[k - 1];
                    c[k] = d[k] + t * d[k - 1];
                    d[k] = e[k] + t * e[k - 1];
                    e[k] = 0.0;
                    b[k] += t * b[k - 1];
                }

                break;
            }
        }

        switch (c[n - 1])
        {
            case 0.0:
                info = n;
                return info;
        }

        //
        //  Back solve.
        //
        b[n - 1] /= c[n - 1];

        switch (n)
        {
            case > 1:
            {
                b[n - 2] = (b[n - 2] - d[n - 2] * b[n - 1]) / c[n - 2];

                for (k = n - 2; 1 <= k; k--)
                {
                    b[k - 1] = (b[k - 1] - d[k - 1] * b[k] - e[k - 1] * b[k + 1]) / c[k - 1];
                }

                break;
            }
        }

        return info;
    }
}