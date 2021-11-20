using System.Numerics;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZGTSL
{
    public static int zgtsl(int n, ref Complex[] c, ref Complex[] d,
            ref Complex[] e, ref Complex[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGTSL solves a complex general tridiagonal system.
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
        //    Input/output, Complex C[N]; on input, the subdiagonal
        //    of the tridiagonal matrix in entries C(2:N).  On output, C has
        //    been overwritten.
        //
        //    Input/output, Complex D[N]; on input, the diagonal of
        //    the tridiagonal matrix.  On output, D has been overwritten.
        //
        //    Input/output, Complex E[N]; on input, the superdiagonal
        //    of the tridiagonal matrix in entries E(1:N-1).  On output, E
        //    has been overwritten.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Output, int ZGTSL.
        //    0, normal value.
        //    K, if the K-th element of the diagonal becomes exactly zero.  The
        //    subroutine returns when this is detected.
        //
    {
        int k;

        int info = 0;
        c[0] = d[0];

        switch (n - 1)
        {
            case >= 1:
            {
                d[0] = e[0];
                e[0] = new Complex(0.0, 0.0);
                e[n - 1] = new Complex(0.0, 0.0);

                for (k = 1; k <= n - 1; k++)
                {
                    Complex t;
                    if (typeMethods.zabs1(c[k - 1]) <= typeMethods.zabs1(c[k]))
                    {
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

                    if (typeMethods.zabs1(c[k - 1]) == 0.0)
                    {
                        info = k;
                        return info;
                    }

                    t = -c[k] / c[k - 1];
                    c[k] = d[k] + t * d[k - 1];
                    d[k] = e[k] + t * e[k - 1];
                    e[k] = new Complex(0.0, 0.0);
                    b[k] += t * b[k - 1];
                }

                break;
            }
        }

        if (typeMethods.zabs1(c[n - 1]) == 0.0)
        {
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