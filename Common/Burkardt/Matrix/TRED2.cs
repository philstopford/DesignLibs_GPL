using System;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static class TRED2
{
    public static void tred2(int n, double[] a, ref double[] d, ref double[] e, ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
        //
        //  Discussion:
        //
        //    This subroutine reduces a real symmetric matrix to a
        //    symmetric tridiagonal matrix using and accumulating
        //    orthogonal similarity transformations.
        //
        //    A and Z may coincide, in which case a single storage area is used
        //    for the input of A and the output of Z.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
        //    Klema, Moler.
        //    C version by John Burkardt.
        //
        //  Reference:
        //
        //    Martin, Reinsch, Wilkinson,
        //    TRED2,
        //    Numerische Mathematik,
        //    Volume 11, pages 181-195, 1968.
        //
        //    James Wilkinson, Christian Reinsch,
        //    Handbook for Automatic Computation,
        //    Volume II, Linear Algebra, Part 2,
        //    Springer, 1971,
        //    ISBN: 0387054146,
        //    LC: QA251.W67.
        //
        //    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
        //    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
        //    Matrix Eigensystem Routines, EISPACK Guide,
        //    Lecture Notes in Computer Science, Volume 6,
        //    Springer Verlag, 1976,
        //    ISBN13: 978-3540075462,
        //    LC: QA193.M37.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the real symmetric input matrix.  Only the
        //    lower triangle of the matrix need be supplied.
        //
        //    Output, double D[N], the diagonal elements of the tridiagonal
        //    matrix.
        //
        //    Output, double E[N], contains the subdiagonal elements of the
        //    tridiagonal matrix in E(2:N).  E(1) is set to zero.
        //
        //    Output, double Z[N*N], the orthogonal transformation matrix
        //    produced in the reduction.
        //
    {
        double g;
        double h;
        int i;
        int j;
        int k;
        int l;

        for (j = 0; j < n; j++)
        {
            for (i = j; i < n; i++)
            {
                z[i + j * n] = a[i + j * n];
            }
        }

        for (j = 0; j < n; j++)
        {
            d[j] = a[n - 1 + j * n];
        }

        for (i = n - 1; 1 <= i; i--)
        {
            l = i - 1;
            h = 0.0;
            //
            //  Scale row.
            //
            double scale = 0.0;
            for (k = 0; k <= l; k++)
            {
                scale += Math.Abs(d[k]);
            }

            switch (scale)
            {
                case 0.0:
                {
                    e[i] = d[l];

                    for (j = 0; j <= l; j++)
                    {
                        d[j] = z[l + j * n];
                        z[i + j * n] = 0.0;
                        z[j + i * n] = 0.0;
                    }

                    d[i] = 0.0;
                    continue;
                }
            }

            for (k = 0; k <= l; k++)
            {
                d[k] /= scale;
            }

            h = 0.0;
            for (k = 0; k <= l; k++)
            {
                h += d[k] * d[k];
            }

            double f = d[l];
            g = -Math.Sqrt(h) * typeMethods.r8_sign(f);
            e[i] = scale * g;
            h -= f * g;
            d[l] = f - g;
            //
            //  Form A*U.
            //
            for (k = 0; k <= l; k++)
            {
                e[k] = 0.0;
            }

            for (j = 0; j <= l; j++)
            {
                f = d[j];
                z[j + i * n] = f;
                g = e[j] + z[j + j * n] * f;

                for (k = j + 1; k <= l; k++)
                {
                    g += z[k + j * n] * d[k];
                    e[k] += z[k + j * n] * f;
                }

                e[j] = g;
            }

            //
            //  Form P.
            //
            for (k = 0; k <= l; k++)
            {
                e[k] /= h;
            }

            f = 0.0;
            for (k = 0; k <= l; k++)
            {
                f += e[k] * d[k];
            }

            double hh = 0.5 * f / h;
            //
            //  Form Q.
            //
            for (k = 0; k <= l; k++)
            {
                e[k] -= hh * d[k];
            }

            //
            //  Form reduced A.
            //
            for (j = 0; j <= l; j++)
            {
                f = d[j];
                g = e[j];

                for (k = j; k <= l; k++)
                {
                    z[k + j * n] = z[k + j * n] - f * e[k] - g * d[k];
                }

                d[j] = z[l + j * n];
                z[i + j * n] = 0.0;
            }

            d[i] = h;
        }

        //
        //  Accumulation of transformation matrices.
        //
        for (i = 1; i < n; i++)
        {
            l = i - 1;
            z[n - 1 + l * n] = z[l + l * n];
            z[l + l * n] = 1.0;
            h = d[i];

            if (h != 0.0)
            {
                for (k = 0; k <= l; k++)
                {
                    d[k] = z[k + i * n] / h;
                }

                for (j = 0; j <= l; j++)
                {
                    g = 0.0;
                    for (k = 0; k <= l; k++)
                    {
                        g += z[k + i * n] * z[k + j * n];
                    }

                    for (k = 0; k <= l; k++)
                    {
                        z[k + j * n] -= g * d[k];
                    }
                }
            }

            for (k = 0; k <= l; k++)
            {
                z[k + i * n] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            d[j] = z[n - 1 + j * n];
        }

        for (j = 0; j < n - 1; j++)
        {
            z[n - 1 + j * n] = 0.0;
        }

        z[n - 1 + (n - 1) * n] = 1.0;

        e[0] = 0.0;
    }

}