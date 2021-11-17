using System;
using Burkardt.Types;

namespace Burkardt.Praxis;

public static class MINFIT
{
    public static void minfit(int n, double tol, ref double[] a, ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MINFIT computes the singular value decomposition of an N by N array.
        //
        //  Discussion:
        //
        //    This is an improved version of the EISPACK routine MINFIT
        //    restricted to the case M = N and P = 0.
        //
        //    The singular values of the array A are returned in Q.  A is
        //    overwritten with the orthogonal matrix V such that U * diag(Q) = A * V,
        //    where U is another orthogonal matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 August 2016
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Richard Brent.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Brent,
        //    Algorithms for Minimization with Derivatives,
        //    Prentice Hall, 1973,
        //    Reprinted by Dover, 2002.
        //
        //    James Wilkinson, Christian Reinsch,
        //    Handbook for Automatic Computation,
        //    Volume II, Linear Algebra, Part 2,
        //    Springer Verlag, 1971.
        //
        //    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, Yasuhiko Ikebe, 
        //    Virginia Klema, Cleve Moler,
        //    Matrix Eigensystem Routines, EISPACK Guide,
        //    Lecture Notes in Computer Science, Volume 6,
        //    Springer Verlag, 1976,
        //    ISBN13: 978-3540075462,
        //    LC: QA193.M37.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, double TOL, a tolerance which determines when a vector
        //    (a column or part of a column of the matrix) may be considered
        //    "essentially" equal to zero.
        //
        //    Input/output, double A[N,N].  On input, an N by N array whose
        //    singular value decomposition is desired.  On output, the
        //    SVD orthogonal matrix factor V.
        //
        //    Input/output, double Q[N], the singular values.
        //
    {
        double c;
        double[] e;
        double eps;
        double f = 0;
        double g;
        double h;
        int i;
        int ii;
        int j;
        int jj;
        int k;
        int kt;
        const int kt_max = 30;
        int l;
        int l2;
        double s;
        bool skip;
        double temp;
        double x;
        double y;
        double z;
        switch (n)
        {
            //
            //  Householder's reduction to bidiagonal form.
            //
            case 1:
                q[0] = a[0 + 0 * n];
                a[0 + 0 * n] = 1.0;
                return;
        }

        e = new double[n];

        eps = typeMethods.r8_epsilon();
        g = 0.0;
        x = 0.0;

        for (i = 1; i <= n; i++)
        {
            e[i - 1] = g;
            l = i + 1;

            s = 0.0;
            for (ii = i; ii <= n; ii++)
            {
                s += a[ii - 1 + (i - 1) * n] * a[ii - 1 + (i - 1) * n];
            }

            g = 0.0;

            if (tol <= s)
            {
                f = a[i - 1 + (i - 1) * n];

                g = f switch
                {
                    >= 0.0 => -g,
                    _ => Math.Sqrt(s)
                };

                h = f * g - s;
                a[i - 1 + (i - 1) * n] = f - g;

                for (j = l; j <= n; j++)
                {
                    f = 0.0;
                    for (ii = i; ii <= n; ii++)
                    {
                        f += a[ii - 1 + (i - 1) * n] * a[ii - 1 + (j - 1) * n];
                    }

                    f /= h;

                    for (ii = i; ii <= n; ii++)
                    {
                        a[ii - 1 + (j - 1) * n] += f * a[ii - 1 + (i - 1) * n];
                    }
                }
            }

            q[i - 1] = g;

            s = 0.0;
            for (j = l; j <= n; j++)
            {
                s += a[i - 1 + (j - 1) * n] * a[i - 1 + (j - 1) * n];
            }

            g = 0.0;

            if (tol <= s)
            {
                if (i < n)
                {
                    f = a[i - 1 + i * n];
                }

                g = f switch
                {
                    >= 0.0 => -g,
                    _ => Math.Sqrt(s)
                };

                h = f * g - s;

                if (i < n)
                {
                    a[i - 1 + i * n] = f - g;
                    for (jj = l; jj <= n; jj++)
                    {
                        e[jj - 1] = a[i - 1 + (jj - 1) * n] / h;
                    }

                    for (j = l; j <= n; j++)
                    {
                        s = 0.0;
                        for (jj = l; jj <= n; jj++)
                        {
                            s += a[j - 1 + (jj - 1) * n] * a[i - 1 + (jj - 1) * n];
                        }

                        for (jj = l; jj <= n; jj++)
                        {
                            a[j - 1 + (jj - 1) * n] += s * e[jj - 1];
                        }
                    }
                }
            }

            y = Math.Abs(q[i - 1]) + Math.Abs(e[i - 1]);

            x = Math.Max(x, y);
        }

        //
        //  Accumulation of right-hand transformations.
        //
        a[n - 1 + (n - 1) * n] = 1.0;
        g = e[n - 1];
        l = n;

        for (i = n - 1; 1 <= i; i--)
        {
            if (g != 0.0)
            {
                h = a[i - 1 + i * n] * g;

                for (ii = l; ii <= n; ii++)
                {
                    a[ii - 1 + (i - 1) * n] = a[i - 1 + (ii - 1) * n] / h;
                }

                for (j = l; j <= n; j++)
                {
                    s = 0.0;
                    for (jj = l; jj <= n; jj++)
                    {
                        s += a[i - 1 + (jj - 1) * n] * a[jj - 1 + (j - 1) * n];
                    }

                    for (ii = l; ii <= n; ii++)
                    {
                        a[ii - 1 + (j - 1) * n] += s * a[ii - 1 + (i - 1) * n];
                    }
                }
            }

            for (jj = l; jj <= n; jj++)
            {
                a[i - 1 + (jj - 1) * n] = 0.0;
            }

            for (ii = l; ii <= n; ii++)
            {
                a[ii - 1 + (i - 1) * n] = 0.0;
            }

            a[i - 1 + (i - 1) * n] = 1.0;

            g = e[i - 1];

            l = i;
        }

        //
        //  Diagonalization of the bidiagonal form.
        //
        eps *= x;

        for (k = n; 1 <= k; k--)
        {
            kt = 0;

            for (;;)
            {
                kt += 1;

                switch (kt)
                {
                    case > kt_max:
                        e[k - 1] = 0.0;
                        Console.WriteLine("");
                        Console.WriteLine("MINFIT - Fatal error!");
                        Console.WriteLine("  The QR algorithm failed to converge.");
                        return;
                }

                skip = false;

                for (l2 = k; 1 <= l2; l2--)
                {
                    l = l2;

                    if (Math.Abs(e[l - 1]) <= eps)
                    {
                        skip = true;
                        break;
                    }

                    if (1 < l)
                    {
                        if (Math.Abs(q[l - 2]) <= eps)
                        {
                            break;
                        }
                    }
                }

                switch (skip)
                {
                    //
                    //  Cancellation of E(L) if 1 < L.
                    //
                    case false:
                    {
                        c = 0.0;
                        s = 1.0;

                        for (i = l; i <= k; i++)
                        {
                            f = s * e[i - 1];
                            e[i - 1] = c * e[i - 1];
                            if (Math.Abs(f) <= eps)
                            {
                                break;
                            }

                            g = q[i - 1];
                            //
                            //  q(i) = h = sqrt(g*g + f*f).
                            //
                            h = typeMethods.r8_hypot(f, g);

                            q[i - 1] = h;

                            switch (h)
                            {
                                case 0.0:
                                    g = 1.0;
                                    h = 1.0;
                                    break;
                            }

                            c = g / h;
                            s = -f / h;
                        }

                        break;
                    }
                }

                //
                //  Test for convergence for this index K.
                //
                z = q[k - 1];

                if (l == k)
                {
                    switch (z)
                    {
                        case < 0.0:
                        {
                            q[k - 1] = -z;
                            for (i = 1; i <= n; i++)
                            {
                                a[i - 1 + (k - 1) * n] = -a[i - 1 + (k - 1) * n];
                            }

                            break;
                        }
                    }

                    break;
                }

                //
                //  Shift from bottom 2*2 minor.
                //
                x = q[l - 1];
                y = q[k - 2];
                g = e[k - 2];
                h = e[k - 1];
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);

                g = typeMethods.r8_hypot(f, 1.0);

                temp = f switch
                {
                    < 0.0 => f - g,
                    _ => f + g
                };

                f = ((x - z) * (x + z) + h * (y / temp - h)) / x;
                //
                //  Next QR transformation.
                //
                c = 1.0;
                s = 1.0;

                for (i = l + 1; i <= k; i++)
                {
                    g = e[i - 1];
                    y = q[i - 1];
                    h = s * g;
                    g *= c;

                    z = typeMethods.r8_hypot(f, h);

                    e[i - 2] = z;

                    switch (z)
                    {
                        case 0.0:
                            f = 1.0;
                            z = 1.0;
                            break;
                    }

                    c = f / z;
                    s = h / z;
                    f = x * c + g * s;
                    g = -x * s + g * c;
                    h = y * s;
                    y *= c;

                    for (j = 1; j <= n; j++)
                    {
                        x = a[j - 1 + (i - 2) * n];
                        z = a[j - 1 + (i - 1) * n];
                        a[j - 1 + (i - 2) * n] = x * c + z * s;
                        a[j - 1 + (i - 1) * n] = -x * s + z * c;
                    }

                    z = typeMethods.r8_hypot(f, h);

                    q[i - 2] = z;

                    switch (z)
                    {
                        case 0.0:
                            f = 1.0;
                            z = 1.0;
                            break;
                    }

                    c = f / z;
                    s = h / z;
                    f = c * g + s * y;
                    x = -s * g + c * y;
                }

                e[l - 1] = 0.0;
                e[k - 1] = f;
                q[k - 1] = x;
            }
        }

    }
}