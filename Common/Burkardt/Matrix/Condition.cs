﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double condition_hager(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_HAGER estimates the L1 condition number of a matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Condition Estimates,
        //    SIAM Journal on Scientific and Statistical Computing,
        //    Volume 5, Number 2, June 1984, pages 311-316.
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Output, double CONDITION_HAGER, the estimated L1 condition.
        //
    {
        double c2;
        int i;

        int i1 = -1;
        double c1 = 0.0;
        //
        //  Factor the matrix.
        //
        double[] a_lu = typeMethods.r8mat_copy_new(n, n, a);

        int[] pivot = new int[n];

        typeMethods.r8ge_fa(n, ref a_lu, ref pivot);

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 1.0 / n;
        }

        while (true)
        {
            int job = 0;
            typeMethods.r8ge_sl(n, a_lu, pivot, ref b, job);

            c2 = typeMethods.r8vec_norm_l1(n, b);

            for (i = 0; i < n; i++)
            {
                b[i] = typeMethods.r8_sign(b[i]);
            }

            job = 1;
            typeMethods.r8ge_sl(n, a_lu, pivot, ref b, job);

            int i2 = typeMethods.r8vec_max_abs_index(n, b);

            if (0 <= i1)
            {
                if (i1 == i2 || c2 <= c1)
                {
                    break;
                }
            }

            i1 = i2;
            c1 = c2;

            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            b[i1] = 1.0;
        }

        double cond = c2 * typeMethods.r8mat_norm_l1(n, n, a);

        return cond;
    }

    public static double condition_linpack(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_LINPACK estimates the L1 condition number.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    For the system A * X = B, relative perturbations in A and B
        //    of size EPSILON may cause relative perturbations in X of size
        //    EPSILON * COND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input/output, double A[N*N].  On input, a matrix to be factored.
        //    On output, the LU factorization of the matrix.
        //
        //    Output, double CONDITION_LINPACK, the estimated L1 condition.
        //
    {
        double cond;
        int i;
        int j;
        int k;
        int l;
        double s;
        double t;
        //
        //  Compute the L1 norm of A.
        //
        double anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            s = 0.0;
            for (i = 0; i < n; i++)
            {
                s += Math.Abs(a[i + j * n]);
            }

            anorm = Math.Max(anorm, s);
        }

        //
        //  Compute the LU factorization.
        //
        int[] pivot = new int[n];

        int info = typeMethods.r8ge_fa(n, ref a, ref pivot);

        if (info != 0)
        {
            cond = 0.0;
            return cond;
        }

        //
        //  COND = norm(A) * (estimate of norm(inverse(A)))
        //
        //  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
        //
        //  where
        //    A * Z = Y
        //  and
        //    A' * Y = E
        //
        //  The components of E are chosen to cause maximum local growth in the
        //  elements of W, where U'*W = E.  The vectors are frequently rescaled
        //  to avoid overflow.
        //
        //  Solve U' * W = E.
        //
        double ek = 1.0;
        double[] z = new double[n];
        for (i = 0; i < n; i++)
        {
            z[i] = 0.0;
        }

        for (k = 0; k < n; k++)
        {
            if (z[k] != 0.0)
            {
                ek = -typeMethods.r8_sign2(ek, z[k]);
            }

            if (Math.Abs(a[k + k * n]) < Math.Abs(ek - z[k]))
            {
                s = Math.Abs(a[k + k * n]) / Math.Abs(ek - z[k]);
                for (i = 0; i < n; i++)
                {
                    z[i] = s * z[i];
                }

                ek = s * ek;
            }

            double wk = ek - z[k];
            double wkm = -ek - z[k];
            s = Math.Abs(wk);
            double sm = Math.Abs(wkm);

            if (a[k + k * n] != 0.0)
            {
                wk /= a[k + k * n];
                wkm /= a[k + k * n];
            }
            else
            {
                wk = 1.0;
                wkm = 1.0;
            }

            if (k + 2 <= n)
            {
                for (j = k + 1; j < n; j++)
                {
                    sm += Math.Abs(z[j] + wkm * a[k + j * n]);
                    z[j] += wk * a[k + j * n];
                    s += Math.Abs(z[j]);
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    for (j = k + 1; j < n; j++)
                    {
                        z[j] += t * a[k + j * n];
                    }
                }
            }

            z[k] = wk;
        }

        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= s;
        }

        //
        //  Solve L' * Y = W
        //
        for (k = n - 1; 0 <= k; k--)
        {
            for (i = k + 1; i < n; i++)
            {
                z[k] += z[i] * a[i + k * n];
            }

            t = Math.Abs(z[k]);
            switch (t)
            {
                case > 1.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        z[i] /= t;
                    }

                    break;
                }
            }

            l = pivot[k] - 1;

            t = z[l];
            z[l] = z[k];
            z[k] = t;
        }

        t = 0.0;
        for (i = 0; i < n; i++)
        {
            t += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= t;
        }

        double ynorm = 1.0;
        //
        //  Solve L * V = Y.
        //
        for (k = 0; k < n; k++)
        {
            l = pivot[k] - 1;

            t = z[l];
            z[l] = z[k];
            z[k] = t;

            for (i = k + 1; i < n; i++)
            {
                z[i] += t * a[i + k * n];
            }

            t = Math.Abs(z[k]);

            switch (t)
            {
                case > 1.0:
                {
                    ynorm /= t;
                    for (i = 0; i < n; i++)
                    {
                        z[i] /= t;
                    }

                    break;
                }
            }
        }

        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= s;
        }

        ynorm /= s;
        //
        //  Solve U * Z = V.
        //
        for (k = n - 1; 0 <= k; k--)
        {
            if (Math.Abs(a[k + k * n]) < Math.Abs(z[k]))
            {
                s = Math.Abs(a[k + k * n]) / Math.Abs(z[k]);
                for (i = 0; i < n; i++)
                {
                    z[i] = s * z[i];
                }

                ynorm = s * ynorm;
            }

            if (a[k + k * n] != 0.0)
            {
                z[k] /= a[k + k * n];
            }
            else
            {
                z[k] = 1.0;
            }

            for (i = 0; i < k; i++)
            {
                z[i] -= a[i + k * n] * z[k];
            }
        }

        //
        //  Normalize Z in the L1 norm.
        //
        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        s = 1.0 / s;

        for (i = 0; i < n; i++)
        {
            z[i] = s * z[i];
        }

        ynorm = s * ynorm;

        cond = anorm / ynorm;

        return cond;
    }

    public static double condition_sample1(int n, double[] a, int m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONDITION_SAMPLE1 estimates the L1 condition number of a matrix.
        //
        //  Discussion:
        //
        //    A naive sampling method is used.
        //
        //    Only "forward" sampling is used, that is, we only look at results
        //    of the form y=A*x.
        //
        //    Presumably, solving systems A*y=x would give us a better idea of 
        //    the inverse matrix.
        //
        //    Moreover, a power sequence y1 = A*x, y2 = A*y1, ... and the same for
        //    the inverse might work better too.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, int M, the number of samples to use.
        //
        //    Output, double CONDITION_SAMPLE1, the estimated L1 condition.
        //
    {
        double cond;
        int i;

        double a_norm = 0.0;
        double ainv_norm = 0.0;
        int seed = 123456789;

        typeMethods.r8vecNormalData data = new();

        for (i = 1; i <= m; i++)
        {
            double[] x = UniformRNG.r8vec_uniform_unit_new(n, ref data, ref seed);
            double x_norm = typeMethods.r8vec_norm_l1(n, x);
            double[] ax = typeMethods.r8mat_mv_new(n, n, a, x);
            double ax_norm = typeMethods.r8vec_norm_l1(n, ax);

            switch (ax_norm)
            {
                case 0.0:
                    cond = 0.0;
                    return cond;
                default:
                    a_norm = Math.Max(a_norm, ax_norm / x_norm);
                    ainv_norm = Math.Max(ainv_norm, x_norm / ax_norm);
                    break;
            }
        }

        cond = a_norm * ainv_norm;

        return cond;
    }
}