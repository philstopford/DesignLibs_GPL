using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt
{
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
            double[] a_lu;
            double[] b;
            double c1;
            double c2;
            double cond;
            int i;
            int i1;
            int i2;
            int job;
            int[] pivot;

            i1 = -1;
            c1 = 0.0;
            //
            //  Factor the matrix.
            //
            a_lu = typeMethods.r8mat_copy_new(n, n, a);

            pivot = new int[n];

            typeMethods.r8ge_fa(n, ref a_lu, ref pivot);

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = 1.0 / (double) (n);
            }

            while (true)
            {
                job = 0;
                typeMethods.r8ge_sl(n, a_lu, pivot, ref b, job);

                c2 = typeMethods.r8vec_norm_l1(n, b);

                for (i = 0; i < n; i++)
                {
                    b[i] = typeMethods.r8_sign(b[i]);
                }

                job = 1;
                typeMethods.r8ge_sl(n, a_lu, pivot, ref b, job);

                i2 = typeMethods.r8vec_max_abs_index(n, b);

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

            cond = c2 * typeMethods.r8mat_norm_l1(n, n, a);

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
            double anorm;
            double cond;
            double ek;
            int i;
            int info;
            int j;
            int k;
            int l;
            int[] pivot;
            double s;
            double sm;
            double t;
            double wk;
            double wkm;
            double ynorm;
            double[] z;
            //
            //  Compute the L1 norm of A.
            //
            anorm = 0.0;
            for (j = 0; j < n; j++)
            {
                s = 0.0;
                for (i = 0; i < n; i++)
                {
                    s = s + Math.Abs(a[i + j * n]);
                }

                anorm = Math.Max(anorm, s);
            }

            //
            //  Compute the LU factorization.
            //
            pivot = new int[n];

            info = typeMethods.r8ge_fa(n, ref a, ref pivot);

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
            ek = 1.0;
            z = new double[n];
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

                wk = ek - z[k];
                wkm = -ek - z[k];
                s = Math.Abs(wk);
                sm = Math.Abs(wkm);

                if (a[k + k * n] != 0.0)
                {
                    wk = wk / a[k + k * n];
                    wkm = wkm / a[k + k * n];
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
                        sm = sm + Math.Abs(z[j] + wkm * a[k + j * n]);
                        z[j] = z[j] + wk * a[k + j * n];
                        s = s + Math.Abs(z[j]);
                    }

                    if (s < sm)
                    {
                        t = wkm - wk;
                        wk = wkm;
                        for (j = k + 1; j < n; j++)
                        {
                            z[j] = z[j] + t * a[k + j * n];
                        }
                    }
                }

                z[k] = wk;
            }

            s = 0.0;
            for (i = 0; i < n; i++)
            {
                s = s + Math.Abs(z[i]);
            }

            for (i = 0; i < n; i++)
            {
                z[i] = z[i] / s;
            }

            //
            //  Solve L' * Y = W
            //
            for (k = n - 1; 0 <= k; k--)
            {
                for (i = k + 1; i < n; i++)
                {
                    z[k] = z[k] + z[i] * a[i + k * n];
                }

                t = Math.Abs(z[k]);
                if (1.0 < t)
                {
                    for (i = 0; i < n; i++)
                    {
                        z[i] = z[i] / t;
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
                t = t + Math.Abs(z[i]);
            }

            for (i = 0; i < n; i++)
            {
                z[i] = z[i] / t;
            }

            ynorm = 1.0;
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
                    z[i] = z[i] + t * a[i + k * n];
                }

                t = Math.Abs(z[k]);

                if (1.0 < t)
                {
                    ynorm = ynorm / t;
                    for (i = 0; i < n; i++)
                    {
                        z[i] = z[i] / t;
                    }
                }
            }

            s = 0.0;
            for (i = 0; i < n; i++)
            {
                s = s + Math.Abs(z[i]);
            }

            for (i = 0; i < n; i++)
            {
                z[i] = z[i] / s;
            }

            ynorm = ynorm / s;
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
                    z[k] = z[k] / a[k + k * n];
                }
                else
                {
                    z[k] = 1.0;
                }

                for (i = 0; i < k; i++)
                {
                    z[i] = z[i] - a[i + k * n] * z[k];
                }
            }

            //
            //  Normalize Z in the L1 norm.
            //
            s = 0.0;
            for (i = 0; i < n; i++)
            {
                s = s + Math.Abs(z[i]);
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
            double a_norm;
            double ainv_norm;
            double[] ax;
            double ax_norm;
            double cond;
            int i;
            int seed;
            double[] x;
            double x_norm;

            a_norm = 0.0;
            ainv_norm = 0.0;
            seed = 123456789;

            for (i = 1; i <= m; i++)
            {
                x = UniformRNG.r8vec_uniform_unit_new(n, ref seed);
                x_norm = typeMethods.r8vec_norm_l1(n, x);
                ax = typeMethods.r8mat_mv_new(n, n, a, x);
                ax_norm = typeMethods.r8vec_norm_l1(n, ax);

                if (ax_norm == 0.0)
                {
                    cond = 0.0;
                    return cond;
                }

                a_norm = Math.Max(a_norm, ax_norm / x_norm);
                ainv_norm = Math.Max(ainv_norm, x_norm / ax_norm);

            }

            cond = a_norm * ainv_norm;

            return cond;
        }
    }
}