﻿using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt
{
    public static partial class Matrix
    {
        public static Complex[] c8mat_expm1(int n, Complex[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 March 2013
            //
            //  Author:
            //
            //    Cleve Moler, Charles Van Loan
            //
            //  Reference:
            //
            //    Cleve Moler, Charles VanLoan,
            //    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
            //    Twenty-Five Years Later,
            //    SIAM Review,
            //    Volume 45, Number 1, March 2003, pages 3-49.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the matrix.
            //
            //    Input, Complex A[N*N], the matrix.
            //
            //    Output, Complex C8MAT_EXPM1[N*N], the estimate for exp(A).
            //
        {
            Complex[] a2;
            double a_norm;
            double c;
            Complex[] d;
            Complex[] e;
            int ee;
            int k;
            const double one = 1.0;
            bool p;
            const int q = 6;
            int s;
            double t;
            Complex[] x;

            a2 = typeMethods.c8mat_copy_new(n, n, a);

            a_norm = typeMethods.c8mat_norm_li(n, n, a2);

            ee = (int) (typeMethods.r8_log_2(a_norm)) + 1;

            s = Math.Max(0, ee + 1);

            t = 1.0 / Math.Pow(2.0, s);

            typeMethods.c8mat_scale_r8(n, n, t, ref a2);

            x = typeMethods.c8mat_copy_new(n, n, a2);

            c = 0.5;

            e = typeMethods.c8mat_identity_new(n);

            typeMethods.c8mat_add_r8(n, n, one, e, c, a2, ref e);

            d = typeMethods.c8mat_identity_new(n);

            typeMethods.c8mat_add_r8(n, n, one, d, -c, a2, ref d);

            p = true;

            for (k = 2; k <= q; k++)
            {
                c = c * (double) (q - k + 1) / (double) (k * (2 * q - k + 1));

                typeMethods.c8mat_mm(n, n, n, a2, x, ref x);

                typeMethods.c8mat_add_r8(n, n, c, x, one, e, ref e);

                if (p)
                {
                    typeMethods.c8mat_add_r8(n, n, c, x, one, d, ref d);
                }
                else
                {
                    typeMethods.c8mat_add_r8(n, n, -c, x, one, d, ref d);
                }

                p = !p;
            }

            //
            //  E -> inverse(D) * E
            //
            typeMethods.c8mat_minvm(n, n, d, e, ref e);
            //
            //  E -> E^(2*S)
            //
            for (k = 1; k <= s; k++)
            {
                typeMethods.c8mat_mm(n, n, n, e, e, ref e);
            }

            return e;
        }

        public static double[] r8mat_expm1(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_EXPM1 is essentially MATLAB's built-in matrix exponential algorithm.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 December 2011
            //
            //  Author:
            //
            //    Cleve Moler, Charles Van Loan
            //
            //  Reference:
            //
            //    Cleve Moler, Charles VanLoan,
            //    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
            //    Twenty-Five Years Later,
            //    SIAM Review,
            //    Volume 45, Number 1, March 2003, pages 3-49.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the matrix.
            //
            //    Input, double A[N*N], the matrix.
            //
            //    Output, double R8MAT_EXPM1[N*N], the estimate for exp(A).
            //
        {
            double[] a2;
            double a_norm;
            double c;
            double[] d;
            double[] e;
            int ee;
            int k;
            const double one = 1.0;
            bool p;
            const int q = 6;
            int s;
            double t;
            double[] x;

            a2 = typeMethods.r8mat_copy_new(n, n, a);

            a_norm = typeMethods.r8mat_norm_li(n, n, a2);

            ee = (int) (typeMethods.r8_log_2(a_norm)) + 1;

            s = Math.Max(0, ee + 1);

            t = 1.0 / Math.Pow(2.0, s);

            typeMethods.r8mat_scale(n, n, t, ref a2);

            x = typeMethods.r8mat_copy_new(n, n, a2);

            c = 0.5;

            e = typeMethods.r8mat_identity_new(n);

            typeMethods.r8mat_add(n, n, one, e, c, a2, e);

            d = typeMethods.r8mat_identity_new(n);

            typeMethods.r8mat_add(n, n, one, d, -c, a2, d);

            p = true;

            for (k = 2; k <= q; k++)
            {
                c = c * (double) (q - k + 1) / (double) (k * (2 * q - k + 1));

                typeMethods.r8mat_mm(n, n, n, a2, x, ref x);

                typeMethods.r8mat_add(n, n, c, x, one, e, e);

                if (p)
                {
                    typeMethods.r8mat_add(n, n, c, x, one, d, d);
                }
                else
                {
                    typeMethods.r8mat_add(n, n, -c, x, one, d, d);
                }

                p = !p;
            }

            //
            //  E -> inverse(D) * E
            //
            typeMethods.r8mat_minvm(n, n, d, e, ref e);
            //
            //  E -> E^(2*S)
            //
            for (k = 1; k <= s; k++)
            {
                typeMethods.r8mat_mm(n, n, n, e, e, ref e);
            }

            return e;
        }

        public static double[] r8mat_expm2(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_EXPM2 uses the Taylor series for the matrix exponential.
            //
            //  Discussion:
            //
            //    Formally,
            //
            //      exp ( A ) = I + A + 1/2 A^2 + 1/3! A^3 + ...
            //
            //    This function sums the series until a tolerance is satisfied.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 December 2011
            //
            //  Author:
            //
            //    Cleve Moler, Charles Van Loan
            //
            //  Reference:
            //
            //    Cleve Moler, Charles VanLoan,
            //    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
            //    Twenty-Five Years Later,
            //    SIAM Review,
            //    Volume 45, Number 1, March 2003, pages 3-49.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the matrix.
            //
            //    Input, double A[N*N], the matrix.
            //
            //    Output, double R8MAT_EXPM2[N*N], the estimate for exp(A).
            //
        {
            double[] e;
            double[] f;
            int k;
            const double one = 1.0;
            double s;

            e = typeMethods.r8mat_zeros_new(n, n);

            f = typeMethods.r8mat_identity_new(n);

            k = 1;

            while (typeMethods.r8mat_is_significant(n, n, e, f))
            {
                typeMethods.r8mat_add(n, n, one, e, one, f, e);

                typeMethods.r8mat_mm(n, n, n, a, f, ref f);

                s = 1.0 / (double) (k);

                typeMethods.r8mat_scale(n, n, s, ref f);

                k = k + 1;
            }

            return e;
        }

        public static double[] r8mat_expm3(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_EXPM3 approximates the matrix exponential using an eigenvalue approach.
            //
            //  Discussion:
            //
            //    exp(A) = V * D * inv(V)
            //
            //    where V is the matrix of eigenvectors of A, and D is the diagonal matrix
            //    whose i-th diagonal entry is exp(lambda(i)), for lambda(i) an eigenvalue
            //    of A.
            //
            //    This function is accurate for matrices which are symmetric, orthogonal,
            //    or normal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 December 2011
            //
            //  Author:
            //
            //    Cleve Moler, Charles Van Loan
            //
            //  Reference:
            //
            //    Cleve Moler, Charles VanLoan,
            //    Nineteen Dubious Ways to Compute the Exponential of a Matrix,
            //    Twenty-Five Years Later,
            //    SIAM Review,
            //    Volume 45, Number 1, March 2003, pages 3-49.
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the matrix.
            //
            //    Input, double A[N*N], the matrix.
            //
            //    Output, double R8MAT_EXPM3[N*N], the estimate for exp(A).
            //
        {
            double[] e = null;
            //
            //  [ V, D ] = eig ( A );
            //  E = V * diag ( exp ( diag ( D ) ) ) / V;
            //
            return e;
        }
    }
}