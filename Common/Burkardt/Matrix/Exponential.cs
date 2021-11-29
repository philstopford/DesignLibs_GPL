using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

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
        int k;
        const double one = 1.0;
        const int q = 6;

        Complex[] a2 = typeMethods.c8mat_copy_new(n, n, a);

        double a_norm = typeMethods.c8mat_norm_li(n, n, a2);

        int ee = (int) typeMethods.r8_log_2(a_norm) + 1;

        int s = Math.Max(0, ee + 1);

        double t = 1.0 / Math.Pow(2.0, s);

        typeMethods.c8mat_scale_r8(n, n, t, ref a2);

        Complex[] x = typeMethods.c8mat_copy_new(n, n, a2);

        double c = 0.5;

        Complex[] e = typeMethods.c8mat_identity_new(n);

        typeMethods.c8mat_add_r8(n, n, one, e, c, a2, ref e);

        Complex[] d = typeMethods.c8mat_identity_new(n);

        typeMethods.c8mat_add_r8(n, n, one, d, -c, a2, ref d);

        bool p = true;

        for (k = 2; k <= q; k++)
        {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));

            typeMethods.c8mat_mm(n, n, n, a2, x, ref x);

            typeMethods.c8mat_add_r8(n, n, c, x, one, e, ref e);

            switch (p)
            {
                case true:
                    typeMethods.c8mat_add_r8(n, n, c, x, one, d, ref d);
                    break;
                default:
                    typeMethods.c8mat_add_r8(n, n, -c, x, one, d, ref d);
                    break;
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
        int k;
        const double one = 1.0;
        const int q = 6;

        double[] a2 = typeMethods.r8mat_copy_new(n, n, a);

        double a_norm = typeMethods.r8mat_norm_li(n, n, a2);

        int ee = (int) typeMethods.r8_log_2(a_norm) + 1;

        int s = Math.Max(0, ee + 1);

        double t = 1.0 / Math.Pow(2.0, s);

        typeMethods.r8mat_scale(n, n, t, ref a2);

        double[] x = typeMethods.r8mat_copy_new(n, n, a2);

        double c = 0.5;

        double[] e = typeMethods.r8mat_identity_new(n);

        typeMethods.r8mat_add(n, n, one, e, c, a2, e);

        double[] d = typeMethods.r8mat_identity_new(n);

        typeMethods.r8mat_add(n, n, one, d, -c, a2, d);

        bool p = true;

        for (k = 2; k <= q; k++)
        {
            c = c * (q - k + 1) / (k * (2 * q - k + 1));

            typeMethods.r8mat_mm(n, n, n, a2, x, ref x);

            typeMethods.r8mat_add(n, n, c, x, one, e, e);

            switch (p)
            {
                case true:
                    typeMethods.r8mat_add(n, n, c, x, one, d, d);
                    break;
                default:
                    typeMethods.r8mat_add(n, n, -c, x, one, d, d);
                    break;
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
        const double one = 1.0;

        double[] e = typeMethods.r8mat_zeros_new(n, n);

        double[] f = typeMethods.r8mat_identity_new(n);

        int k = 1;

        while (typeMethods.r8mat_is_significant(n, n, e, f))
        {
            typeMethods.r8mat_add(n, n, one, e, one, f, e);

            typeMethods.r8mat_mm(n, n, n, a, f, ref f);

            double s = 1.0 / k;

            typeMethods.r8mat_scale(n, n, s, ref f);

            k += 1;
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
        //
        //  [ V, D ] = eig ( A );
        //  E = V * diag ( exp ( diag ( D ) ) ) / V;
        //
        return null;
    }
}