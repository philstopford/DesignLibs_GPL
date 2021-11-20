using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZTRSL
{
    public static int ztrsl(Complex[] t, int ldt, int n, ref Complex[] b,
            int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZTRSL solves triangular systems T*X=B or Hermitian(T)*X=B.
        //
        //  Discussion:
        //
        //    Hermitian ( T ) denotes the Complex.Conjugateugate transpose of the matrix T.
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
        //    Input, Complex T[LDT*N], the matrix of the system.  The zero
        //    elements of the matrix are not referenced, and the corresponding
        //    elements of the array can be used to store other information.
        //
        //    Input, int LDT, the leading dimension of T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Input, int JOB, specifies what kind of system is to be solved.
        //    00, solve T*X=B, T lower triangular,
        //    01, solve T*X=B, T upper triangular,
        //    10, solve hermitian(T)*X=B, T lower triangular,
        //    11, solve hermitian(T)*X=B, T upper triangular.
        //
        //    Output, int ZTRSL.
        //    0, the system is nonsingular.
        //    K, the index of the first zero diagonal element of T.
        //
    {
        int i;
        int info;
        int j;
        int jj;
        Complex temp;
        //
        //  Check for zero diagonal elements.
        //
        for (i = 0; i < n; i++)
        {
            if (typeMethods.zabs1(t[i + i * ldt]) != 0.0)
            {
                continue;
            }

            info = i + 1;
            return info;
        }

        info = 0;
        //
        //  Determine the task and go to it.
        //
        int kase = 1;

        if (job % 10 != 0)
        {
            kase = 2;
        }

        if (job % 100 / 10 != 0)
        {
            kase += 2;
        }

        switch (kase)
        {
            //
            //  Solve T * X = B for T lower triangular.
            //
            case 1:
            {
                b[0] /= t[0 + 0 * ldt];

                for (j = 2; j <= n; j++)
                {
                    temp = -b[j - 2];
                    BLAS1Z.zaxpy(n - j + 1, temp, t, 1, ref b, 1, xIndex: +j - 1 + (j - 2) * ldt, yIndex: +j - 1);
                    b[j - 1] /= t[j - 1 + (j - 1) * ldt];
                }

                break;
            }
            //
            //  Solve T * X = B for T upper triangular.
            //
            case 2:
            {
                b[n - 1] /= t[n - 1 + (n - 1) * ldt];

                for (jj = 2; jj <= n; jj++)
                {
                    j = n - jj + 1;
                    temp = -b[j];
                    BLAS1Z.zaxpy(j, temp, t, 1, ref b, 1, xIndex: +0 + j * ldt);
                    b[j - 1] /= t[j - 1 + (j - 1) * ldt];
                }

                break;
            }
            //
            //  Solve hermitian(T) * X = B for T lower triangular.
            //
            case 3:
            {
                b[n - 1] /= Complex.Conjugate(t[n - 1 + (n - 1) * ldt]);

                for (jj = 2; jj <= n; jj++)
                {
                    j = n - jj + 1;
                    b[j - 1] -= BLAS1Z.zdotc(jj - 1, t, 1, b, 1, xIndex: +j + (j - 1) * ldt, yIndex: +j);
                    b[j - 1] /= Complex.Conjugate(t[j - 1 + (j - 1) * ldt]);
                }

                break;
            }
            //
            //  Solve hermitian(T) * X = B for T upper triangular.
            //
            case 4:
            {
                b[0] /= Complex.Conjugate(t[0 + 0 * ldt]);

                for (j = 2; j <= n; j++)
                {
                    b[j - 1] -= BLAS1Z.zdotc(j - 1, t, 1, b, 1, xIndex: +0 + (j - 1) * ldt);
                    b[j - 1] /= Complex.Conjugate(t[j - 1 + (j - 1) * ldt]);
                }

                break;
            }
        }

        return info;
    }

}