using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DTRSL
{
    public static int dtrsl(double[] t, int ldt, int n, ref double[] b, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRSL solves triangular linear systems.
        //
        //  Discussion:
        //
        //    DTRSL can solve T * X = B or T' * X = B where T is a triangular
        //    matrix of order N.
        //
        //    Here T' denotes the transpose of the matrix T.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 May 2005
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
        //    Input, double T[LDT*N], the matrix of the system.  The zero
        //    elements of the matrix are not referenced, and the corresponding
        //    elements of the array can be used to store other information.
        //
        //    Input, int LDT, the leading dimension of the array T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Input, int JOB, specifies what kind of system is to be solved:
        //    00, solve T * X = B, T lower triangular,
        //    01, solve T * X = B, T upper triangular,
        //    10, solve T'* X = B, T lower triangular,
        //    11, solve T'* X = B, T upper triangular.
        //
        //    Output, int DTRSL, singularity indicator.
        //    0, the system is nonsingular.
        //    nonzero, the index of the first zero diagonal element of T.
        //
    {
        int info;
        int j;
        int jj;
        int kase;
        double temp;
        //
        //  Check for zero diagonal elements.
        //
        for (j = 1; j <= n; j++)
        {
            switch (t[j - 1 + (j - 1) * ldt])
            {
                case 0.0:
                    info = j;
                    return info;
            }
        }

        info = 0;
        kase = (job % 10) switch
        {
            //
            //  Determine the task and go to it.
            //
            0 => 1,
            _ => 2
        };

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
                    BLAS1D.daxpy(n - j + 1, temp, t, 1, ref b, 1, xIndex: +(j - 1) + (j - 2) * ldt, yIndex: +j - 1);
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
                    BLAS1D.daxpy(j, temp, t, 1, ref b, 1, xIndex: +0 + j * ldt);
                    b[j - 1] /= t[j - 1 + (j - 1) * ldt];
                }

                break;
            }
            //
            //  Solve T' * X = B for T lower triangular.
            //
            case 3:
            {
                b[n - 1] /= t[n - 1 + (n - 1) * ldt];
                for (jj = 2; jj <= n; jj++)
                {
                    j = n - jj + 1;
                    b[j - 1] -= BLAS1D.ddot(jj - 1, t, 1, b, 1, xIndex: +j + (j - 1) * ldt, yIndex: +j);
                    b[j - 1] /= t[j - 1 + (j - 1) * ldt];
                }

                break;
            }
            //
            //  Solve T' * X = B for T upper triangular.
            //
            case 4:
            {
                b[0] /= t[0 + 0 * ldt];
                for (j = 2; j <= n; j++)
                {
                    b[j - 1] -= BLAS1D.ddot(j - 1, t, 1, b, 1, xIndex: +0 + (j - 1) * ldt);
                    b[j - 1] /= t[j - 1 + (j - 1) * ldt];
                }

                break;
            }
        }

        return info;
    }

}