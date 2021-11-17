using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DTRDI
{
    public static int dtrdi(ref double[] t, int ldt, int n, ref double[] det, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRDI computes the determinant and inverse of a real triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 March 2005
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
        //    Input/output, double T[LDT*N].
        //    On input, T contains the triangular matrix.  The zero elements of the
        //    matrix are not referenced, and the corresponding elements of the array
        //    can be used to store other information.
        //    On output, T contains the inverse matrix, if it was requested.
        //
        //    Input, int LDT, the leading dimension of T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double DET[2], the determinant of the matrix, if
        //    requested.  The determinant = DET[0] * 10.0**DET[1], with
        //    1.0 <= abs ( DET[0] ) < 10.0, or DET[0] == 0.
        //
        //    Input, int JOB, specifies the shape of T, and the task.
        //    010, inverse of lower triangular matrix.
        //    011, inverse of upper triangular matrix.
        //    100, determinant only.
        //    110, determinant and inverse of lower triangular.
        //    111, determinant and inverse of upper triangular.
        //
        //    Output, int DTRDI.
        //    If the inverse was requested, then
        //    0, if the system was nonsingular;
        //    nonzero, if the system was singular.
        //
    {
        int i;
        int info;
        int j;
        int k;
        double temp;
        //
        //  Determinant.
        //
        info = 0;

        if (job / 100 != 0)
        {
            det[0] = 1.0;
            det[1] = 0.0;

            for (i = 1; i <= n; i++)
            {
                det[0] *= t[i - 1 + (i - 1) * ldt];

                if (det[0] == 0.0)
                {
                    break;
                }

                while (Math.Abs(det[0]) < 1.0)
                {
                    det[0] *= 10.0;
                    det[1] -= 1.0;
                }

                while (10.0 <= Math.Abs(det[0]))
                {
                    det[0] /= 10.0;
                    det[1] += 1.0;
                }
            }
        }

        switch (job / 10 % 10)
        {
            case 0:
                return info;
        }

        //
        //  Inverse of an upper triangular matrix.
        //
        if (job % 10 != 0)
        {
            info = 0;

            for (k = 1; k <= n; k++)
            {
                if (t[k - 1 + (k - 1) * ldt] == 0.0)
                {
                    info = k;
                    break;
                }

                t[k - 1 + (k - 1) * ldt] = 1.0 / t[k - 1 + (k - 1) * ldt];
                temp = -t[k - 1 + (k - 1) * ldt];
                BLAS1D.dscal(k - 1, temp, ref t, 1, index: +0 + (k - 1) * ldt);

                for (j = k + 1; j <= n; j++)
                {
                    temp = t[k - 1 + (j - 1) * ldt];
                    t[k - 1 + (j - 1) * ldt] = 0.0;
                    BLAS1D.daxpy(k, temp, t, 1, ref t, 1, xIndex: +0 + (k - 1) * ldt, yIndex: +0 + (j - 1) * ldt);
                }
            }
        }
        //
        //  Inverse of a lower triangular matrix.
        //
        else
        {
            info = 0;

            for (k = n; 1 <= k; k--)
            {
                if (t[k - 1 + (k - 1) * ldt] == 0.0)
                {
                    info = k;
                    break;
                }

                t[k - 1 + (k - 1) * ldt] = 1.0 / t[k - 1 + (k - 1) * ldt];
                temp = -t[k - 1 + (k - 1) * ldt];

                if (k != n)
                {
                    BLAS1D.dscal(n - k, temp, ref t, 1, index: +k + (k - 1) * ldt);
                }

                for (j = 1; j <= k - 1; j++)
                {
                    temp = t[k - 1 + (j - 1) * ldt];
                    t[k - 1 + (j - 1) * ldt] = 0.0;
                    BLAS1D.daxpy(n - k + 1, temp, t, 1, ref t, 1, xIndex: +k - 1 + (k - 1) * ldt,
                        yIndex: +k - 1 + (j - 1) * ldt);
                }
            }
        }

        return info;
    }

}