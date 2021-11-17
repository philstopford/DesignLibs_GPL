using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZTRDI
{
    public static int ztrdi(ref Complex[] t, int ldt, int n, ref Complex[] det,
            int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZTRDI computes the determinant and inverse of a complex triangular matrix.
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
        //    Input/output, Complex T[LDT*N], the triangular matrix.  The zero
        //    elements of the matrix are not referenced, and the corresponding
        //    elements of the array can be used to store other information.
        //    On output, if an inverse was requested, then T has been overwritten
        //    by its inverse.
        //
        //    Input, int LDT, the leading dimension of T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int JOB.
        //    010, no determinant,    inverse, matrix is lower triangular.
        //    011, no determinant,    inverse, matrix is upper triangular.
        //    100,    determinant, no inverse.
        //    110,    determinant,    inverse, matrix is lower triangular.
        //    111,    determinant,    inverse, matrix is upper triangular.
        //
        //    Output, Complex DET[2], the determinant of the original matrix,
        //    if requested.  Otherwise not referenced.
        //    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= typeMethods.zabs1 ( DET(1) ) < 10.0
        //    or DET(1) == 0.0.  Also, DET(2) is strictly real.
        //
        //    Output, int ZTRDI.
        //    0, an inverse was requested and the matrix is nonsingular.
        //    K, an inverse was requested, but the K-th diagonal element
        //    of T is zero.
        //
    {
        int i;
        int info = 0;
        int j;
        int k;
        Complex temp;

        if (job / 100 != 0)
        {
            det[0] = new Complex(1.0, 0.0);
            det[1] = new Complex(0.0, 0.0);

            for (i = 0; i < n; i++)
            {
                det[0] *= t[i + i * ldt];

                if (typeMethods.zabs1(det[0]) == 0.0)
                {
                    break;
                }

                while (typeMethods.zabs1(det[0]) < 1.0)
                {
                    det[0] *= new Complex(10.0, 0.0);
                    det[1] -= new Complex(1.0, 0.0);
                }

                while (10.0 <= typeMethods.zabs1(det[0]))
                {
                    det[0] /= new Complex(10.0, 0.0);
                    det[1] += new Complex(1.0, 0.0);
                }
            }
        }

        //
        //  Compute inverse of upper triangular matrix.
        //
        if (job / 10 % 10 != 0)
        {
            if (job % 10 != 0)
            {
                info = 0;

                for (k = 0; k < n; k++)
                {
                    if (typeMethods.zabs1(t[k + k * ldt]) == 0.0)
                    {
                        info = k + 1;
                        break;
                    }

                    t[k + k * ldt] = new Complex(1.0, 0.0) / t[k + k * ldt];
                    temp = -t[k + k * ldt];
                    BLAS1Z.zscal(k, temp, ref t, 1, index: +0 + k * ldt);

                    for (j = k + 1; j < n; j++)
                    {
                        temp = t[k + j * ldt];
                        t[k + j * ldt] = new Complex(0.0, 0.0);
                        BLAS1Z.zaxpy(k + 1, temp, t, 1, ref t, 1, xIndex: +0 + k * ldt, yIndex: +0 + j * ldt);
                    }
                }
            }
            //
            //  Compute inverse of lower triangular matrix.
            //
            else
            {
                info = 0;

                for (k = n - 1; 0 <= k; k--)
                {
                    if (typeMethods.zabs1(t[k + k * ldt]) == 0.0)
                    {
                        info = k + 1;
                        break;
                    }

                    t[k + k * ldt] = new Complex(1.0, 0.0) / t[k + k * ldt];

                    if (k != n - 1)
                    {
                        temp = -t[k + k * ldt];
                        BLAS1Z.zscal(n - k - 1, temp, ref t, 1, index: +k + 1 + k * ldt);
                    }

                    for (j = 0; j < k; j++)
                    {
                        temp = t[k + j * ldt];
                        t[k + j * ldt] = new Complex(0.0, 0.0);
                        BLAS1Z.zaxpy(n - k, temp, t, 1, ref t, 1, xIndex: +k + k * ldt, yIndex: +k + j * ldt);
                    }
                }
            }
        }

        return info;
    }

}