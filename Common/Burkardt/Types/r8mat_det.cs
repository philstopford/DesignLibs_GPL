using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8mat_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET computes the determinant of an R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 October 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Helmut Spaeth
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Helmut Spaeth,
        //    Cluster Analysis Algorithms
        //    for Data Reduction and Classification of Objects,
        //    Ellis Horwood, 1980, page 125-127.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix whose determinant is desired.
        //
        //    Output, double R8MAT_DET, the determinant of the matrix.
        //
    {
        int i;
        int j;
        int k;

        double[] b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = a[i + j * n];
            }
        }

        double det = 1.0;

        for (k = 1; k <= n; k++)
        {
            int m = k;
            int kk;
            for (kk = k + 1; kk <= n; kk++)
            {
                if (Math.Abs(b[m - 1 + (k - 1) * n]) < Math.Abs(b[kk - 1 + (k - 1) * n]))
                {
                    m = kk;
                }
            }

            double temp;
            if (m != k)
            {
                det = -det;

                temp = b[m - 1 + (k - 1) * n];
                b[m - 1 + (k - 1) * n] = b[k - 1 + (k - 1) * n];
                b[k - 1 + (k - 1) * n] = temp;
            }

            det *= b[k - 1 + (k - 1) * n];

            if (b[k - 1 + (k - 1) * n] == 0.0)
            {
                continue;
            }

            for (i = k + 1; i <= n; i++)
            {
                b[i - 1 + (k - 1) * n] = -b[i - 1 + (k - 1) * n] / b[k - 1 + (k - 1) * n];
            }

            for (j = k + 1; j <= n; j++)
            {
                if (m != k)
                {
                    temp = b[m - 1 + (j - 1) * n];
                    b[m - 1 + (j - 1) * n] = b[k - 1 + (j - 1) * n];
                    b[k - 1 + (j - 1) * n] = temp;
                }

                for (i = k + 1; i <= n; i++)
                {
                    b[i - 1 + (j - 1) * n] += b[i - 1 + (k - 1) * n] * b[k - 1 + (j - 1) * n];
                }
            }
        }

        return det;
    }

    public static double r8mat_det_2d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Discussion:
        //
        //    The determinant of a 2 by 2 matrix is
        //
        //      a11 * a22 - a12 * a21.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2*2], the matrix whose determinant is desired.
        //
        //    Output, double R8MAT_DET_2D, the determinant of the matrix.
        //
    {
        double det = a[0 + 0 * 2] * a[1 + 1 * 2] - a[0 + 1 * 2] * a[1 + 0 * 2];

        return det;
    }

    public static double r8mat_det_3d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The determinant of a 3 by 3 matrix is
        //
        //        a11 * a22 * a33 - a11 * a23 * a32
        //      + a12 * a23 * a31 - a12 * a21 * a33
        //      + a13 * a21 * a32 - a13 * a22 * a31
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3*3], the matrix whose determinant is desired.
        //
        //    Output, double R8MAT_DET_3D, the determinant of the matrix.
        //
    {
        double det = a[0 + 0 * 3] * (a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3])
                     + a[0 + 1 * 3] * (a[1 + 2 * 3] * a[2 + 0 * 3] - a[1 + 0 * 3] * a[2 + 2 * 3])
                     + a[0 + 2 * 3] * (a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]);

        return det;
    }

    public static double r8mat_det_4d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[4*4], the matrix whose determinant is desired.
        //
        //    Output, double R8MAT_DET_4D, the determinant of the matrix.
        //
    {
        double det = a[0 + 0 * 4] * (
                         a[1 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                         - a[1 + 2 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4])
                         + a[1 + 3 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4]))
                     - a[0 + 1 * 4] * (
                         a[1 + 0 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                         - a[1 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 0 * 4])
                         + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4]))
                     + a[0 + 2 * 4] * (
                         a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4])
                         - a[1 + 1 * 4] * (a[2 + 0 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 0 * 4])
                         + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4]))
                     - a[0 + 3 * 4] * (
                         a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])
                         - a[1 + 1 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4])
                         + a[1 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4]));

        return det;
    }

    public static double r8mat_det_5d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[5*5], the matrix whose determinant is desired.
        //
        //    Output, double R8MAT_DET_5D, the determinant of the matrix.
        //
    {
        double[] b = new double[4 * 4];
        int k;
        //
        //  Expand the determinant into the sum of the determinants of the
        //  five 4 by 4 matrices created by dropping row 1, and column k.
        //
        double det = 0.0;
        double sign = 1.0;

        for (k = 0; k < 5; k++)
        {
            int i;
            for (i = 0; i < 4; i++)
            {
                int j;
                for (j = 0; j < 4; j++)
                {
                    int inc;
                    inc = j < k ? 0 : 1;

                    b[i + j * 4] = a[i + 1 + (j + inc) * 5];
                }
            }

            det += sign * a[0 + k * 5] * r8mat_det_4d(b);

            sign = -sign;
        }

        return det;
    }
        
}