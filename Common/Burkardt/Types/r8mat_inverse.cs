namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_l1_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_L1_INVERSE inverts a unit lower triangular R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    A unit lower triangular matrix is a matrix with only 1's on the main
        //    diagonal, and only 0's above the main diagonal.
        //
        //    The inverse of a unit lower triangular matrix is also
        //    a unit lower triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6.
        //
        //  Parameters:
        //
        //    Input, int N, number of rows and columns in the matrix.
        //
        //    Input, double A[N*N], the unit lower triangular matrix.
        //
        //    Output, double R8MAT_L1_INVERSE[N*N], the inverse matrix.
        //
    {
        int i;
        int j;

        double[] b = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i < j)
                {
                    b[i + j * n] = 0.0;
                }
                else if (j == i)
                {
                    b[i + j * n] = 1.0;
                }
                else
                {
                    b[i + j * n] = 0.0;
                    int k;
                    for (k = 0; k < i; k++)
                    {
                        b[i + j * n] -= a[i + k * n] * b[k + j * n];
                    }
                }
            }
        }

        return b;
    }
        
    public static double[] r8mat_inverse_2d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_INVERSE_2D inverts a 2 by 2 matrix using Cramer's rule.
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
        //    23 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2*2], the matrix to be inverted.
        //
        //    Output, double R8MAT_INVERSE_2D[2*2], the inverse of the matrix A.
        //
    {
        //
        //  Compute the determinant of A.
        //
        double det = a[0 + 0 * 2] * a[1 + 1 * 2] - a[0 + 1 * 2] * a[1 + 0 * 2];
        switch (det)
        {
            //
            //  If the determinant is zero, bail out.
            //
            case 0.0:
                return null;
        }

        //
        //  Compute the entries of the inverse matrix using an explicit formula.
        //
        double[] b = new double[2 * 2];

        b[0 + 0 * 2] = +a[1 + 1 * 2] / det;
        b[0 + 1 * 2] = -a[0 + 1 * 2] / det;
        b[1 + 0 * 2] = -a[1 + 0 * 2] / det;
        b[1 + 1 * 2] = +a[0 + 0 * 2] / det;

        return b;
    }

    public static double[] r8mat_inverse_3d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_INVERSE_3D inverts a 3 by 3 matrix using Cramer's rule.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    If the determinant is zero, A is singular, and does not have an
        //    inverse.  In that case, the output is set to NULL.
        //
        //    If the determinant is nonzero, its value is an estimate
        //    of how nonsingular the matrix A is.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3*3], the matrix to be inverted.
        //
        //    Output, double R8MAT_INVERSE_3D[3*3], the inverse of the matrix A.
        //
    {
        //
        //  Compute the determinant of A.
        //
        double det =
            a[0 + 0 * 3] * (a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3])
            + a[0 + 1 * 3] * (a[1 + 2 * 3] * a[2 + 0 * 3] - a[1 + 0 * 3] * a[2 + 2 * 3])
            + a[0 + 2 * 3] * (a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]);

        switch (det)
        {
            case 0.0:
                return null;
        }

        double[] b = new double[3 * 3];

        b[0 + 0 * 3] = (a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3]) / det;
        b[0 + 1 * 3] = -(a[0 + 1 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 1 * 3]) / det;
        b[0 + 2 * 3] = (a[0 + 1 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 1 * 3]) / det;

        b[1 + 0 * 3] = -(a[1 + 0 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 0 * 3]) / det;
        b[1 + 1 * 3] = (a[0 + 0 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 0 * 3]) / det;
        b[1 + 2 * 3] = -(a[0 + 0 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 0 * 3]) / det;

        b[2 + 0 * 3] = (a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]) / det;
        b[2 + 1 * 3] = -(a[0 + 0 * 3] * a[2 + 1 * 3] - a[0 + 1 * 3] * a[2 + 0 * 3]) / det;
        b[2 + 2 * 3] = (a[0 + 0 * 3] * a[1 + 1 * 3] - a[0 + 1 * 3] * a[1 + 0 * 3]) / det;

        return b;
    }

    public static double[] r8mat_inverse_4d(double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_INVERSE_4D inverts a 4 by 4 matrix using Cramer's rule.
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
        //    18 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[4][4], the matrix to be inverted.
        //
        //    Output, double R8MAT_INVERSE_4D[4][4], the inverse of the matrix A.
        //
    {
        //
        //  Compute the determinant of A.
        //
        double det = r8mat_det_4d(a);
        switch (det)
        {
            //
            //  If the determinant is zero, bail out.
            //
            case 0.0:
                return null;
        }

        //
        //  Compute the entries of the inverse matrix using an explicit formula.
        //
        double[] b = new double[4 * 4];

        b[0 + 0 * 4] =
            +(
                +a[1 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                + a[1 + 2 * 4] * (a[2 + 3 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 3 * 4])
                + a[1 + 3 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])
            ) / det;

        b[1 + 0 * 4] =
            -(
                +a[1 + 0 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                + a[1 + 2 * 4] * (a[2 + 3 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 3 * 4])
                + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4])
            ) / det;

        b[2 + 0 * 4] =
            +(
                +a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4])
                + a[1 + 1 * 4] * (a[2 + 3 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 3 * 4])
                + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[3 + 0 * 4] =
            -(
                +a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])
                + a[1 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 2 * 4])
                + a[1 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[0 + 1 * 4] =
            -(
                +a[0 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[2 + 3 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])
            ) / det;

        b[1 + 1 * 4] =
            +(
                +a[0 + 0 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[2 + 3 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4])
            ) / det;

        b[2 + 1 * 4] =
            -(
                +a[0 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4])
                + a[0 + 1 * 4] * (a[2 + 3 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[3 + 1 * 4] =
            +(
                +a[0 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])
                + a[0 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 0 * 4] - a[2 + 0 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[0 + 2 * 4] =
            +(
                +a[0 + 1 * 4] * (a[1 + 2 * 4] * a[3 + 3 * 4] - a[1 + 3 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 3 * 4] * a[3 + 1 * 4] - a[1 + 1 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 1 * 4] * a[3 + 2 * 4] - a[1 + 2 * 4] * a[3 + 1 * 4])
            ) / det;

        b[1 + 2 * 4] =
            -(
                +a[0 + 0 * 4] * (a[1 + 2 * 4] * a[3 + 3 * 4] - a[1 + 3 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 3 * 4] * a[3 + 0 * 4] - a[1 + 0 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 0 * 4] * a[3 + 2 * 4] - a[1 + 2 * 4] * a[3 + 0 * 4])
            ) / det;

        b[2 + 2 * 4] =
            +(
                +a[0 + 0 * 4] * (a[1 + 1 * 4] * a[3 + 3 * 4] - a[1 + 3 * 4] * a[3 + 1 * 4])
                + a[0 + 1 * 4] * (a[1 + 3 * 4] * a[3 + 0 * 4] - a[1 + 0 * 4] * a[3 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 0 * 4] * a[3 + 1 * 4] - a[1 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[3 + 2 * 4] =
            -(
                +a[0 + 0 * 4] * (a[1 + 1 * 4] * a[3 + 2 * 4] - a[1 + 2 * 4] * a[3 + 1 * 4])
                + a[0 + 1 * 4] * (a[1 + 2 * 4] * a[3 + 0 * 4] - a[1 + 0 * 4] * a[3 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 0 * 4] * a[3 + 1 * 4] - a[1 + 1 * 4] * a[3 + 0 * 4])
            ) / det;

        b[0 + 3 * 4] =
            -(
                +a[0 + 1 * 4] * (a[1 + 2 * 4] * a[2 + 3 * 4] - a[1 + 3 * 4] * a[2 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 3 * 4] * a[2 + 1 * 4] - a[1 + 1 * 4] * a[2 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 1 * 4] * a[2 + 2 * 4] - a[1 + 2 * 4] * a[2 + 1 * 4])
            ) / det;

        b[1 + 3 * 4] =
            +(
                +a[0 + 0 * 4] * (a[1 + 2 * 4] * a[2 + 3 * 4] - a[1 + 3 * 4] * a[2 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 3 * 4] * a[2 + 0 * 4] - a[1 + 0 * 4] * a[2 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 0 * 4] * a[2 + 2 * 4] - a[1 + 2 * 4] * a[2 + 0 * 4])
            ) / det;

        b[2 + 3 * 4] =
            -(
                +a[0 + 0 * 4] * (a[1 + 1 * 4] * a[2 + 3 * 4] - a[1 + 3 * 4] * a[2 + 1 * 4])
                + a[0 + 1 * 4] * (a[1 + 3 * 4] * a[2 + 0 * 4] - a[1 + 0 * 4] * a[2 + 3 * 4])
                + a[0 + 3 * 4] * (a[1 + 0 * 4] * a[2 + 1 * 4] - a[1 + 1 * 4] * a[2 + 0 * 4])
            ) / det;

        b[3 + 3 * 4] =
            +(
                +a[0 + 0 * 4] * (a[1 + 1 * 4] * a[2 + 2 * 4] - a[1 + 2 * 4] * a[2 + 1 * 4])
                + a[0 + 1 * 4] * (a[1 + 2 * 4] * a[2 + 0 * 4] - a[1 + 0 * 4] * a[2 + 2 * 4])
                + a[0 + 2 * 4] * (a[1 + 0 * 4] * a[2 + 1 * 4] - a[1 + 1 * 4] * a[2 + 0 * 4])
            ) / det;

        return b;
    }

    public static double[] r8mat_l_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_L_INVERSE inverts a lower triangular R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    A lower triangular matrix is a matrix whose only nonzero entries
        //    occur on or below the diagonal.
        //
        //    The inverse of a lower triangular matrix is a lower triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6.
        //
        //  Parameters:
        //
        //    Input, int N, number of rows and columns in the matrix.
        //
        //    Input, double A[N*N], the lower triangular matrix.
        //
        //    Output, double R8MAT_L_INVERSE[N*N], the inverse matrix.
        //
    {
        int j;

        double[] b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i < j)
                {
                    b[i + j * n] = 0.0;
                }
                else if (j == i)
                {
                    b[i + j * n] = 1.0 / a[i + j * n];
                }
                else
                {
                    double temp = 0.0;
                    int k;
                    for (k = 0; k < i; k++)
                    {
                        temp += a[i + k * n] * b[k + j * n];
                    }

                    b[i + j * n] = -temp / a[i + i * n];
                }
            }
        }

        return b;
    }
       
}