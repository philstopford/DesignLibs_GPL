using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8mat_givens_post(int n, double[] a, int row, int col)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Givens post-multiplier matrix G(ROW,COL) has the property that
            //    the (ROW,COL)-th entry of A*G is zero.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices A and G.
            //
            //    Input, double A[N*N], the matrix to be operated upon.
            //
            //    Input, int ROW, COL, the row and column of the
            //    entry of A*G which is to be zeroed out.
            //
            //    Output, double R8MAT_GIVENS_POST[N*N], the Givens rotation matrix.
            //    G is an orthogonal matrix, that is, the inverse of
            //    G is the transpose of G.
            //
        {
            double[] g;
            double theta;

            g = r8mat_identity_new(n);

            theta = Math.Atan2(a[row - 1 + (col - 1) * n], a[row - 1 + (row - 1) * n]);

            g[row - 1 + (row - 1) * n] = Math.Cos(theta);
            g[row - 1 + (col - 1) * n] = -Math.Sin(theta);
            g[col - 1 + (row - 1) * n] = Math.Sin(theta);
            g[col - 1 + (col - 1) * n] = Math.Cos(theta);

            return g;
        }

        public static double[] r8mat_givens_pre(int n, double[] a, int row, int col)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    The Givens premultiplier rotation matrix G(ROW,COL) has the
            //    property that the (ROW,COL)-th entry of G*A is zero.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices A and G.
            //
            //    Input, double A[N*N], the matrix to be operated upon.
            //
            //    Input, int ROW, COL, the row and column of the
            //    entry of the G*A which is to be zeroed out.
            //
            //    Output, double R8MAT_GIVENS_PRE[N*N], the Givens rotation matrix.
            //    G is an orthogonal matrix, that is, the inverse of
            //    G is the transpose of G.
            //
        {
            double[] g;
            double theta;

            g = r8mat_identity_new(n);

            theta = Math.Atan2(a[row - 1 + (col - 1) * n], a[col - 1 + (col - 1) * n]);

            g[row - 1 + (row - 1) * n] = Math.Cos(theta);
            g[row - 1 + (col - 1) * n] = -Math.Sin(theta);
            g[col - 1 + (row - 1) * n] = Math.Sin(theta);
            g[col - 1 + (col - 1) * n] = Math.Cos(theta);

            return g;
        }
        
    }
}