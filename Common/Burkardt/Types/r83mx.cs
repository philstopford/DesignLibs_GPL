namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r83_mxv_new(int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_MXV_NEW multiplies a R83 matrix times a vector.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83_MXV_NEW[N], the product A * x.
        //
        {
            double[] b;
            int i;

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = a[1 + i * 3] * x[i];
            }

            for (i = 0; i < n - 1; i++)
            {
                b[i] = b[i] + a[0 + (i + 1) * 3] * x[i + 1];
            }

            for (i = 1; i < n; i++)
            {
                b[i] = b[i] + a[2 + (i - 1) * 3] * x[i - 1];
            }

            return b;
        }
    }
}