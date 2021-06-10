namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8mat_sum(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_SUM returns the sum of an R8MAT.
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
            //    03 January 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, double A[M*N], the array.
            //
            //    Output, double R8MAT_SUM, the sum of the entries.
            //
        {
            int i;
            int j;
            double value;

            value = 0.0;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    value = value + a[i + j * m];
                }
            }

            return value;
        }

        public static double[] r8mat_sum_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_SUM_COLUMNS returns the column sums of an R8MAT.
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
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_SUM_COLUMNS[N], the column sums.
            //
        {
            int i;
            int j;
            double[] sum_columns;
            double value;

            sum_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                value = 0.0;
                for (i = 0; i < m; i++)
                {
                    value = value + a[i + j * m];
                }

                sum_columns[j] = value;
            }

            return sum_columns;
        }

        public static double[] r8mat_sum_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_SUM_ROWS returns the row sums of an R8MAT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows in A.
            //
            //    Input, int N, the number of columns in A.
            //
            //    Input, double A[M*N], the M by N matrix.
            //
            //    Output, double R8MAT_SUM_ROWS[M], the row sums.
            //
        {
            int i;
            int j;
            double[] sum_rows;
            double value;

            sum_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                value = 0.0;
                for (j = 0; j < n; j++)
                {
                    value = value + a[i + j * m];
                }

                sum_rows[i] = value;
            }

            return sum_rows;
        }
        
    }
}