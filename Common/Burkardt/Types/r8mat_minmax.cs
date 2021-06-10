using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8mat_amax ( int m, int n, double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_AMAX returns the maximum absolute value entry of an R8MAT.
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
            //    21 April 2012
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
            //    Output, double R8MAT_AMAX, the maximum absolute value entry of A.
            //
        {
            int i;
            int j;
            double value;

            value = Math.Abs ( a[0+0*m] );

            for ( j = 0; j < n; j++ )
            {
                for ( i = 0; i < m; i++ )
                {
                    value = Math.Max ( value, Math.Abs ( a[i+j*m] ) );
                }
            }
            return value;
        }
        
        public static double r8mat_max(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAX returns the maximum entry of an R8MAT.
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
            //    21 May 2011
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
            //    Output, double R8MAT_MAX, the maximum entry of A.
            //
        {
            int i;
            int j;
            double value;

            value = a[0 + 0 * m];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (value < a[i + j * m])
                    {
                        value = a[i + j * m];
                    }
                }
            }

            return value;
        }

        public static double[] r8mat_max_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAX_COLUMNS returns the column maximums of an R8MAT.
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
            //    Output, double R8MAT_MAX_COLUMNS[N], the column maximums.
            //
        {
            int i;
            int j;
            double[] max_columns;
            double value;

            max_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                value = a[0 + j * m];
                for (i = 1; i < m; i++)
                {
                    if (value < a[i + j * m])
                    {
                        value = a[i + j * m];
                    }
                }

                max_columns[j] = value;
            }

            return max_columns;
        }

        public static void r8mat_max_index(int m, int n, double[] a, ref int i_max, ref int j_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.
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
            //    18 September 2005
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
            //    Output, int &I_MAX, &J_MAX, the indices of the maximum entry of A.
            //
        {
            int i;
            int i2;
            int j;
            int j2;

            i2 = -1;
            j2 = -1;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (i2 == -1 && j2 == -1)
                    {
                        i2 = i;
                        j2 = j;
                    }
                    else if (a[i2 + j2 * m] < a[i + j * m])
                    {
                        i2 = i;
                        j2 = j;
                    }
                }
            }

            i_max = i2 + 1;
            j_max = j2 + 1;

        }

        public static double[] r8mat_max_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAX_ROWS returns the row maximums of an R8MAT.
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
            //    Output, double R8MAT_MAX_ROWS[M], the row maximums.
            //
        {
            int i;
            int j;
            double[] max_rows;
            double value;

            max_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                value = a[i + 0 * m];
                for (j = 1; j < n; j++)
                {
                    if (value < a[i + j * m])
                    {
                        value = a[i + j * m];
                    }
                }

                max_rows[i] = value;
            }

            return max_rows;
        }

        public static double r8mat_maxcol_minrow(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    R8MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )
            //
            //    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 October 2005
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
            //    Input, double A[M*N], the matrix.
            //
            //    Output, double R8MAT_MAXCOL_MINROW, the maximum column
            //    minimum row entry of A.
            //
        {
            int i;
            int j;
            double minrow;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            value = -r8_huge;

            for (i = 0; i < m; i++)
            {
                minrow = r8_huge;

                for (j = 0; j < n; j++)
                {
                    minrow = Math.Min(minrow, a[i + j * m]);
                }

                value = Math.Max(value, minrow);
            }

            return value;
        }

        public static double r8mat_maxrow_mincol(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    R8MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )
            //
            //    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 October 2005
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
            //    Input, double A[M*N], the matrix.
            //
            //    Output, double R8MAT_MAXROW_MINCOL, the maximum row
            //    minimum column entry of A.
            //
        {
            int i;
            int j;
            double mincol;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            value = -r8_huge;

            for (j = 0; j < n; j++)
            {
                mincol = r8_huge;
                for (i = 0; i < m; i++)
                {
                    mincol = Math.Min(mincol, a[i + j * m]);
                }

                value = Math.Max(value, mincol);
            }

            return value;
        }

        public static double r8mat_mean(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MEAN returns the mean of an R8MAT.
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
            //    03 September 2013
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
            //    Output, double R8MAT_MEAN, the mean of A.
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

            value = value / (double) (m * n);

            return value;
        }

        public static double[] r8mat_mean_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MEAN_COLUMNS returns the column means of an R8MAT.
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
            //    Output, double R8MAT_MEAN_COLUMNS[N], the column means.
            //
        {
            int i;
            int j;
            double[] mean_columns;
            double value;

            mean_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                value = 0.0;
                for (i = 0; i < m; i++)
                {
                    value = value + a[i + j * m];
                }

                mean_columns[j] = value / (double) m;
            }

            return mean_columns;
        }

        public static double[] r8mat_mean_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MEAN_ROWS returns the row means of an R8MAT.
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
            //    Output, double R8MAT_MEAN_ROWS[M], the row means.
            //
        {
            int i;
            int j;
            double[] mean_rows;
            double value;

            mean_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                value = 0.0;
                for (j = 0; j < n; j++)
                {
                    value = value + a[i + j * m];
                }

                mean_rows[i] = value / (double) n;
            }

            return mean_rows;
        }

        public static double r8mat_min(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MIN returns the minimum entry of an R8MAT.
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
            //    21 May 2011
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
            //    Output, double R8MAT_MIN, the minimum entry of A.
            //
        {
            int i;
            int j;
            double value;

            value = a[0 + 0 * m];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (a[i + j * m] < value)
                    {
                        value = a[i + j * m];
                    }
                }
            }

            return value;
        }

        public static double[] r8mat_min_columns(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MIN_COLUMNS returns the column minimums of an R8MAT.
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
            //    Output, double R8MAT_MIN_COLUMNS[N], the column minimums.
            //
        {
            int i;
            int j;
            double[] min_columns;
            double value;

            min_columns = new double[n];

            for (j = 0; j < n; j++)
            {
                value = a[0 + j * m];
                for (i = 1; i < m; i++)
                {
                    if (a[i + j * m] < value)
                    {
                        value = a[i + j * m];
                    }
                }

                min_columns[j] = value;
            }

            return min_columns;
        }

        public static void r8mat_min_index(int m, int n, double[] a, ref int i_min, ref int j_min)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.
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
            //    18 September 2005
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
            //    Output, int &I_MIN, &J_MIN, the indices of the minimum entry of A.
            //
        {
            int i;
            int i2;
            int j;
            int j2;

            i2 = -1;
            j2 = -1;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    if (i2 == -1 && j2 == -1)
                    {
                        i2 = i;
                        j2 = j;
                    }
                    else if (a[i + j * m] < a[i2 + j2 * m])
                    {
                        i2 = i;
                        j2 = j;
                    }
                }
            }

            i_min = i2 + 1;
            j_min = j2 + 1;
        }

        public static double[] r8mat_min_rows(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MIN_ROWS returns the row minimums of an R8MAT.
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
            //    Output, double R8MAT_MIN_ROWS[M], the row minimums.
            //
        {
            int i;
            int j;
            double[] min_rows;
            double value;

            min_rows = new double[m];

            for (i = 0; i < m; i++)
            {
                value = a[i + 0 * m];
                for (j = 1; j < n; j++)
                {
                    if (a[i + j * m] < value)
                    {
                        value = a[i + j * m];
                    }
                }

                min_rows[i] = value;
            }

            return min_rows;
        }

        public static double r8mat_mincol_maxrow(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    R8MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )
            //
            //    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 October 2005
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
            //    Input, double A(M,N), the matrix.
            //
            //    Output, double R8MAT_MINCOL_MAXROW, the minimum column
            //    maximum row entry of A.
            //
        {
            int i;
            int j;
            double maxrow;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            value = r8_huge;

            for (i = 0; i < m; i++)
            {
                maxrow = -r8_huge;
                for (j = 0; j < n; j++)
                {
                    maxrow = Math.Max(maxrow, a[i + j * m]);
                }

                value = Math.Min(value, maxrow);
            }

            return value;
        }

        public static double r8mat_minrow_maxcol(int m, int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N matrix.
            //
            //  Discussion:
            //
            //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
            //    in column-major order.
            //
            //    R8MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )
            //
            //    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 October 2005
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
            //    Input, double A[M*N], the matrix.
            //
            //    Output, double R8MAT_MINROW_MAXCOL, the minimum row
            //    maximum column entry of A.
            //
        {
            int i;
            int j;
            double maxcol;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            value = r8_huge;

            for (j = 0; j < n; j++)
            {
                maxcol = -r8_huge;
                for (i = 0; i < m; i++)
                {
                    maxcol = Math.Max(maxcol, a[i + j * m]);
                }

                value = Math.Min(value, maxcol);
            }

            return value;
        }

        public static void r8mat_minvm(int n1, int n2, double[] a, double[] b, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MINVM computes inverse(A) * B for R8MAT's.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8 values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N1, N2, the order of the matrices.
            //
            //    Input, double A[N1*N1], B[N1*N2], the matrices.
            //
            //    Output, double C[N1*N2], the result, C = inverse(A) * B.
            //
        {
            double[] alu;
            double[] d;

            alu = r8mat_copy_new(n1, n1, a);

            d = r8mat_fss_new(n1, ref alu, n2, b);

            r8mat_copy(n1, n2, d, ref c);
        }

        public static double[] r8mat_minvm_new(int n1, int n2, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_MINVM_NEW returns inverse(A) * B for R8MAT's.
            //
            //  Discussion:
            //
            //    An R8MAT is an array of R8 values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N1, N2, the order of the matrices.
            //
            //    Input, double A[N1*N1], B[N1*N2], the matrices.
            //
            //    Output, double R8MAT_MINVM_NEW[N1*N2], the result, C = inverse(A) * B.
            //
        {
            double[] alu;
            double[] c;

            alu = r8mat_copy_new(n1, n1, a);
            c = r8mat_fss_new(n1, ref alu, n2, b);

            return c;
        }

        
    }
}