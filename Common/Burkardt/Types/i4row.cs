using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int i4row_compare(int m, int n, int[] a, int i, int j)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_COMPARE compares two rows of a integer array.
            //
            //  Discussion:
            //
            //    The two dimensional information is stored in a one dimensional array,
            //    by columns.  The entry A(I,J) is stored in A[I+J*M].
            //
            //    The input arguments I and J are row indices.  They DO NOT use the
            //    C convention of starting at 0, but rather start at 1.
            //
            //  Example:
            //
            //    Input:
            //
            //      M = 3, N = 4, I = 2, J = 3
            //
            //      A = (
            //        1  2  3  4
            //        5  6  7  8
            //        9 10 11 12 )
            //
            //    Output:
            //
            //      I4ROW_COMPARE = -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, int  A[M*N], the array of data.
            //
            //    Input, int I, J, the rows to be compared.
            //    I and J must be between 1 and M.
            //
            //    Output, int I4ROW_COMPARE, the results of the comparison:
            //    -1, row I < row J,
            //     0, row I = row J,
            //    +1, row J < row I.
            //
        {
            int k;
            //
            //  Check that I and J are legal.
            //
            if (i < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_COMPARE - Fatal error!");
                Console.WriteLine("  Row index I is less than 1.");
                Console.WriteLine("  I = " + i + "");
                return (1);
            }
            else if (m < i)
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_COMPARE - Fatal error!");
                Console.WriteLine("  Row index I is out of bounds.");
                Console.WriteLine("  I = " + i + "");
                Console.WriteLine("  Maximum legal value is M = " + m + "");
                return (1);
            }

            if (j < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_COMPARE - Fatal error!");
                Console.WriteLine("  Row index J is less than 1.");
                Console.WriteLine("  J = " + j + "");
                return (1);
            }
            else if (m < j)
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_COMPARE - Fatal error!");
                Console.WriteLine("  Row index J is out of bounds.");
                Console.WriteLine("  J = " + j + "");
                Console.WriteLine("  Maximum legal value is M = " + m + "");
                return (1);
            }

            if (i == j)
            {
                return 0;
            }

            for (k = 0; k < n; k++)
            {
                if (a[(i - 1) + k * m] < a[(j - 1) + k * m])
                {
                    return -1;
                }
                else if (a[(j - 1) + k * m] < a[(i - 1) + k * m])
                {
                    return +1;
                }
            }

            return 0;
        }

        public static int[] i4row_max(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MAX returns the maximums of an I4ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MAX =
            //      3
            //      7
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, int I4ROW_AMAX[M], the maximums of the rows.
            //
        {
            int i;
            int j;
            int[] amax = new int[m];

            for (i = 0; i < m; i++)
            {
                amax[i] = a[i + 0 * m];

                for (j = 1; j < n; j++)
                {
                    if (amax[i] < a[i + j * m])
                    {
                        amax[i] = a[i + j * m];
                    }
                }
            }

            return amax;
        }

        public static double[] i4row_mean(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MEAN returns the means of an I4ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MEAN =
            //      2
            //      5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, double I4ROW_MEAN[M], the means, or averages, of the rows.
            //
        {
            int i;
            int j;
            double[] mean = new double[m];

            for (i = 0; i < m; i++)
            {
                mean[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean[i] = mean[i] + (double) a[i + j * m];
                }

                mean[i] = mean[i] / (double) (n);
            }

            return mean;
        }

        public static int[] i4row_min(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_MIN returns the minimums of an I4ROW.
            //
            //  Example:
            //
            //    A =
            //      1  2  3
            //      2  6  7
            //
            //    MIN =
            //      1
            //      2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array to be examined.
            //
            //    Output, int I4ROW_AMIN[M], the minimums of the rows.
            //
        {
            int i;
            int j;
            int[] amin = new int[m];

            for (i = 0; i < m; i++)
            {
                amin[i] = a[i + 0 * m];
                for (j = 1; j < n; j++)
                {
                    if (a[i + j * m] < amin[i])
                    {
                        amin[i] = a[i + j * m];
                    }
                }
            }

            return amin;
        }

        public static void i4row_sort_a ( int m, int n, ref int[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_SORT_A ascending sorts the rows of an I4ROW.
            //
            //  Definition:
            //
            //    In lexicographic order, the statement "X < Y", applied to two
            //    vectors X and Y of length M, means that there is some index I, with
            //    1 <= I <= M, with the property that
            //
            //      X(J) = Y(J) for J < I,
            //    and
            //      X(I) < Y(I).
            //
            //    In other words, X is less than Y if, at the first index where they
            //    differ, the X value is less than the Y value.
            //
            //  Example:
            //
            //    Input:
            //
            //      M = 5, N = 3
            //
            //      A =
            //        3  2  1
            //        2  4  3
            //        3  1  8
            //        2  4  2
            //        1  9  9
            //
            //    Output:
            //
            //      A =
            //        1  9  9
            //        2  4  2
            //        2  4  3
            //        3  1  8
            //        3  2  1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of A.
            //
            //    Input, int N, the number of columns of A.
            //
            //    Input/output, int A[M*N].
            //    On input, the array of M rows of N-vectors.
            //    On output, the rows of A have been sorted in ascending
            //    lexicographic order.
            //
        {
            int i;
            int indx;
            int isgn;
            int j;
            //
            //  Initialize.
            //
            i = 0;
            indx = 0;
            isgn = 0;
            j = 0;
            //
            //  Call the external heap sorter.
            //
            SortHeapExternalData data = new SortHeapExternalData();
            for ( ; ; )
            {
                Helpers.sort_heap_external (ref data, m, ref indx, ref i, ref j, isgn );
                //
                //  Interchange the I and J objects.
                //
                if ( 0 < indx )
                {
                    i4row_swap ( m, n, ref a, i, j );
                }
                //
                //  Compare the I and J objects.
                //
                else if ( indx < 0 )
                {
                    isgn = i4row_compare ( m, n, a, i, j );
                }
                else if ( indx == 0 )
                {
                    break;
                }

            }
        }
        
        public static void i4row_swap ( int m, int n, ref int[] a, int irow1, int irow2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4ROW_SWAP swaps two rows of an I4ROW.
        //
        //  Discussion:
        //
        //    The two dimensional information is stored as a one dimensional
        //    array, by columns.
        //
        //    The row indices are 1 based, NOT 0 based!  However, a preprocessor
        //    variable, called OFFSET, can be reset from 1 to 0 if you wish to
        //    use 0-based indices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, int A[M*N], an array of data.
        //
        //    Input, int IROW1, IROW2, the two rows to swap.
        //    These indices should be between 1 and M.
        //
        {
            int OFFSET = 1;

            int j;
            int t;
            //
            //  Check.
            //
            if ( irow1 < 0+OFFSET || m-1+OFFSET < irow1 )
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_SWAP - Fatal error!");
                Console.WriteLine("  IROW1 is out of range.");
                return;
            }

            if ( irow2 < 0+OFFSET || m-1+OFFSET < irow2 )
            {
                Console.WriteLine("");
                Console.WriteLine("I4ROW_SWAP - Fatal error!");
                Console.WriteLine("  IROW2 is out of range.");
                return;
            }

            if ( irow1 == irow2 )
            {
                return;
            }
            for ( j = 0; j < n; j++ )
            {
                t                   = a[irow1-OFFSET+j*m];
                a[irow1-OFFSET+j*m] = a[irow2-OFFSET+j*m];
                a[irow2-OFFSET+j*m] = t;
            }
        }
        
        public static double[] i4row_variance(int m, int n, int[] a)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROW_VARIANCE returns the variances of an I4ROW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 October 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the array.
            //
            //    Input, int A[M*N], the array whose variances are desired.
            //
            //    Output, double I4ROW_VARIANCE[M], the variances of the rows.
            //
        {
            int i;
            int j;
            double mean;
            double[] variance = new double[m];

            for (i = 0; i < m; i++)
            {
                mean = 0.0;
                for (j = 0; j < n; j++)
                {
                    mean = mean + (double) a[i + j * m];
                }

                mean = mean / (double) (n);

                variance[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    variance[i] = variance[i] + Math.Pow(((double) a[i + j * m] - mean), 2);
                }

                if (1 < n)
                {
                    variance[i] = variance[i] / (double) (n - 1);
                }
                else
                {
                    variance[i] = 0.0;
                }
            }

            return variance;
        }
        
        public static int[] i4rows_copy_new ( int m, int n, int[] a1 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4MROWS_COPY_NEW copies an I4ROWS to a "new" I4MAT.
            //
            //  Discussion:
            //
            //    An I4ROWS is an MxN array of I4's, stored by (I,J) -> [I*N+J].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, int A1[M*N], the matrix to be copied.
            //
            //    Output, int I4ROWS_COPY_NEW[M*N], the copy of A1.
            //
        {
            int[] a2;
            int i;
            int j;

            a2 = new int[m*n];

            for ( i = 0; i < m; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    a2[i*n+j] = a1[i*n+j];
                }
            }
            return a2;
        }
        
        public static int[] i4rows_zeros_new ( int m, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4ROWS_ZEROS_NEW returns a new zeroed I4ROWS.
            //
            //  Discussion:
            //
            //    An I4ROWS is a doubly dimensioned array of I4 values, stored as a vector
            //    in row-major order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 April 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns.
            //
            //    Output, int I4ROWS_ZEROS_NEW[M*N], the new zeroed matrix.
            //
        {
            int[] a;
            int i;
            int j;

            a = new int[m*n];

            for ( i = 0; i < m; i++ )
            {
                for ( j = 0; j < n; j++ )
                {
                    a[i*n+j] = 0;
                }
            }
            return a;
        }
        
    }
}