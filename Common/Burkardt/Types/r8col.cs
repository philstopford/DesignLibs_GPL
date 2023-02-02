using System;
using System.Globalization;
using Burkardt.SortNS;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int r8col_compare ( int m, int n, double[] a, int i, int j )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_COMPARE compares two columns in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 3, N = 4, I = 2, J = 4
        //
        //      A = (
        //        1.  2.  3.  4.
        //        5.  6.  7.  8.
        //        9. 10. 11. 12. )
        //
        //    Output:
        //
        //      R8COL_COMPARE = -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the M by N array.
        //
        //    Input, int I, J, the columns to be compared.
        //    I and J must be between 1 and N.
        //
        //    Output, int R8COL_COMPARE, the results of the comparison:
        //    -1, column I < column J,
        //     0, column I = column J,
        //    +1, column J < column I.
        //
    {
        //
        //  Check.
        //
        if (i < 1 || n < i)
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_COMPARE - Fatal error!");
            Console.WriteLine("  Column index I is out of bounds.");
            Console.WriteLine("  I = " + i + "");
            return 1;
        }

        if (j < 1 || n < j)
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_COMPARE - Fatal error!");
            Console.WriteLine("  Column index J is out of bounds.");
            Console.WriteLine("  J = " + j + "");
            return 1;
        }

        int value = 0;

        if (i == j)
        {
            return value;
        }

        int k = 0;

        while (k < m)
        {
            if (a[k + (i - 1) * m] < a[k + (j - 1) * m])
            {
                value = -1;
                return value;
            }

            if (a[k + (j - 1) * m] < a[k + (i - 1) * m])
            {
                value = +1;
                return value;
            }

            k += 1;
        }

        return value;
    }

    public static double[] r8col_duplicates(int m, int n, int n_unique, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_DUPLICATES generates an R8COL with some duplicate columns.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    This routine generates a random R8COL with a specified number of
        //    duplicate columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in each column of A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, int N_UNIQUE, the number of unique columns in A.
        //    1 <= N_UNIQUE <= N.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double R8COL_DUPLICATES[M*N], the array.
        //
    {
        int i;
        int j1;
        int j2;

        if (n_unique < 1 || n < n_unique)
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_DUPLICATES - Fatal error!");
            Console.WriteLine("  1 <= N_UNIQUE <= N is required.");
            return null;
        }

        double[] a = UniformRNG.r8col_uniform_01_new(m, n_unique, ref seed);
        //
        //  Randomly copy unique columns.
        //
        for (j1 = n_unique; j1 < n; j1++)
        {
            j2 = UniformRNG.i4_uniform_ab(0, n_unique - 1, ref seed);
            for (i = 0; i < m; i++)
            {
                a[i + j1 * m] = a[i + j2 * m];
            }
        }

        //
        //  Permute the columns.
        //
        for (j1 = 0; j1 < n; j1++)
        {
            j2 = UniformRNG.i4_uniform_ab(j1, n - 1, ref seed);
            for (i = 0; i < m; i++)
            {
                (a[i + j1 * m], a[i + j2 * m]) = (a[i + j2 * m], a[i + j1 * m]);
            }
        }

        return a;
    }

    public static int r8col_find(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_FIND seeks a column value in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 3,
        //      N = 4,
        //
        //      A = (
        //        1.  2.  3.  4.
        //        5.  6.  7.  8.
        //        9. 10. 11. 12. )
        //
        //      x = ( 3.,
        //            7.,
        //           11. )
        //
        //    Output:
        //
        //      R8COL_FIND = 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], a table of numbers, regarded as
        //    N columns of vectors of length M.
        //
        //    Input, double X[M], a vector to be matched with a column of A.
        //
        //    Output, int R8COL_FIND, the (one-based) index of the first column of A
        //    which exactly matches every entry of X, or -1 if no match
        //    could be found.
        //
    {
        int j;

        int col = -1;

        for (j = 1; j <= n; j++)
        {
            col = j;

            int i;
            for (i = 1; i <= m; i++)
            {
                if (!(Math.Abs(x[i - 1] - a[i - 1 + (j - 1) * m]) > typeMethods.r8_epsilon()))
                {
                    continue;
                }

                col = -1;
                break;
            }

            if (col != -1)
            {
                return col;
            }
        }

        return col;
    }

    public static int[] r8col_first_index(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_FIRST_INDEX indexes the first occurrence of values in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
        //    the first column whose entries are equal to A(1:M,J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of A.
        //    The length of an "element" of A, and the number of "elements".
        //
        //    Input, double A[M*N], the array.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8COL_FIRST_INDEX[N], the first occurrence index.
        //
    {
        int j1;

        int[] first_index = new int[n];

        for (j1 = 0; j1 < n; j1++)
        {
            first_index[j1] = -1;
        }

        for (j1 = 0; j1 < n; j1++)
        {
            switch (first_index[j1])
            {
                case -1:
                {
                    first_index[j1] = j1;

                    int j2;
                    for (j2 = j1 + 1; j2 < n; j2++)
                    {
                        double diff = 0.0;
                        int i;
                        for (i = 0; i < m; i++)
                        {
                            diff = Math.Max(diff, Math.Abs(a[i + j1 * m] - a[i + j2 * m]));
                        }

                        if (diff <= tol)
                        {
                            first_index[j2] = j1;
                        }
                    }

                    break;
                }
            }
        }

        return first_index;
    }

    public static void r8col_flip(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_FLIP flips the entries in each column of an R8COL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the array to be flipped.
        //
    {
        int j;

        int ihi = m / 2;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < ihi; i++)
            {
                double t = a[i + j * m];
                a[i + j * m] = a[m + 1 - i + j * m];
                a[m - 1 - j + j * m] = t;
            }
        }
    }

    public static double[] r8col_indicator_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_INDICATOR_NEW sets up an "indicator" R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The value of each entry suggests its location, as in:
        //
        //      11  12  13  14
        //      21  22  23  24
        //      31  32  33  34
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8COL_INDICATOR_NEW[M*N], the table.
        //
    {
        int i;

        double[] a = new double[m * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return a;
    }

    public static int r8col_insert(int n_max, int m, ref int n, ref double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_INSERT inserts a column into an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    Input:
        //
        //      N_MAX = 10,
        //      M = 3,
        //      N = 4,
        //
        //      A = (
        //        1.  2.  3.  4.
        //        5.  6.  7.  8.
        //        9. 10. 11. 12. )
        //
        //      X = ( 3., 4., 18. )
        //
        //    Output:
        //
        //      N = 5,
        //
        //      A = (
        //        1.  2.  3.  3.  4.
        //        5.  6.  4.  7.  8.
        //        9. 10. 18. 11. 12. )
        //
        //      R8COL_INSERT = 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N_MAX, the maximum number of columns in A.
        //
        //    Input, int M, the number of rows.
        //
        //    Input/output, int N, the number of columns.
        //    If the new column is inserted into the table, then the output
        //    value of N will be increased by 1.
        //
        //    Input/output, double A[M*N_MAX], a table of numbers, regarded
        //    as an array of columns.  The columns must have been sorted
        //    lexicographically.
        //
        //    Input, double X[M], a vector of data which will be inserted
        //    into the table if it does not already occur.
        //
        //    Output, int R8COL_INSERT.
        //    I, X was inserted into column I.
        //    -I, column I was already equal to X.
        //    0, N = N_MAX.
        //
    {
        int col;
        int i;
        int j;
        //
        //  Refuse to work if N_MAX <= N.
        //
        if (n_max <= n)
        {
            col = 0;
            return col;
        }

        //
        //  Stick X temporarily in column N+1, just so it's easy to use R8COL_COMPARE.
        //
        for (i = 0; i < m; i++)
        {
            a[i + n * m] = x[i];
        }

        //
        //  Do a binary search.
        //
        int low = 1;
        int high = n;

        for (;;)
        {
            if (high < low)
            {
                col = low;
                break;
            }

            int mid = (low + high) / 2;

            int isgn = r8col_compare(m, n + 1, a, mid, n + 1);

            switch (isgn)
            {
                case 0:
                    col = -mid;
                    return col;
                case -1:
                    low = mid + 1;
                    break;
                case +1:
                    high = mid - 1;
                    break;
            }
        }

        //
        //  Shift part of the table up to make room.
        //
        for (j = n - 1; col - 1 <= j; j--)
        {
            for (i = 0; i < m; i++)
            {
                a[i + (j + 1) * m] = a[i + j * m];
            }
        }

        //
        //  Insert the new column.
        //
        for (i = 0; i < m; i++)
        {
            a[i + (col - 1) * m] = x[i];
        }

        n += 1;

        return col;
    }

    public static double[] r8col_max(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MAX returns the column maximums of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8COL_MAX[N], the maximums of the columns.
        //
    {
        int j;

        double[] amax = new double[n];

        for (j = 0; j < n; j++)
        {
            amax[j] = a[0 + j * m];
            int i;
            for (i = 0; i < m; i++)
            {
                amax[j] = Math.Max(amax[j], a[i + j * m]);
            }
        }

        return amax;
    }

    public static int[] r8col_max_index(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MAX_INDEX returns the indices of column maximums in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, int R8COL_MAX_INDEX[N]; entry I is the row of A in which
        //    the maximum for column I occurs.
        //
    {
        int j;

        int[] imax = new int[n];

        for (j = 0; j < n; j++)
        {
            imax[j] = 1;
            double amax = a[0 + j * m];

            int i;
            for (i = 1; i < m; i++)
            {
                if (!(amax < a[i + j * m]))
                {
                    continue;
                }

                imax[j] = i + 1;
                amax = a[i + j * m];
            }
        }

        return imax;
    }

    public static void r8col_max_one(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MAX_ONE rescales an R8COL so each column maximum is 1.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the array to be rescaled.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i_big = 0;
            int i;
            for (i = 1; i < m; i++)
            {
                if (Math.Abs(a[i_big + j * m]) < Math.Abs(a[i + j * m]))
                {
                    i_big = i;
                }
            }

            double temp = a[i_big + j * m];

            if (temp == 0.0)
            {
                continue;
            }

            for (i = 0; i < m; i++)
            {
                a[i + j * m] /= temp;
            }
        }
    }

    public static void r8col_sort_heap_a ( int m, int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORT_HEAP_A ascending heapsorts an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    In lexicographic order, the statement "X < Y", applied to two real
        //    vectors X and Y of length M, means that there is some index I, with
        //    1 <= I <= M, with the property that
        //
        //      X(J) = Y(J) for J < I,
        //    and
        //      X(I) < Y(I).
        //
        //    In other words, the first time they differ, X is smaller.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N].
        //    On input, the array of N columns of M-vectors.
        //    On output, the columns of A have been sorted in lexicographic order.
        //
    {
        SortHeapExternalData data = new();

        switch (m)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 1:
                return;
        }
        //
        //  Initialize.
        //
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;
        //
        //  Call the external heap sorter.
        //
        for ( ; ; )
        {
            Sort.sort_heap_external ( ref data, n, ref indx, ref i, ref j, isgn );
            //
            //  Interchange the I and J objects.
            //
            if ( 0 < indx )
            {
                r8col_swap ( m, n, a, i, j );
            }
            //
            //  Compare the I and J objects.
            //
            else if ( indx < 0 )
            {
                isgn = r8col_compare ( m, n, a, i, j );
            }
            else
            {
                break;
            }
        }
    }

    public static int[] r8col_sort_heap_index_a(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      A(*,INDX(*)) is sorted,
        //
        //    Note that the index vector is 0-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in each column of A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the array.
        //
        //    Output, int R8COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
        //    I-th column of the sorted array is A(*,INDX(I)).
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                return null;
        }

        int[] indx = new int[n];

        for (i = 0; i < n; i++)
        {
            indx[i] = i;
        }

        switch (n)
        {
            case 1:
                indx[0] = indx[0];
                return indx;
        }

        double[] column = new double[m];

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            int indxt;
            int k;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                for (k = 0; k < m; k++)
                {
                    column[k] = a[k + indxt * m];
                }
            }
            else
            {
                indxt = indx[ir - 1];
                for (k = 0; k < m; k++)
                {
                    column[k] = a[k + indxt * m];
                }

                indx[ir - 1] = indx[0];
                ir -= 1;

                if (ir == 1)
                {
                    indx[0] = indxt;
                    break;
                }
            }

            i = l;
            int j = l + l;

            while (j <= ir)
            {
                int isgn;
                if (j < ir)
                {
                    isgn = r8vec_compare(m, a, a, aIndex: + indx[j - 1] * m, bIndex: + indx[j] * m);

                    switch (isgn)
                    {
                        case < 0:
                            j += 1;
                            break;
                    }
                }

                isgn = r8vec_compare(m, column, a, bIndex: + indx[j - 1] * m);

                switch (isgn)
                {
                    case < 0:
                        indx[i - 1] = indx[j - 1];
                        i = j;
                        j += j;
                        break;
                    default:
                        j = ir + 1;
                        break;
                }
            }

            indx[i - 1] = indxt;
        }

        return indx;
    }

    public static void r8col_sort_quick_a(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORT_QUICK_A ascending quick sorts an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the row order of A, and the length of a column.
        //
        //    Input, int N, the number of columns of A.
        //
        //    Input/output, double A[M*N].
        //    On input, the array to be sorted.
        //    On output, the array has been sorted.
        //
    {
        const int LEVEL_MAX = 30;

        int l_segment = 0;
        int[] rsave = new int[LEVEL_MAX];
        int r_segment = 0;

        switch (m)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8COL_SORT_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                Console.WriteLine("  N = " + n + "");
                return;
            case 1:
                return;
        }

        int level = 1;
        rsave[level - 1] = n + 1;
        int base_ = 1;
        int n_segment = n;

        for (;;)
        {
            //
            //  Partition the segment.
            //
            r8col_part_quick_a(m, n_segment, ref a, ref l_segment, ref r_segment, aIndex: + (base_ - 1) * m);
            switch (l_segment)
            {
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                case > 1 when LEVEL_MAX < level:
                    Console.WriteLine("");
                    Console.WriteLine("R8COL_SORT_QUICK_A - Fatal error!");
                    Console.WriteLine("  Exceeding recursion maximum of " + LEVEL_MAX + "");
                    return;
                case > 1:
                    level += 1;
                    n_segment = l_segment;
                    rsave[level - 1] = r_segment + base_ - 1;
                    break;
                //
                default:
                {
                    if (r_segment < n_segment)
                    {
                        n_segment = n_segment + 1 - r_segment;
                        base_ = base_ + r_segment - 1;
                    }
                    //
                    //  Otherwise, we back up a level if there is an earlier one.
                    //
                    else
                    {
                        for (;;)
                        {
                            switch (level)
                            {
                                case <= 1:
                                    return;
                            }

                            base_ = rsave[level - 1];
                            n_segment = rsave[level - 2] - rsave[level - 1];
                            level -= 1;

                            if (0 < n_segment)
                            {
                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }
    }

    public static void r8col_sorted_tol_undex ( int m, int n, double[] a, int unique_num,
            double tol, ref int[] undx, ref int[] xdnu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_TOL_UNDEX: index tolerably unique entries of a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points, to the unique elements of A, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of A, the index of
        //    the unique sorted element of A.
        //
        //    This is all done with index vectors, so that the elements of
        //    A are never moved.
        //
        //    Assuming A is already sorted, we examine the entries of A in order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector A could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of X in sorted order, using the formula
        //
        //      XU(*) = A(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector A, or
        //    any element of it, by index, as follows:
        //
        //      A(I) = XU(XDNU(I)).
        //
        //    We could then replace A by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of A, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector A, the unique sort and inverse unique
        //    sort vectors and the compressed unique sorted vector.
        //
        //      I      A      XU  Undx  Xdnu
        //    ----+------+------+-----+-----+
        //      0 | 11.0 |  11.0    0     0
        //      1 | 11.0 |  22.0    4     0
        //      2 | 11.0 |  33.0    7     0
        //      3 | 11.0 |  55.0    8     0
        //      4 | 22.0 |                1
        //      5 | 22.0 |                1
        //      6 | 22.0 |                1
        //      7 | 33.0 |                2
        //      8 | 55.0 |                3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values,
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, int UNIQUE_NUM, the number of unique values in A.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[N], the XDNU vector.
        //
    {
        //
        //  Consider entry I = 0.
        //  It is unique, so set the number of unique items to K.
        //  Set the K-th unique item to I.
        //  Set the representative of item I to the K-th unique item.
        //
        int i = 0;
        int k = 0;
        undx[k] = i;
        xdnu[i] = k;
        //
        //  Consider entry I.
        //
        //  If it is unique, increase the unique count K, set the
        //  K-th unique item to I, and set the representative of I to K.
        //
        //  If it is not unique, set the representative of item I to a
        //  previously determined unique item that is close to it.
        //
        for ( i = 1; i < n; i++ )
        {
            bool unique = true;

            int j;
            for ( j = 0; j <= k; j++ )
            {
                int i2 = undx[j];
                double diff = 0.0;
                int i3;
                for ( i3 = 0; i3 < m; i3++ )
                {
                    diff = Math.Max ( diff, Math.Abs ( a[i3+i*m] - a[i3+i2*m] ) );
                }

                if (!(diff <= tol))
                {
                    continue;
                }

                unique = false;
                xdnu[i] = j;
                break;
            }
            switch (unique)
            {
                case true:
                    k += 1;
                    undx[k] = i;
                    xdnu[i] = k;
                    break;
            }
        }
    }

    public static int r8col_sorted_tol_unique(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_TOL_UNIQUE keeps tolerably unique elements in a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The columns of the array can be ascending or descending sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A(M,N).
        //    On input, the sorted array of N columns of M-vectors.
        //    On output, a sorted array of columns of M-vectors.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8COL_SORTED_TOL_UNIQUE, the number of unique columns.
        //
    {
        int i;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }

        unique_num = 1;

        for (i = 1; i < n; i++)
        {
            bool unique = true;
            int k;
            int j;
            for (j = 0; j < unique_num; j++)
            {
                double diff = 0.0;
                for (k = 0; k < m; k++)
                {
                    diff = Math.Max(diff, Math.Abs(a[k + i * m] - a[k + j * m]));
                }

                if (!(diff < tol))
                {
                    continue;
                }

                unique = false;
                break;
            }

            switch (unique)
            {
                case true:
                {
                    for (k = 0; k < m; k++)
                    {
                        a[k + unique_num * m] = a[k + i * m];
                    }

                    unique_num += 1;
                    break;
                }
            }
        }

        return unique_num;
    }

    public static int r8col_sorted_tol_unique_count(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_TOL_UNIQUE_COUNT counts tolerably unique elements in a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The columns of the array may be ascending or descending sorted.
        //
        //    If the tolerance is large enough, then the concept of uniqueness
        //    can become ambiguous.  If we have a tolerance of 1.5, then in the
        //    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
        //    one unique entry?  That would be because 1 may be regarded as unique,
        //    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
        //    be unique and so on.
        //
        //    This seems wrongheaded.  So I prefer the idea that an item is not
        //    unique under a tolerance only if it is close to something that IS unique.
        //    Thus, the unique items are guaranteed to cover the space if we include
        //    a disk of radius TOL around each one.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], a sorted array, containing
        //    N columns of data.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8COL_SORTED_UNIQUE_COUNT, the number of unique columns.
        //
    {
        int[] undx = new int[n];
        //
        //  Consider entry I = 0.
        //  It is unique, so set the number of unique items to K.
        //  Set the K-th unique item to I.
        //  Set the representative of item I to the K-th unique item.
        //
        int i = 0;
        int k = 0;
        undx[k] = i;
        //
        //  Consider entry I.
        //
        //  If it is unique, increase the unique count K, set the
        //  K-th unique item to I, and set the representative of I to K.
        //
        //  If it is not unique, set the representative of item I to a
        //  previously determined unique item that is close to it.
        //
        for (i = 1; i < n; i++)
        {
            bool unique = true;

            int j;
            for (j = 0; j <= k; j++)
            {
                int i2 = undx[j];
                double diff = 0.0;
                int i3;
                for (i3 = 0; i3 < m; i3++)
                {
                    diff = Math.Max(diff, Math.Abs(a[i3 + i * m] - a[i3 + i2 * m]));
                }

                if (!(diff <= tol))
                {
                    continue;
                }

                unique = false;
                break;
            }

            switch (unique)
            {
                case true:
                    k += 1;
                    undx[k] = i;
                    break;
            }
        }

        k += 1;

        return k;
    }

    public static void r8col_sorted_undex(int m, int n, double[] a, int unique_num,
            int[] undx, int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_UNDEX returns unique sorted indexes for a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points, to the unique elements of A, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of A, the index of
        //    the unique sorted element of A.
        //
        //    This is all done with index vectors, so that the elements of
        //    A are never moved.
        //
        //    Assuming A is already sorted, we examine the entries of A in order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector A could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of X in sorted order, using the formula
        //
        //      XU(*) = A(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector A, or
        //    any element of it, by index, as follows:
        //
        //      A(I) = XU(XDNU(I)).
        //
        //    We could then replace A by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of A, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector A, the unique sort and inverse unique
        //    sort vectors and the compressed unique sorted vector.
        //
        //      I      A      XU  Undx  Xdnu
        //    ----+------+------+-----+-----+
        //      0 | 11.0 |  11.0    0     0
        //      1 | 11.0 |  22.0    4     0
        //      2 | 11.0 |  33.0    7     0
        //      3 | 11.0 |  55.0    8     0
        //      4 | 22.0 |                1
        //      5 | 22.0 |                1
        //      6 | 22.0 |                1
        //      7 | 33.0 |                2
        //      8 | 55.0 |                3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values,
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, int UNIQUE_NUM, the number of unique values in A.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[N], the XDNU vector.
        //
    {
        //
        //  Walk through the sorted array.
        //
        int i = 0;

        int j = 0;
        undx[j] = i;

        xdnu[i] = j;

        for (i = 1; i < n; i++)
        {
            double diff = 0.0;
            int k;
            for (k = 0; k < m; k++)
            {
                diff = Math.Max(diff, Math.Abs(a[k + i * m] - a[k + undx[j] * m]));
            }

            switch (diff)
            {
                case > 0.0:
                    j += 1;
                    undx[j] = i;
                    break;
            }

            xdnu[i] = j;
        }
    }

    public static int r8col_sorted_unique(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_UNIQUE keeps unique elements in a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The columns of the array can be ascending or descending sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A(M,N).
        //    On input, the sorted array of N columns of M-vectors.
        //    On output, a sorted array of columns of M-vectors.
        //
        //    Output, int UNIQUE_NUM, the number of unique columns.
        //
    {
        int j2;
        int unique_num;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return unique_num;
        }

        int j1 = 0;

        for (j2 = 1; j2 < n; j2++)
        {
            bool equal = true;
            int i;
            for (i = 0; i < m; i++)
            {
                if (!(Math.Abs(a[i + j1 * m] - a[i + j2 * m]) > typeMethods.r8_epsilon()))
                {
                    continue;
                }

                equal = false;
                break;
            }

            switch (equal)
            {
                case false:
                {
                    j1 += 1;
                    for (i = 0; i < m; i++)
                    {
                        a[i + j1 * m] = a[i + j2 * m];
                    }

                    break;
                }
            }
        }

        unique_num = j1 + 1;

        return unique_num;
    }

    public static int r8col_sorted_unique_count(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The columns of the array may be ascending or descending sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], a sorted array, containing
        //    N columns of data.
        //
        //    Output, int R8COL_SORTED_UNIQUE_COUNT, the number of unique columns.
        //
    {
        int j2;

        int unique_num = 0;

        switch (n)
        {
            case <= 0:
                return unique_num;
        }

        unique_num = 1;
        int j1 = 0;

        for (j2 = 1; j2 < n; j2++)
        {
            bool equal = true;
            int i;
            for (i = 0; i < m; i++)
            {
                if (!(Math.Abs(a[i + j1 * m] - a[i + j2 * m]) > typeMethods.r8_epsilon()))
                {
                    continue;
                }

                equal = false;
                break;
            }

            switch (equal)
            {
                case false:
                    unique_num += 1;
                    j1 = j2;
                    break;
            }
        }

        return unique_num;
    }

    public static void r8col_sortr_a(int m, int n, ref double[] a, int key)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SORTR_A ascending sorts one column of an R8COL, adjusting all entries.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N].
        //    On input, an unsorted M by N array.
        //    On output, rows of the array have been shifted in such
        //    a way that column KEY of the array is in nondecreasing order.
        //
        //    Input, int KEY, the column in which the "key" value
        //    is stored.  On output, column KEY of the array will be
        //    in nondecreasing order.
        //
    {
        SortHeapExternalData data = new();

        switch (m)
        {
            case <= 0:
                return;
        }

        if (key < 1 || n < key)
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_SORTR_A - Fatal error!");
            Console.WriteLine("  The value of KEY is not a legal column index.");
            Console.WriteLine("  KEY = " + key + "");
            Console.WriteLine("  N = " + n + "");
            return;
        }

        //
        //  Initialize.
        //
        int i = 0;
        int indx = 0;
        int isgn = 0;
        int j = 0;
        //
        //  Call the external heap sorter.
        //
        for (;;)
        {
            Sort.sort_heap_external(ref data, m, ref indx, ref i, ref j, isgn);
            //
            //  Interchange the I and J objects.
            //
            if (0 < indx)
            {
                int k;
                for (k = 0; k < n; k++)
                {
                    (a[i - 1 + k * m], a[j - 1 + k * m]) = (a[j - 1 + k * m], a[i - 1 + k * m]);
                }
            }
            //
            //  Compare the I and J objects.
            //
            else if (indx < 0)
            {
                if (a[i - 1 + (key - 1) * m] < a[j - 1 + (key - 1) * m])
                {
                    isgn = -1;
                }
                else
                {
                    isgn = +1;
                }
            }
            else
            {
                break;
            }
        }

    }

    public static double[] r8col_sum(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SUM sums the columns of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8COL_SUM[N], the sums of the columns.
        //
    {
        int j;

        double[] colsum = new double[n];

        for (j = 0; j < n; j++)
        {
            colsum[j] = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                colsum[j] += a[i + j * m];
            }
        }

        return colsum;
    }

    public static void r8col_swap(int m, int n, double[] a, int j1, int j2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SWAP swaps columns J1 and J2 of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 3, N = 4, J1 = 2, J2 = 4
        //
        //      A = (
        //        1.  2.  3.  4.
        //        5.  6.  7.  8.
        //        9. 10. 11. 12. )
        //
        //    Output:
        //
        //      A = (
        //        1.  4.  3.  2.
        //        5.  8.  7.  6.
        //        9. 12. 11. 10. )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the M by N array.
        //
        //    Input, int J1, J2, the columns to be swapped.
        //    These columns are 1-based.
        //
    {
        int i;

        if (j1 < 1 || n < j1 || j2 < 1 || n < j2)
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_SWAP - Fatal error!");
            Console.WriteLine("  J1 or J2 is out of bounds.");
            Console.WriteLine("  J1 =   " + j1 + "");
            Console.WriteLine("  J2 =   " + j2 + "");
            Console.WriteLine("  NCOL = " + n + "");
            return;
        }

        if (j1 == j2)
        {
            return;
        }

        for (i = 0; i < m; i++)
        {
            (a[i + (j1 - 1) * m], a[i + (j2 - 1) * m]) = (a[i + (j2 - 1) * m], a[i + (j1 - 1) * m]);
        }

    }

    public static double[] r8col_to_r8vec(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TO_R8VEC converts an R8COL to an R8VEC.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    This routine is not really useful in our C++ implementation, since
        //    we actually store an M by N matrix exactly as a vector already.
        //
        //  Example:
        //
        //    M = 3, N = 4
        //
        //    A =
        //      11 12 13 14
        //      21 22 23 24
        //      31 32 33 34
        //
        //    R8COL_TO_R8VEC = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the M by N array.
        //
        //    Output, double X[M*N], a vector containing the N columns of A.
        //
    {
        int j;

        double[] x = new double[m * n];

        int k = 0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                x[k] = a[i + j * m];
                k += 1;
            }
        }

        return x;
    }

    public static void r8col_tol_undex(int m, int n, double[] a, int unique_num, double tol,
            ref int[] undx, ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TOL_UNDEX indexes tolerably unique entries of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points to the unique elements of A, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of A, the index of
        //    the unique sorted element of A.
        //
        //    This is all done with index vectors, so that the elements of
        //    A are never moved.
        //
        //    The first step of the algorithm requires the indexed sorting
        //    of A, which creates arrays INDX and XDNI.  (If all the entries
        //    of A are unique, then these arrays are the same as UNDX and XDNU.)
        //
        //    We then use INDX to examine the entries of A in sorted order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector A could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of A in sorted order, using the formula
        //
        //      XU(*) = A(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector A, or
        //    any element of it, by index, as follows:
        //
        //      A(I) = XU(XDNU(I)).
        //
        //    We could then replace A by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of A, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector A, the sort and inverse sort
        //    index vectors, and the unique sort and inverse unique sort vectors
        //    and the compressed unique sorted vector.
        //
        //      I     A  Indx  Xdni       XU  Undx  Xdnu
        //    ----+-----+-----+-----+--------+-----+-----+
        //      0 | 11.     0     0 |    11.     0     0
        //      1 | 22.     2     4 |    22.     1     1
        //      2 | 11.     5     1 |    33.     3     0
        //      3 | 33.     8     7 |    55.     4     2
        //      4 | 55.     1     8 |                  3
        //      5 | 11.     6     2 |                  0
        //      6 | 22.     7     5 |                  1
        //      7 | 22.     3     6 |                  1
        //      8 | 11.     4     3 |                  0
        //
        //    INDX(2) = 3 means that sorted item(2) is A(3).
        //    XDNI(2) = 5 means that A(2) is sorted item(5).
        //
        //    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
        //    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
        //
        //    XU(XDNU(I))) = X(I).
        //    XU(I)        = X(UNDX(I)).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values,
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, int UNIQUE_NUM, the number of unique values in A.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[N], the XDNU vector.
        //
    {
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8col_sort_heap_index_a(m, n, a);
        //
        //  Consider entry I = 0.
        //  It is unique, so set the number of unique items to K.
        //  Set the K-th unique item to I.
        //  Set the representative of item I to the K-th unique item.
        //
        int i = 0;
        int k = 0;
        undx[k] = indx[i];
        xdnu[indx[i]] = k;
        //
        //  Consider entry I.
        //
        //  If it is unique, increase the unique count K, set the
        //  K-th unique item to I, and set the representative of I to K.
        //
        //  If it is not unique, set the representative of item I to a
        //  previously determined unique item that is close to it.
        //
        for (i = 1; i < n; i++)
        {
            bool unique = true;
            int j;
            for (j = 0; j <= k; j++)
            {
                double diff = 0.0;
                int i2;
                for (i2 = 0; i2 < m; i2++)
                {
                    diff = Math.Max(diff, Math.Abs(a[i2 + indx[i] * m] - a[i2 + undx[j] * m]));
                }

                if (!(diff <= tol))
                {
                    continue;
                }

                unique = false;
                xdnu[indx[i]] = j;
                break;
            }

            switch (unique)
            {
                case true:
                    k += 1;
                    undx[k] = indx[i];
                    xdnu[indx[i]] = k;
                    break;
            }
        }
    }

    public static int r8col_tol_unique_count(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The columns of the array may be ascending or descending sorted.
        //
        //    If the tolerance is large enough, then the concept of uniqueness
        //    can become ambiguous.  If we have a tolerance of 1.5, then in the
        //    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
        //    one unique entry?  That would be because 1 may be regarded as unique,
        //    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
        //    be unique and so on.
        //
        //    This seems wrongheaded.  So I prefer the idea that an item is not
        //    unique under a tolerance only if it is close to something that IS unique.
        //    Thus, the unique items are guaranteed to cover the space if we include
        //    a disk of radius TOL around each one.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array of N columns of data.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8COL_TOL_UNIQUE_COUNT, the number of unique columns.
        //
    {
        int[] undx = new int[n];
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8col_sort_heap_index_a(m, n, a);
        //
        //  Consider entry I = 0.
        //  It is unique, so set the number of unique items to K.
        //  Set the K-th unique item to I.
        //  Set the representative of item I to the K-th unique item.
        //
        int i = 0;
        int k = 0;
        undx[k] = indx[i];
        //
        //  Consider entry I.
        //
        //  If it is unique, increase the unique count K, set the
        //  K-th unique item to I, and set the representative of I to K.
        //
        //  If it is not unique, set the representative of item I to a
        //  previously determined unique item that is close to it.
        //
        for (i = 1; i < n; i++)
        {
            bool unique = true;
            int j;
            for (j = 0; j <= k; j++)
            {
                double diff = 0.0;
                int i2;
                for (i2 = 0; i2 < m; i2++)
                {
                    diff = Math.Max(diff, Math.Abs(a[i2 + indx[i] * m] - a[i2 + undx[j] * m]));
                }

                if (!(diff <= tol))
                {
                    continue;
                }

                unique = false;
                break;
            }

            switch (unique)
            {
                case true:
                    k += 1;
                    undx[k] = indx[i];
                    break;
            }
        }

        k += 1;

        return k;
    }

    public static int[] r8col_tol_unique_index(int m, int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TOL_UNIQUE_INDEX indexes tolerably unique entries in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
        //    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
        //    gathered in order, then
        //
        //      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of A.
        //
        //    Input, double A[M*N], the array.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int R8COL_TOL_UNIQUE_INDEX[N], the unique index.
        //
    {
        int j1;

        int[] unique_index = new int[n];

        for (j1 = 0; j1 < n; j1++)
        {
            unique_index[j1] = -1;
        }

        int unique_num = 0;

        for (j1 = 0; j1 < n; j1++)
        {
            switch (unique_index[j1])
            {
                case -1:
                {
                    unique_index[j1] = unique_num;

                    int j2;
                    for (j2 = j1 + 1; j2 < n; j2++)
                    {
                        double diff = 0.0;
                        int i;
                        for (i = 0; i < m; i++)
                        {
                            diff = Math.Max(diff, Math.Abs(a[i + j1 * m] - a[i + j2 * m]));
                        }

                        if (diff <= tol)
                        {
                            unique_index[j2] = unique_num;
                        }
                    }

                    unique_num += 1;
                    break;
                }
            }
        }

        return unique_index;
    }

    public static void r8col_transpose_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TRANSPOSE_PRINT prints an R8MAT, transposed.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8col_transpose_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8col_transpose_print_some(int m, int n, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, int ILO, JLO, the first row and column to print.
        //
        //    Input, int IHI, JHI, the last row and column to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int i2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        int i2lo_lo = ilo switch
        {
            < 1 => 1,
            _ => ilo
        };

        int i2lo_hi = ihi < m ? m : ihi;

        for (i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo += INCX)
        {
            int i2hi = i2lo + INCX - 1;

            if (m < i2hi)
            {
                i2hi = m;
            }

            if (ihi < i2hi)
            {
                i2hi = ihi;
            }

            int inc = i2hi + 1 - i2lo;

            Console.WriteLine("");
            string cout = "  Row: ";
            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout += (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Col");
            Console.WriteLine("");

            int j2lo = jlo switch
            {
                < 1 => 1,
                _ => jlo
            };

            int j2hi = n < jhi ? n : jhi;

            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout = (j - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ":";
                int i2;
                for (i2 = 1; i2 <= inc; i2++)
                {
                    i = i2lo - 1 + i2;
                    cout += a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8col_undex(int m, int n, double[] a, int unique_num, ref int[] undx,
            ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_UNDEX indexes unique entries in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points to the unique elements of A, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of A, the index of
        //    the unique sorted element of A.
        //
        //    This is all done with index vectors, so that the elements of
        //    A are never moved.
        //
        //    The first step of the algorithm requires the indexed sorting
        //    of A, which creates arrays INDX and XDNI.  (If all the entries
        //    of A are unique, then these arrays are the same as UNDX and XDNU.)
        //
        //    We then use INDX to examine the entries of A in sorted order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector A could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of A in sorted order, using the formula
        //
        //      XU(*) = A(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector A, or
        //    any element of it, by index, as follows:
        //
        //      A(I) = XU(XDNU(I)).
        //
        //    We could then replace A by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of A, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector A, the sort and inverse sort
        //    index vectors, and the unique sort and inverse unique sort vectors
        //    and the compressed unique sorted vector.
        //
        //      I     A  Indx  Xdni       XU  Undx  Xdnu
        //    ----+-----+-----+-----+--------+-----+-----+
        //      0 | 11.     0     0 |    11.     0     0
        //      1 | 22.     2     4 |    22.     1     1
        //      2 | 11.     5     1 |    33.     3     0
        //      3 | 33.     8     7 |    55.     4     2
        //      4 | 55.     1     8 |                  3
        //      5 | 11.     6     2 |                  0
        //      6 | 22.     7     5 |                  1
        //      7 | 22.     3     6 |                  1
        //      8 | 11.     4     3 |                  0
        //
        //    INDX(2) = 3 means that sorted item(2) is A(3).
        //    XDNI(2) = 5 means that A(2) is sorted item(5).
        //
        //    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
        //    XDNU(8) = 2 means that A(8) is at unique sorted item(2).
        //
        //    XU(XDNU(I))) = A(I).
        //    XU(I)        = A(UNDX(I)).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the data values.
        //
        //    Input, int N, the number of data values,
        //
        //    Input, double A[M*N], the data values.
        //
        //    Input, int UNIQUE_NUM, the number of unique values in A.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Output, int UNDX[UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[N], the XDNU vector.
        //
    {
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8col_sort_heap_index_a(m, n, a);
        //
        //  Walk through the implicitly sorted array.
        //
        int i = 0;

        int j = 0;
        undx[j] = indx[i];

        xdnu[indx[i]] = j;

        for (i = 1; i < n; i++)
        {
            double diff = 0.0;
            int k;
            for (k = 0; k < m; k++)
            {
                diff = Math.Max(diff, Math.Abs(a[k + indx[i] * m] - a[k + undx[j] * m]));
            }

            switch (diff)
            {
                case > 0.0:
                    j += 1;
                    undx[j] = indx[i];
                    break;
            }

            xdnu[indx[i]] = j;
        }
    }

    public static double[] r8col_variance(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_VARIANCE returns the variances of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, double A[M*N], the array whose variances are desired.
        //
        //    Output, double R8COL_VARIANCE[N], the variances of the rows.
        //
    {
        double[] variance = new double[n];

        for (int j = 0; j < n; j++)
        {
            double mean = 0.0;
            for (int i = 0; i < m; i++)
            {
                mean += a[i + j * m];
            }

            mean /= m;

            variance[j] = 0.0;
            for (int i = 0; i < m; i++)
            {
                variance[j] += Math.Pow(a[i + j * m] - mean, 2);
            }

            switch (m)
            {
                case > 1:
                    variance[j] /= m - 1;
                    break;
                default:
                    variance[j] = 0.0;
                    break;
            }
        }

        return variance;
    }

    public static double[] r8col_mean(int m, int n, double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MEAN returns the column means of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Example:
        //
        //    A =
        //      1  2  3
        //      2  6  7
        //
        //    R8COL_MEAN =
        //      1.5  4.0  5.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
        //
    {
        double[] mean = new double[n];

        for (int j = 0; j < n; j++)
        {
            mean[j] = 0.0;
            for (int i = 0; i < m; i++)
            {
                mean[j] += a[i + j * m];
            }

            mean[j] /= m;
        }

        return mean;
    }

    public static double[] r8col_min ( int m, int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MIN returns the column minimums of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, double R8COL_MIN[N], the minimums of the columns.
        //
    {
        int j;

        double[] amin = new double[n];

        for ( j = 0; j < n; j++ )
        {
            amin[j] = a[0+j*m];
            int i;
            for ( i = 0; i < m; i++ )
            {
                amin[j] = Math.Min ( amin[j], a[i+j*m] );
            }
        }

        return amin;
    }

    public static int[] r8col_min_index(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_MIN_INDEX returns the indices of column minimums in an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], the array to be examined.
        //
        //    Output, int R8COL_MIN_INDEX[N]; entry I is the row of A in which
        //    the minimum for column I occurs.
        //
    {
        int j;

        int[] imin = new int[n];

        for (j = 0; j < n; j++)
        {
            imin[j] = 1;
            double amin = a[0 + j * m];

            int i;
            for (i = 1; i < m; i++)
            {
                if (!(a[i + j * m] < amin))
                {
                    continue;
                }

                imin[j] = i + 1;
                amin = a[i + j * m];
            }
        }

        return imin;
    }

    public static void r8col_normalize_li(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_NORMALIZE_LI normalizes an R8COL with the column infinity norm.
        //
        //  Discussion:
        //
        //    Each column is scaled so that the entry of maximum norm has the value 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the array to be normalized.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            double c = a[0 + j * m];
            int i;
            for (i = 1; i < m; i++)
            {
                if (Math.Abs(c) < Math.Abs(a[i + j * m]))
                {
                    c = a[i + j * m];
                }
            }

            if (c == 0.0)
            {
                continue;
            }

            for (i = 0; i < m; i++)
            {
                a[i + m * j] /= c;
            }
        }
    }

    public static void r8col_part_quick_a(int m, int n, ref double[] a, ref int l, ref int r, int aIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PART_QUICK_A reorders the columns of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    The routine reorders the columns of A.  Using A(1:M,1) as a
        //    key, all entries of A that are less than or equal to the key will
        //    precede the key, which precedes all entries that are greater than the key.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 2, N = 8
        //      A = ( 2  8  6  0 10 10  0  5
        //            4  8  2  2  6  0  6  8 )
        //
        //    Output:
        //
        //      L = 2, R = 4
        //
        //      A = (  0  0  2  8  6 10 10  4
        //             2  6  4  8  2  6  0  8 )
        //             ----     -------------
        //             LEFT KEY     RIGHT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the row dimension of A, and the length of a column.
        //
        //    Input, int N, the column dimension of A.
        //
        //    Input/output, double A[M*N].  On input, the array to be checked.
        //    On output, A has been reordered as described above.
        //
        //    Output, int &L, &R, the indices of A that define the three segments.
        //    Let KEY = the input value of A(1:M,1).  Then
        //    I <= L                 A(1:M,I) < KEY;
        //         L < I < R         A(1:M,I) = KEY;
        //                 R <= I    KEY < A(1:M,I).
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8COL_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case 1:
                l = 0;
                r = 2;
                return;
        }

        double[] key = new double[m];

        for (i = 0; i < m; i++)
        {
            key[i] = a[(i + 0 * m + aIndex ) % a.Length];
        }

        int k = 1;
        //
        //  The elements of unknown size have indices between L+1 and R-1.
        //
        l = 1;
        r = n + 1;

        for (j = 1; j < n; j++)
        {
            if (r8vec_gt(m, a, key, a1Index: +l * m))
            {
                r -= 1;
                r8vec_swap(m, ref a, ref a, startIndexA1: +(r - 1) * m + aIndex, startIndexA2: +l * m + aIndex);
            }
            else if (r8vec_eq(m, a, key, startIndexA1: +l * m + aIndex))
            {
                k += 1;
                r8vec_swap(m, ref a, ref a, startIndexA1: +(k - 1) * m + aIndex, startIndexA2: +l * m + aIndex);
                l += 1;
            }
            else if (r8vec_lt(m, a, key, startIndexA1: +l * m + aIndex))
            {
                l += 1;
            }
        }

        //
        //  Shift small elements to the left.
        //
        for (j = 0; j < l - k; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[(i + j * m + aIndex ) % a.Length] = a[(i + (j + k) * m + aIndex ) % a.Length];
            }
        }

        //
        //  Shift KEY elements to center.
        //
        for (j = l - k; j < l; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[(i + j * m + aIndex ) % a.Length] = key[i];
            }
        }

        //
        //  Update L.
        //
        l -= k;

    }

    public static void r8col_permute ( int m, int n, int[] p, int base_, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PERMUTE permutes an R8COL in place.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    This routine permutes an array of real "objects", but the same
        //    logic can be used to permute an array of objects of any arithmetic
        //    type, or an array of objects of any complexity.  The only temporary
        //    storage required is enough to store a single object.  The number
        //    of data movements made is N + the number of cycles of order 2 or more,
        //    which is never more than N + N/2.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 2
        //      N = 5
        //      P = (   2,    4,    5,    1,    3 )
        //      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
        //          (11.0, 22.0, 33.0, 44.0, 55.0 )
        //      BASE = 1
        //
        //    Output:
        //
        //      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
        //             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the length of objects.
        //
        //    Input, int N, the number of objects.
        //
        //    Input, int P[N], the permutation.  P(I) = J means
        //    that the I-th element of the output array should be the J-th
        //    element of the input array.
        //
        //    Input, int BASE, is 0 for a 0-based permutation and 1 for a
        //    1-based permutation.
        //
        //    Input/output, double A[M*N], the array to be permuted.
        //
    {
        int i;
        int istart;
        int j;

        if (!perm_check(n, p, base_))
        {
            Console.WriteLine("");
            Console.WriteLine("R8COL_PERMUTE - Fatal error!");
            Console.WriteLine("  PERM_CHECK rejects this permutation.");
            return;
        }

        //
        //  In order for the sign negation trick to work, we need to assume that the
        //  entries of P are strictly positive.  Presumably, the lowest number is BASE.
        //  So temporarily add 1-BASE to each entry to force positivity.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = p[i] + 1 - base_;
        }

        double[] a_temp = new double[m];
        //
        //  Search for the next element of the permutation that has not been used.
        //
        for (istart = 1; istart <= n; istart++)
        {
            switch (p[istart - 1])
            {
                case < 0:
                    break;
                default:
                {
                    if (p[istart - 1] == istart)
                    {
                        p[istart - 1] = -p[istart - 1];
                    }
                    else
                    {
                        for (i = 0; i < m; i++)
                        {
                            a_temp[i] = a[i + (istart - 1) * m];
                        }

                        int iget = istart;
                        //
                        //  Copy the new value into the vacated entry.
                        //
                        for (;;)
                        {
                            int iput = iget;
                            iget = p[iget - 1];

                            p[iput - 1] = -p[iput - 1];

                            if (iget < 1 || n < iget)
                            {
                                Console.WriteLine("");
                                Console.WriteLine("R8COL_PERMUTE - Fatal error!");
                                Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                                Console.WriteLine("  an illegal value IGET = " + iget + ".");
                                return;
                            }

                            if (iget == istart)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    a[i + (iput - 1) * m] = a_temp[i];
                                }

                                break;
                            }

                            for (i = 0; i < m; i++)
                            {
                                a[i + (iput - 1) * m] = a[i + (iget - 1) * m];
                            }
                        }
                    }

                    break;
                }
            }
        }

        //
        //  Restore the signs of the entries.
        //
        for (j = 0; j < n; j++)
        {
            p[j] = -p[j];
        }

        //
        //  Restore the base of the entries.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = p[i] - 1 + base_;
        }
    }

    public static void r8col_permute(int m, int n, int[] p, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PERMUTE permutes an R8COL in place.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of double precision values, regarded
        //    as an array of N columns of length M.
        //
        //    This routine permutes an array of real "objects", but the same
        //    logic can be used to permute an array of objects of any arithmetic
        //    type, or an array of objects of any complexity.  The only temporary
        //    storage required is enough to store a single object.  The number
        //    of data movements made is N + the number of cycles of order 2 or more,
        //    which is never more than N + N/2.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 2
        //      N = 5
        //      P = (   2,    4,    5,    1,    3 )
        //      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
        //          (11.0, 22.0, 33.0, 44.0, 55.0 )
        //
        //    Output:
        //
        //      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
        //             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the length of objects.
        //
        //    Input, int N, the number of objects.
        //
        //    Input/output, double A[M*N], the array to be permuted.
        //
        //    Input, int P[N], the permutation.  P(I) = J means
        //    that the I-th element of the output array should be the J-th
        //    element of the input array.  P must be a legal permutation
        //    of the integers from 1 to N, otherwise the algorithm will
        //    fail catastrophically.
        //
    {
        int i;
        int istart;
        int j;

        double[] a_temp = new double[m];
        //
        //  Need to increment the entries by 1 in order to use the sign trick.
        //
        for (i = 0; i < n; i++)
        {
            p[i] += 1;
        }

        //
        //  Search for the next element of the permutation that has not been used.
        //
        for (istart = 1; istart <= n; istart++)
        {
            switch (p[istart - 1])
            {
                case < 0:
                    break;
                default:
                {
                    if (p[istart - 1] == istart)
                    {
                        p[istart - 1] = -p[istart - 1];
                    }
                    else
                    {
                        for (i = 0; i < m; i++)
                        {
                            a_temp[i] = a[i + (istart - 1) * m];
                        }

                        int iget = istart;
                        //
                        //  Copy the new value into the vacated entry.
                        //
                        for (;;)
                        {
                            int iput = iget;
                            iget = p[iget - 1];

                            p[iput - 1] = -p[iput - 1];

                            if (iget < 1 || n < iget)
                            {
                                Console.WriteLine();
                                Console.WriteLine("R8COL_PERMUTE - Fatal error!");
                                Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                                Console.WriteLine("  an illegal value IGET = " + iget + "");
                                return;
                            }

                            if (iget == istart)
                            {
                                for (i = 0; i < m; i++)
                                {
                                    a[i + (iput - 1) * m] = a_temp[i];
                                }

                                break;
                            }

                            for (i = 0; i < m; i++)
                            {
                                a[i + (iput - 1) * m] = a[i + (iget - 1) * m];
                            }
                        }
                    }

                    break;
                }
            }
        }

        //
        //  Restore the signs of the entries.
        //
        for (j = 0; j < n; j++)
        {
            p[j] = -p[j];
        }

        //
        //  Need to unincrement the entries by 1.
        //
        for (i = 0; i < n; i++)
        {
            p[i] -= 1;
        }
    }

    public static void r8col_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PRINT prints an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
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
        //    Input, string TITLE, a title.
        //
    {
        r8col_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8col_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_PRINT_SOME prints some of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            if (n < j2hi)
            {
                j2hi = n;
            }

            if (jhi < j2hi)
            {
                j2hi = jhi;
            }

            Console.WriteLine("");
            //
            //  For each column J in the current range...
            //
            //  Write the header.
            //
            string cout = "  Col:    ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += (j - 1).ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("");
            int i2lo = ilo switch
            {
                //
                //  Determine the range of the rows in this strip.
                //
                > 1 => ilo,
                _ => 1
            };

            int i2hi = ihi < m ? ihi : m;

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = (i - 1).ToString(CultureInfo.InvariantCulture).PadLeft(5) + ": ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    cout += a[i - 1 + (j - 1) * m].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void r8col_reverse(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_REVERSE reverses the order of the columns of an R8COL.
        //
        //  Discussion:
        //
        //    An R8COL is an M by N array of R8's, regarded as an array of N columns,
        //    each of length M.
        //
        //    To reverse the columns is to start with something like
        //
        //      11 12 13 14 15
        //      21 22 23 24 25
        //      31 32 33 34 35
        //      41 42 43 44 45
        //      51 52 53 54 55
        //
        //    and return
        //
        //      15 14 13 12 11
        //      25 24 23 22 21
        //      35 34 33 32 31
        //      45 44 43 42 41
        //      55 54 53 52 51
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input/output, double A[M*N], the matrix whose columns are to be flipped.
        //
    {
        int i;

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n / 2; j++)
            {
                (a[i + j * m], a[i + (n - 1 - j) * m]) = (a[i + (n - 1 - j) * m], a[i + j * m]);
            }
        }
    }

    public static void r8col_separation(int m, int n, double[] a, ref double d_min, ref double d_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8COL_SEPARATION returns the "separation" of an R8COL.
        //
        //  Discussion:
        //
        //    D_MIN is the minimum distance between two columns,
        //    D_MAX is the maximum distance between two columns.
        //
        //    The distances are measured using the Loo norm.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns 
        //    in the array.  If N < 2, it does not make sense to call this routine.
        //
        //    Input, double A[M*N], the array whose variances are desired.
        //
        //    Output, double &D_MIN, &D_MAX, the minimum and maximum distances.
        //
    {
        int j1;

        d_min = r8_huge();
        d_max = 0.0;

        for (j1 = 0; j1 < n; j1++)
        {
            int j2;
            for (j2 = j1 + 1; j2 < n; j2++)
            {
                double d = 0.0;
                int i;
                for (i = 0; i < m; i++)
                {
                    d = Math.Max(d, Math.Abs(a[i + j1 * m] - a[i + j2 * m]));
                }

                d_min = Math.Min(d_min, d);
                d_max = Math.Max(d_max, d);
            }
        }
    }
}