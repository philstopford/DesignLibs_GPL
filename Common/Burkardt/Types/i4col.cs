using System;
using Burkardt.SortNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int i4col_compare(int m, int n, int[] a, int i, int j )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_COMPARE compares columns I and J of an I4COL.
        //
        //  Example:
        //
        //    Input:
        //
        //      M = 3, N = 4, I = 2, J = 4
        //
        //      A = (
        //        1  2  3  4
        //        5  6  7  8
        //        9 10 11 12 )
        //
        //    Output:
        //
        //      I4COL_COMPARE = -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int A[M*N], an array of N columns of vectors of length M.
        //
        //    Input, int I, J, the columns to be compared.
        //    I and J must be between 1 and N.
        //
        //    Output, int I4COL_COMPARE, the results of the comparison:
        //    -1, column I < column J,
        //     0, column I = column J,
        //    +1, column J < column I.
        //
    {
        switch (i)
        {
            //
            //  Check.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("I4COL_COMPARE - Fatal error!");
                Console.WriteLine("  Column index I = " + i + " is less than 1.");
                return 1;
        }

        if (n < i)
        {
            Console.WriteLine("");
            Console.WriteLine("I4COL_COMPARE - Fatal error!");
            Console.WriteLine("  N = " + n + " is less than column index I = " + i + ".");
            return 1;
        }

        switch (j)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("I4COL_COMPARE - Fatal error!");
                Console.WriteLine("  Column index J = " + j + " is less than 1.");
                return 1;
        }

        if (n < j)
        {
            Console.WriteLine("");
            Console.WriteLine("I4COL_COMPARE - Fatal error!");
            Console.WriteLine("  N = " + n + " is less than column index J = " + j + ".");
            return 1;
        }

        if (i == j)
        {
            return 0;
        }

        int k = 1;

        while (k <= m)
        {
            if (a[k - 1 + (i - 1) * m] < a[k - 1 + (j - 1) * m])
            {
                return -1;
            }

            if (a[k - 1 + (j - 1) * m] < a[k - 1 + (i - 1) * m])
            {
                return 1;
            }

            k += 1;
        }

        return 0;
    }

    public static void i4col_find_item(int m, int n, int[] table, int item,
            ref int row, ref int col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_FIND_ITEM searches an I4COL for a given value.
        //
        //  Discussion:
        //
        //    The two dimensional information in TABLE is stored as a one dimensional
        //    array, by columns.
        //
        //    The values ROW and COL will be one-based indices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns
        //    in the table.
        //
        //    Input, int TABLE[M*N], the table to search.
        //
        //    Input, int ITEM, the value to search for.
        //
        //    Output, int *ROW, *COL, the row and column indices
        //    of the first occurrence of the value ITEM.  The search
        //    is conducted by rows.  If the item is not found, then
        //    ROW = COL = -1.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (table[i + j * m] == item)
                {
                    row = i + 1;
                    col = j + 1;
                    return;
                }
            }
        }

        row = -1;
        col = -1;

    }

    public static void i4col_find_pair_wrap(int m, int n, int[] a, int item1, int item2,
            ref int row, ref int col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_FIND_PAIR_WRAP wrap searches an I4COL for a pair of items.
        //
        //  Discussion:
        //
        //    The items (ITEM1, ITEM2) must occur consecutively.
        //    However, wrapping is allowed, that is, if ITEM1 occurs
        //    in the last row, and ITEM2 "follows" it in the first row
        //    of the same column, a match is declared.
        //
        //    If the pair of items is not found, then ROW = COL = -1.
        //
        //    The values ROW and COL will be one-based indices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns
        //    in the table.
        //
        //    Input, int A[M*N], the table to search.
        //
        //    Input, int ITEM1, ITEM2, the values to search for.
        //
        //    Output, int *ROW, *COL, the row and column indices
        //    of the first occurrence of the value ITEM1 followed immediately
        //    by ITEM2.
        //
    {
        int j;

        for (j = 1; j <= n; j++)
        {
            int i;
            for (i = 1; i <= m; i++)
            {
                if (a[i - 1 + (j - 1) * m] == item1)
                {
                    int i2 = i + 1;

                    if (m < i2)
                    {
                        i2 = 1;
                    }

                    if (a[i2 - 1 + (j - 1) * m] == item2)
                    {
                        row = i;
                        col = j;
                        return;
                    }
                }
            }
        }

        row = -1;
        col = -1;
    }

    public static void i4col_sort2_a ( int m, int n, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_SORT2_A ascending sorts the elements of each column of an I4COL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of A.
        //
        //    Input, int N, the number of columns of A, and the length
        //    of a vector of data.
        //
        //    Input/output, int A[M*N].
        //    On input, the array of N columns of M vectors.
        //    On output, the elements of each column of A have been sorted in ascending
        //    order.
        //
    {
        int col;

        SortHeapExternalData data = new();

        switch (m)
        {
            case <= 1:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }
        //
        //  Initialize.
        //
        for ( col = 0; col < n; col++ )
        {
            int i = 0;
            int indx = 0;
            int isgn = 0;
            int j = 0;
            //
            //  Call the external heap sorter.
            //
            for ( ; ; )
            {
                Sort.sort_heap_external ( ref data, m, ref indx, ref i, ref j, isgn );
                //
                //  Interchange the I and J objects.
                //
                if ( 0 < indx )
                {
                    (a[i-1+col*m], a[j-1+col*m]) = (a[j-1+col*m], a[i-1+col*m]);
                }
                //
                //  Compare the I and J objects.
                //
                else if ( indx < 0 )
                {
                    if ( a[j-1+col*m] < a[i-1+col*m] )
                    {
                        isgn = +1;
                    }
                    else
                    {
                        isgn = -1;
                    }
                }
                else
                {
                    break;
                }
            }
        }
    }
    public static void i4col_sort_a(int m, int n, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_SORT_A ascending sorts an I4COL.
        //
        //  Discussion:
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 June 2005
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
        //    On input, the array of N columns of M vectors;
        //    On output, the columns of A have been sorted in ascending
        //    lexicographic order.
        //
    {
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
        SortHeapExternalData data = new();
        for (;;)
        {
            Sort.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);
            //
            //  Interchange the I and J objects.
            //
            if (0 < indx)
            {
                i4col_swap(m, n, ref a, i, j);
            }
            //
            //  Compare the I and J objects.
            //
            else if (indx < 0)
            {
                isgn = i4col_compare(m, n, a, i, j);
            }
            else
            {
                break;
            }
        }
    }

    public static int i4col_sorted_unique_count(int m, int n, int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
        //
        //  Discussion:
        //
        //    The columns of the array may be ascending or descending sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int A[M*N], a sorted array, containing
        //    N columns of data.
        //
        //    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
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

        unique_num = 1;
        int j1 = 0;

        for (j2 = 1; j2 < n; j2++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (a[i + j1 * m] == a[i + j2 * m])
                {
                    continue;
                }

                unique_num += 1;
                j1 = j2;
                break;
            }
        }

        return unique_num;
    }

    public static void i4col_swap(int m, int n, ref int[] a, int icol1, int icol2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4COL_SWAP swaps two columns of an I4COL.
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
        //    03 April 2005
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
        //    Input, int ICOL1, ICOL2, the two columns to swap.
        //    These indices should be between 1 and N.
        //
    {
        const int OFFSET = 1;

        int i;
        //
        //  Check.
        //
        if (icol1 - OFFSET < 0 || n - 1 < icol1 - OFFSET)
        {
            Console.WriteLine("");
            Console.WriteLine("I4COL_SWAP - Fatal error!");
            Console.WriteLine("  ICOL1 is out of range.");
            return;
        }

        if (icol2 - OFFSET < 0 || n - 1 < icol2 - OFFSET)
        {
            Console.WriteLine("");
            Console.WriteLine("I4COL_SWAP - Fatal error!");
            Console.WriteLine("  ICOL2 is out of range.");
            return;
        }

        if (icol1 == icol2)
        {
            return;
        }

        for (i = 0; i < m; i++)
        {
            (a[i + (icol1 - OFFSET) * m], a[i + (icol2 - OFFSET) * m]) = (a[i + (icol2 - OFFSET) * m], a[i + (icol1 - OFFSET) * m]);
        }
    }
}