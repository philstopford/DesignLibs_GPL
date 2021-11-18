using System;
using Burkardt.SortNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_undex(int x_num, double[] x_val, int x_unique_num, double tol,
            ref int[] undx, ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_UNDEX returns unique sorted indexes for an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points, to the unique elements of X, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of X, the index of
        //    the unique sorted element of X.
        //
        //    This is all done with index vectors, so that the elements of
        //    X are never moved.
        //
        //    The first step of the algorithm requires the indexed sorting
        //    of X, which creates arrays INDX and XDNI.  (If all the entries
        //    of X are unique, then these arrays are the same as UNDX and XDNU.)
        //
        //    We then use INDX to examine the entries of X in sorted order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector X could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of X in sorted order, using the formula
        //
        //      XU(*) = X(UNDX(*)).
        //
        //    We could then, if we wished, reconstruct the entire vector X, or
        //    any element of it, by index, as follows:
        //
        //      X(I) = XU(XDNU(I)).
        //
        //    We could then replace X by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of X, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector X, the sort and inverse sort
        //    index vectors, and the unique sort and inverse unique sort vectors
        //    and the compressed unique sorted vector.
        //
        //      I     X  Indx  Xdni       XU  Undx  Xdnu
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
        //    INDX(2) = 3 means that sorted item(2) is X(3).
        //    XDNI(2) = 5 means that X(2) is sorted item(5).
        //
        //    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
        //    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
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
        //    03 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of data values.
        //
        //    Input, double X_VAL[X_NUM], the data values.
        //
        //    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[X_NUM], the XDNU vector.
        //
    {
        //
        //  Implicitly sort the array.
        //
        int[] indx = r8vec_sort_heap_index_a_new(x_num, x_val);
        //
        //  Walk through the implicitly sorted array X.
        //
        int i = 0;

        int j = 0;
        undx[j] = indx[i];

        xdnu[indx[i]] = j;

        for (i = 1; i < x_num; i++)
        {
            if (tol < Math.Abs(x_val[indx[i]] - x_val[undx[j]]))
            {
                j += 1;
                undx[j] = indx[i];
            }

            xdnu[indx[i]] = j;
        }
    }

    public static void r8vec_sort_heap_a ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 1999
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
        //    Input, int N, the number of entries in the array.
        //
        //    Input/output, double A[N].
        //    On input, the array to be sorted;
        //    On output, the array has been sorted.
        //
    {
        int n1;

        switch (n)
        {
            case <= 1:
                return;
        }
        //
        //  1: Put A into descending heap form.
        //
        r8vec_heap_d ( n, ref a );
        //
        //  2: Sort A.
        //
        //  The largest object in the heap is in A[0].
        //  Move it to position A[N-1].
        //
        double temp = a[0];
        a[0] = a[n-1];
        a[n-1] = temp;
        //
        //  Consider the diminished heap of size N1.
        //
        for ( n1 = n-1; 2 <= n1; n1-- )
        {
            //
            //  Restore the heap structure of the initial N1 entries of A.
            //
            r8vec_heap_d ( n1, ref a );
            //
            //  Take the largest object from A[0] and move it to A[N1-1].
            //
            temp = a[0];
            a[0] = a[n1-1];
            a[n1-1] = temp;
        }
    }

    public static void r8vec_sort_insert_a(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Kreher, Douglas Simpson,
        //    Algorithm 1.1,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998, page 11.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in the vector.
        //    N must be positive.
        //
        //    Input/output, double A[N].
        //
        //    On input, A contains data to be sorted.
        //    On output, the entries of A have been sorted in ascending order.
        //
    {
        int i;

        for (i = 1; i < n; i++)
        {
            double x = a[i];

            int j = i;

            while (1 <= j && x < a[j - 1])
            {
                a[j] = a[j - 1];
                j -= 1;
            }

            a[j] = x;
        }
    }
        
    public static int[] r8vec_sort_insert_index_a ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998, page 11.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in the vector.
        //    N must be positive.
        //
        //    Input, double A[N], the array to be sorted.
        //
        //    Output, int R8VEC_SORT_INSERT_INDEX_A[N], the sorted indices.  The array
        //    is sorted when listed from A(INDX(1)) through A(INDX(N)).
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                return null;
        }

        int[] indx = i4vec_indicator0_new ( n );

        for ( i = 1; i < n; i++ )
        {
            double x = a[i];
            int j = i - 1;

            while ( 0 <= j )
            {
                if ( a[indx[j]] <= x )
                {
                    break;
                }

                indx[j+1] = indx[j];
                j -= 1;
            }
            indx[j+1] = i;
        }

        return indx;
    }
        
    public static void r8vec_sort_heap_d(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2005
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
        //    Input, int N, the number of entries in the array.
        //
        //    Input/output, double A[N].
        //    On input, the array to be sorted;
        //    On output, the array has been sorted.
        //
    {
        int n1;

        switch (n)
        {
            case <= 1:
                return;
        }

        //
        //  1: Put A into ascending heap form.
        //
        r8vec_heap_a(n, ref a);
        //
        //  2: Sort A.
        //
        //  The smallest object in the heap is in A[0].
        //  Move it to position A[N-1].
        //
        double temp = a[0];
        a[0] = a[n - 1];
        a[n - 1] = temp;
        //
        //  Consider the diminished heap of size N1.
        //
        for (n1 = n - 1; 2 <= n1; n1--)
        {
            //
            //  Restore the heap structure of the initial N1 entries of A.
            //
            r8vec_heap_a(n1, ref a);
            //
            //  Take the largest object from A[0] and move it to A[N1-1].
            //
            temp = a[0];
            a[0] = a[n1 - 1];
            a[n1 - 1] = temp;
        }

    }

    public static int[] r8vec_sort_heap_index_a(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
        //
        //  Discussion:
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      A(INDX(I)), I = 1 to N is sorted,
        //
        //    after which A(I), I = 1 to N is sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], an array to be index-sorted.
        //
        //    Output, int R8VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
        //    I-th element of the sorted array is A(INDX(I)).
        //
    {
        int i;

        int[] indx = new int[n];

        for (i = 1; i <= n; i++)
        {
            indx[i - 1] = i;
        }

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            double aval;
            int indxt;
            switch (l)
            {
                case > 1:
                    l -= 1;
                    indxt = indx[l - 1];
                    aval = a[indxt - 1];
                    break;
                default:
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt - 1];
                    indx[ir - 1] = indx[0];
                    ir -= 1;

                    switch (ir)
                    {
                        case 1:
                        {
                            indx[0] = indxt;
                            for (i = 0; i < n; i++)
                            {
                                indx[i] -= 1;
                            }

                            return indx;
                        }
                    }

                    break;
                }
            }

            i = l;
            int j = l + l;

            while (j <= ir)
            {
                if (j < ir)
                {
                    if (a[indx[j - 1] - 1] < a[indx[j] - 1])
                    {
                        j += 1;
                    }
                }

                if (aval < a[indx[j - 1] - 1])
                {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i - 1] = indxt;
        }
    }

    public static int[] r8vec_sort_heap_index_a_new(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC.
        //
        //  Discussion:
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      A(INDX(I)), I = 1 to N is sorted,
        //
        //    after which A(I), I = 1 to N is sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], an array to be index-sorted.
        //
        //    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
        //    I-th element of the sorted array is A(INDX(I)).
        //
    {
        int[] indx = new int[n];

        for (int i = 1; i <= n; i++)
        {
            indx[i - 1] = i;
        }

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            double aval;
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval = a[indxt - 1];
            }
            else
            {
                indxt = indx[ir - 1];
                aval = a[indxt - 1];
                indx[ir - 1] = indx[0];
                ir -= 1;

                if (ir == 1)
                {
                    indx[0] = indxt;
                    for (int i = 0; i < n; i++)
                    {
                        indx[i] -= 1;
                    }

                    break;
                }
            }

            int i2 = l;
            int j = l + l;

            while (j <= ir)
            {
                if (j < ir)
                {
                    if (a[indx[j - 1] - 1] < a[indx[j] - 1])
                    {
                        j += 1;
                    }
                }

                if (aval < a[indx[j - 1] - 1])
                {
                    indx[i2 - 1] = indx[j - 1];
                    i2 = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i2 - 1] = indxt;
        }

        return indx;
    }

    public static void r8vec_sort_heap_index_d(int n, double[] a, ref int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_INDEX_D_NEW: indexed heap descending sort of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      a(indx(*))
        //
        //    or explicitly, by the call
        //
        //      r8vec_permute ( n, indx, 0, a )
        //
        //    after which a(*) is sorted.
        //
        //    Note that the index vector is 0-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], an array to be index-sorted.
        //
        //    Output, int INDX[N], contains the sort index.  The
        //    I-th element of the sorted array is A(INDX(I)).
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                return;
        }

        for (i = 0; i < n; i++)
        {
            indx[i] = i;
        }

        switch (n)
        {
            case 1:
                return;
        }

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            double aval;
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval = a[indxt];
            }
            else
            {
                indxt = indx[ir - 1];
                aval = a[indxt];
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
                if (j < ir)
                {
                    if (a[indx[j]] < a[indx[j - 1]])
                    {
                        j += 1;
                    }
                }

                if (a[indx[j - 1]] < aval)
                {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i - 1] = indxt;
        }
    }

    public static int[] r8vec_sort_heap_index_d_new(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_INDEX_D_NEW: indexed heap descending sort of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      a(indx(*))
        //
        //    or explicitly, by the call
        //
        //      r8vec_permute ( n, indx, 0, a )
        //
        //    after which a(*) is sorted.
        //
        //    Note that the index vector is 0-based.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], an array to be index-sorted.
        //
        //    Output, int R8VEC_SORT_HEAP_INDEX_D[N], contains the sort index.  The
        //    I-th element of the sorted array is A(INDX(I)).
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
                return indx;
        }

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            double aval;
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval = a[indxt];
            }
            else
            {
                indxt = indx[ir - 1];
                aval = a[indxt];
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
                if (j < ir)
                {
                    if (a[indx[j]] < a[indx[j - 1]])
                    {
                        j += 1;
                    }
                }

                if (a[indx[j - 1]] < aval)
                {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i - 1] = indxt;
        }

        return indx;
    }

    public static int[] r8vec_sort_heap_mask_a(int n, double[] a, int mask_num, int[] mask)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An array A is given.  An array MASK of indices into A is given.
        //    The routine produces a vector INDX, which is a permutation of the
        //    entries of MASK, so that:
        //
        //      A(MASK(INDX(I)) <= A(MASK(INDX(J))
        //
        //    whenever
        //
        //      I <= J
        //
        //    In other words, only the elements of A that are indexed by MASK
        //    are to be considered, and the only thing that happens is that
        //    a rearrangment of the indices in MASK is returned that orders the
        //    masked elements.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[N], an array to be index-sorted.
        //
        //    Input, int MASK_NUM, the number of mask elements.
        //
        //    Input, int MASK[MASK_NUM], the mask array.  This is
        //    simply a list of indices of A.  The entries of MASK should
        //    be unique, and each one should be between 1 and N.
        //
        //    Output, int INDX[MASK_NUM], the sort index.  There are MASK_NUM
        //    elements of A selected by MASK.  If we want to list those elements
        //    in order, then the I-th element is A(MASK(INDX(I))).
        //
    {
        int i;
        int[] indx;

        switch (n)
        {
            case < 1:
                return null;
        }

        switch (mask_num)
        {
            case < 1:
                return null;
            case 1:
                indx = new int[1];
                indx[0] = 0;
                return indx;
        }

        indx = i4vec_indicator1_new(mask_num);

        int l = mask_num / 2 + 1;
        int ir = mask_num;

        for (;;)
        {
            double aval;
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval = a[mask[indxt - 1] - 1];
            }
            else
            {
                indxt = indx[ir - 1];
                aval = a[mask[indxt - 1] - 1];
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
                if (j < ir)
                {
                    if (a[mask[indx[j - 1] - 1] - 1] < a[mask[indx[j] - 1] - 1])
                    {
                        j += 1;
                    }
                }

                if (aval < a[mask[indx[j - 1] - 1] - 1])
                {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i - 1] = indxt;
        }

        for (i = 0; i < mask_num; i++)
        {
            indx[i] -= 1;
        }

        return indx;
    }


    public static void r8vec_sort_quick_a(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 7
        //
        //      A = ( 6, 7, 3, 2, 9, 1, 8 )
        //
        //    Output:
        //
        //      A = ( 1, 2, 3, 6, 7, 8, 9 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of A.
        //
        //    Input/output, double A[N].  On input, the array to be sorted.
        //    On output, A has been reordered into ascending order.
        //
    {
        const int LEVEL_MAX = 30;

        int l_segment = 0;
        int[] rsave = new int[LEVEL_MAX];
        int r_segment = 0;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_SORT_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case 1:
                return;
        }

        int level = 1;
        rsave[0] = n + 1;
        int base_ = 1;
        int n_segment = n;

        while (0 < n_segment)
        {
            //
            //  Partition the segment.
            //
            r8vec_part_quick_a(n_segment, ref a, base_ - 1, ref l_segment, ref r_segment);
            switch (l_segment)
            {
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                case > 1 when LEVEL_MAX < level:
                    Console.WriteLine("");
                    Console.WriteLine("R8VEC_SORT_QUICK_A - Fatal error!");
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
                            if (1 < level)
                            {
                                base_ = rsave[level - 1];
                                n_segment = rsave[level - 2] - rsave[level - 1];
                                level -= 1;
                                if (0 < n_segment)
                                {
                                    break;
                                }
                            }
                            else
                            {
                                n_segment = 0;
                                break;
                            }
                        }
                    }

                    break;
                }
            }
        }

    }

    public static void r8vec2_sort_a(int n, ref double[] a1, ref double[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SORT_A ascending sorts an R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    Each item to be sorted is a pair of reals (X,Y), with the X
        //    and Y values stored in separate vectors A1 and A2.
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
        //    Input, int N, the number of items of data.
        //
        //    Input/output, double A1[N], A2[N], the data to be sorted.
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
                double temp = a1[i - 1];
                a1[i - 1] = a1[j - 1];
                a1[j - 1] = temp;

                temp = a2[i - 1];
                a2[i - 1] = a2[j - 1];
                a2[j - 1] = temp;
            }
            //
            //  Compare the I and J objects.
            //
            else if (indx < 0)
            {
                isgn = r8vec2_compare(n, a1, a2, i, j);
            }
            else
            {
                break;
            }
        }
    }

    public static void r8vec2_sort_d(int n, double[] a1, double[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SORT_D descending sorts an R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    Each item to be sorted is a pair of reals (X,Y), with the X
        //    and Y values stored in separate vectors A1 and A2.
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
        //    Input, int N, the number of items of data.
        //
        //    Input/output, double A1[N], A2[N], the data to be sorted.
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
                double temp = a1[i - 1];
                a1[i - 1] = a1[j - 1];
                a1[j - 1] = temp;

                temp = a2[i - 1];
                a2[i - 1] = a2[j - 1];
                a2[j - 1] = temp;
            }
            //
            //  Compare the I and J objects.
            //
            else if (indx < 0)
            {
                isgn = -r8vec2_compare(n, a1, a2, i, j);
            }
            else
            {
                break;
            }
        }
    }

    public static int[] r8vec2_sort_heap_index_a(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:
        //
        //    * X(I) < X(J), or
        //
        //    * X(I) = X(J), and Y(I) < Y(J).
        //
        //    Once the index array is computed, the sorting can be carried out
        //    implicitly:
        //
        //      ( x(indx(*)), y(indx(*) )
        //
        //    or explicitly, by the calls
        //
        //      r8vec_permute ( n, indx, 0, x )
        //      r8vec_permute ( n, indx, 0, y )
        //
        //    after which ( x(*), y(*) ), is sorted.
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
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double X[N], Y[N], pairs of X, Y coordinates of points.
        //
        //    Output, int INDX[N], the sort index.  The
        //    I-th element of the sorted array has coordinates
        //    ( X(INDX(I)), Y(INDX(I) ).
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

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            int indxt;
            double yval;
            double xval;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                xval = x[indxt];
                yval = y[indxt];
            }
            else
            {
                indxt = indx[ir - 1];
                xval = x[indxt];
                yval = y[indxt];
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
                if (j < ir)
                {
                    if (x[indx[j - 1]] < x[indx[j]] ||
                        Math.Abs(x[indx[j - 1]] - x[indx[j]]) <= double.Epsilon && y[indx[j - 1]] < y[indx[j]])
                    {
                        j += 1;
                    }
                }

                if (xval < x[indx[j - 1]] ||
                    Math.Abs(xval - x[indx[j - 1]]) <= double.Epsilon && yval < y[indx[j - 1]])
                {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                }
                else
                {
                    j = ir + 1;
                }
            }

            indx[i - 1] = indxt;
        }

        return indx;
    }
        
    public static void r8vec2_sorted_unique(int n, double[] a1, double[] a2, ref int unique_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SORTED_UNIQUE keeps the unique elements in an R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    Item I is stored as the pair A1(I), A2(I).
        //
        //    The items must have been sorted, or at least it must be the
        //    case that equal items are stored in adjacent vector locations.
        //
        //    If the items were not sorted, then this routine will only
        //    replace a string of equal values by a single representative.
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
        //    Input, int N, the number of items.
        //
        //    Input/output, double A1[N], A2[N].
        //    On input, the array of N items.
        //    On output, an array of UNIQUE_NUM unique items.
        //
        //    Output, int &UNIQUE_NUM, the number of unique items.
        //
    {
        int itest;

        unique_num = 0;

        switch (n)
        {
            case <= 0:
                return;
        }

        unique_num = 1;

        for (itest = 1; itest < n; itest++)
        {
            if (!(Math.Abs(a1[itest] - a1[unique_num - 1]) > double.Epsilon) &&
                !(Math.Abs(a2[itest] - a2[unique_num - 1]) > double.Epsilon))
            {
                continue;
            }

            a1[unique_num] = a1[itest];
            a2[unique_num] = a2[itest];
            unique_num += 1;
        }
    }

    public static void r8vec2_sorted_unique_index(int n, ref double[] a1, ref double[] a2,
            ref int unique_num, ref int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.
        //
        //  Discussion:
        //
        //    An R8VEC2 is a dataset consisting of N pairs of real values, stored
        //    as two separate vectors A1 and A2.
        //
        //    Item I is stored as the pair A1(I), A2(I).
        //
        //    The items must have been sorted, or at least it should be the
        //    case that equal items are stored in adjacent vector locations.
        //
        //    If the items are not sorted, then this routine will only
        //    replace a string of equal values by a single representative.
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
        //    Input, int N, the number of items.
        //
        //    Input/output, double A1[N], A2[N].
        //    On input, the array of N items.
        //    On output, an array of unique items.
        //
        //    Output, int &UNIQUE_NUM, the number of unique items.
        //
        //    Output, int INDX[N], contains in entries 1 through UNIQUE_NUM an index
        //    array of the unique items.  To build new arrays with no repeated elements:
        //      B1(*) = A1(INDX(*))
        //
    {
        int itest;

        switch (n)
        {
            case <= 0:
                unique_num = 0;
                return;
        }

        i4vec_zeros(n, ref indx);

        unique_num = 1;
        indx[0] = 1;

        for (itest = 2; itest <= n; itest++)
        {
            if (!(Math.Abs(a1[itest - 2] - a1[itest - 1]) > double.Epsilon) &&
                !(Math.Abs(a2[itest - 2] - a2[itest - 1]) > double.Epsilon))
            {
                continue;
            }

            unique_num += 1;
            indx[unique_num - 1] = itest;
        }
    }

    public static int r8vec_sorted_unique_count(int n, double[] a, double tol)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Because the array is sorted, this algorithm is O(N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, double A[N], the sorted array to examine.
        //
        //    Input, double TOL, a tolerance for checking equality.
        //
        //    Output, int R8VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
        //
    {
        int i;

        int unique_num = 0;

        switch (n)
        {
            case < 1:
                return unique_num;
        }

        unique_num = 1;

        for (i = 1; i < n; i++)
        {
            if (tol < Math.Abs(a[i - 1] - a[i]))
            {
                unique_num += 1;
            }
        }

        return unique_num;
    }

    public static void r8vec_sorted_unique_hist(int n, double[] a, double tol, int maxuniq,
            ref int unique_num, double[] auniq, int[] acount)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_UNIQUE_HIST histograms unique elements of a sorted R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, double A[N], the array to examine, which must have been
        //    sorted.
        //
        //    Input, double TOL, a tolerance for checking equality.
        //
        //    Input, int MAXUNIQ, the maximum number of unique elements
        //    that can be handled.  If there are more than MAXUNIQ unique
        //    elements in A, the excess will be ignored.
        //
        //    Output, int &UNIQUE_NUM, the number of unique elements of A.
        //
        //    Output, double AUNIQ[UNIQUE_NUM], the unique elements of A.
        //
        //    Output, int ACOUNT[UNIQUE_NUM], the number of times each element
        //    of AUNIQ occurs in A.
        //
    {
        int i;
        //
        //  Start taking statistics.
        //
        int index = -1;

        for (i = 0; i < n; i++)
        {
            switch (i)
            {
                case 0:
                    index = 0;
                    auniq[index] = a[0];
                    acount[index] = 1;
                    break;
                default:
                {
                    if (Math.Abs(a[i] - auniq[index]) <= tol)
                    {
                        acount[index] += 1;
                    }
                    else if (index + 1 < maxuniq)
                    {
                        index += 1;
                        auniq[index] = a[i];
                        acount[index] = 1;
                    }

                    break;
                }
            }
        }

        unique_num = index + 1;
    }

    public static void r8vec_sort_shell_a(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input/output, double A[N].
        //    On input, an array to be sorted.
        //    On output, the sorted array.
        //
    {
        int ipow;

        switch (n)
        {
            case <= 1:
                return;
        }

        //
        //  Determine the smallest MAXPOW so that
        //    N <= ( 3^MAXPOW - 1 ) / 2
        //
        int maxpow = 1;
        int test = 3;

        while (test < 2 * n + 1)
        {
            maxpow += 1;
            test *= 3;
        }

        switch (maxpow)
        {
            case > 1:
                maxpow -= 1;
                test /= 3;
                break;
        }

        //
        //  Now sort groups of size ( 3^IPOW - 1 ) / 2.
        //
        for (ipow = maxpow; 1 <= ipow; ipow--)
        {
            int inc = (test - 1) / 2;
            test /= 3;
            //
            //  Sort the values with indices equal to K mod INC.
            //
            int k;
            for (k = 1; k <= inc; k++)
            {
                //
                //  Insertion sort of the items with index
                //  INC+K, 2*INC+K, 3*INC+K, ...
                //
                int i;
                for (i = inc + k; i <= n; i += inc)
                {
                    double asave = a[i - 1];
                    int ifree = i;
                    int j = i - inc;

                    for (;;)
                    {
                        if (j < 1)
                        {
                            break;
                        }

                        if (a[j - 1] <= asave)
                        {
                            break;
                        }

                        ifree = j;
                        a[j + inc - 1] = a[j - 1];
                        j -= inc;
                    }

                    a[ifree - 1] = asave;
                }
            }
        }
    }

    public static double[] r8vec_sorted_merge_a(int na, double[] a, int nb, double[] b, ref int nc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The elements of A and B should be sorted in ascending order.
        //
        //    The elements in the output array C will also be in ascending order,
        //    and unique.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NA, the dimension of A.
        //
        //    Input, double A[NA], the first sorted array.
        //
        //    Input, int NB, the dimension of B.
        //
        //    Input, double B[NB], the second sorted array.
        //
        //    Output, int &NC, the number of entries in the merged vector.
        //
        //    Output, double R8VEC_SORTED_MERGE_A[NC], the merged unique sorted array.
        //
    {
        int na2 = na;
        int nb2 = nb;

        int ja = 0;
        int jb = 0;
        nc = 0;
        int nd = 0;
        double[] d = new double[na + nb];

        int order = r8vec_order_type(na2, a);

        switch (order)
        {
            case < 0:
            case > 2:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_SORTED_MERGE_A - Fatal error!");
                Console.WriteLine("  The input array A is not ascending sorted.");
                return null;
        }

        order = r8vec_order_type(nb2, b);

        switch (order)
        {
            case < 0:
            case > 2:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_SORTED_MERGE_A - Fatal error!");
                Console.WriteLine("  The input array B is not ascending sorted.");
                return null;
        }

        for (;;)
        {
            //
            //  If we've used up all the entries of A, stick the rest of B on the end.
            //
            int j;
            if (na2 <= ja)
            {
                for (j = 1; j <= nb2 - jb; j++)
                {
                    jb += 1;
                    switch (nd)
                    {
                        case 0:
                            nd += 1;
                            d[nd - 1] = b[jb - 1];
                            break;
                        default:
                        {
                            if (d[nd - 1] < b[jb - 1])
                            {
                                nd += 1;
                                d[nd - 1] = b[jb - 1];
                            }

                            break;
                        }
                    }
                }

                break;
            }
            //
            //  If we've used up all the entries of B, stick the rest of A on the end.
            //

            if (nb2 <= jb)
            {
                for (j = 1; j <= na2 - ja; j++)
                {
                    ja += 1;
                    switch (nd)
                    {
                        case 0:
                            nd += 1;
                            d[nd - 1] = a[ja - 1];
                            break;
                        default:
                        {
                            if (d[nd - 1] < a[ja - 1])
                            {
                                nd += 1;
                                d[nd - 1] = a[ja - 1];
                            }

                            break;
                        }
                    }
                }

                break;
            }
            //
            //  Otherwise, if the next entry of A is smaller, that's our candidate.
            //
            if (a[ja] <= b[jb])
            {
                ja += 1;
                switch (nd)
                {
                    case 0:
                        nd += 1;
                        d[nd - 1] = a[ja - 1];
                        break;
                    default:
                    {
                        if (d[nd - 1] < a[ja - 1])
                        {
                            nd += 1;
                            d[nd - 1] = a[ja - 1];
                        }

                        break;
                    }
                }
            }
            //
            //  ...or if the next entry of B is the smaller, consider that.
            //
            else
            {
                jb += 1;
                switch (nd)
                {
                    case 0:
                        nd += 1;
                        d[nd - 1] = b[jb - 1];
                        break;
                    default:
                    {
                        if (d[nd - 1] < b[jb - 1])
                        {
                            nd += 1;
                            d[nd - 1] = b[jb - 1];
                        }

                        break;
                    }
                }
            }
        }

        nc = nd;

        double[] c = r8vec_copy_new(nd, d);

        return c;
    }

    public static int r8vec_sorted_nearest(int n, double[] a, double value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of A.
        //
        //    Input, double A[N], a sorted vector.
        //
        //    Input, double VALUE, the value whose nearest vector entry is sought.
        //
        //    Output, int R8VEC_SORTED_NEAREST, the index of the nearest
        //    entry in the vector.
        //
    {
        int hi;
        int lo;
        int mid;

        switch (n)
        {
            case < 1:
                return -1;
            case 1:
                return 1;
        }

        if (a[0] < a[n - 1])
        {
            if (value < a[0])
            {
                return 1;
            }

            if (a[n - 1] < value)
            {
                return n;
            }

            //
            //  Seek an interval containing the value.
            //
            lo = 1;
            hi = n;

            while (lo < hi - 1)
            {
                mid = (lo + hi) / 2;

                if (Math.Abs(value - a[mid - 1]) <= double.Epsilon)
                {
                    return mid;
                }

                if (value < a[mid - 1])
                {
                    hi = mid;
                }
                else
                {
                    lo = mid;
                }
            }

            //
            //  Take the nearest.
            //
            return Math.Abs(value - a[lo - 1]) < Math.Abs(value - a[hi - 1]) ? lo : hi;
        }
        //
        //  A descending sorted vector A.
        //

        if (value < a[n - 1])
        {
            return n;
        }

        if (a[0] < value)
        {
            return 1;
        }

        //
        //  Seek an interval containing the value.
        //
        lo = n;
        hi = 1;

        while (lo < hi - 1)
        {
            mid = (lo + hi) / 2;

            if (Math.Abs(value - a[mid - 1]) <= double.Epsilon)
            {
                return mid;
            }

            if (value < a[mid - 1])
            {
                hi = mid;
            }
            else
            {
                lo = mid;
            }
        }

        //
        //  Take the nearest.
        //
        return Math.Abs(value - a[lo - 1]) < Math.Abs(value - a[hi - 1]) ? lo : hi;
    }

    public static void r8vec_sorted_range(int n, double[] r, double r_lo, double r_hi,
            ref int i_lo, ref int i_hi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_RANGE searches a sorted vector for elements in a range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in the vector.
        //
        //    Input, double R[N], the sorted vector.
        //
        //    Input, double R_LO, R_HI, the limits of the range.
        //
        //    Output, int &I_LO, &I_HI, the range of indices
        //    so that I_LO <= I <= I_HI => R_LO <= R(I) <= R_HI.  If no
        //    values in R lie in the range, then I_HI < I_LO will be returned.
        //
    {
        int i1;
        int i2;
        int j1;
        int j2;
        //
        //  Cases we can handle immediately.
        //
        if (r[n - 1] < r_lo)
        {
            i_lo = -1;
            i_hi = -2;
            return;
        }

        if (r_hi < r[0])
        {
            i_lo = -1;
            i_hi = -2;
            return;
        }

        switch (n)
        {
            //
            //  Are there are least two intervals?
            //
            case 1:
            {
                if (r_lo <= r[0] && r[0] <= r_hi)
                {
                    i_lo = 1;
                    i_hi = 1;
                }
                else
                {
                    i_lo = -1;
                    i_hi = -2;
                }

                return;
            }
        }

        //
        //  Bracket R_LO.
        //
        if (r_lo <= r[0])
        {
            i_lo = 0;
        }
        else
        {
            //
            //  R_LO is in one of the intervals spanned by R(J1) to R(J2).
            //  Examine the intermediate interval [R(I1), R(I1+1)].
            //  Does R_LO lie here, or below or above?
            //
            j1 = 0;
            j2 = n - 1;
            i1 = (j2 - 1) / 2;
            i2 = i1 + 1;

            for (;;)
            {
                if (r_lo < r[i1])
                {
                    j2 = i1;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else if (r[i2] < r_lo)
                {
                    j1 = i2;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else
                {
                    i_lo = i1;
                    break;
                }
            }
        }

        //
        //  Bracket R_HI
        //
        if (r[n - 1] <= r_hi)
        {
            i_hi = n - 1;
        }
        else
        {
            j1 = i_lo;
            j2 = n - 1;
            i1 = (j1 + j2 - 1) / 2;
            i2 = i1 + 1;

            for (;;)
            {
                if (r_hi < r[i1])
                {
                    j2 = i1;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else if (r[i2] < r_hi)
                {
                    j1 = i2;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else
                {
                    i_hi = i2;
                    break;
                }
            }
        }

        //
        //  We expect to have computed the largest I_LO and smallest I_HI such that
        //    R(I_LO) <= R_LO <= R_HI <= R(I_HI)
        //  but what we want is actually
        //    R_LO <= R(I_LO) <= R(I_HI) <= R_HI
        //  which we can usually get simply by incrementing I_LO and decrementing I_HI.
        //
        if (r[i_lo] < r_lo)
        {
            i_lo += 1;
            if (n - 1 < i_lo)
            {
                i_hi = i_lo - 1;
            }
        }

        if (!(r_hi < r[i_hi]))
        {
            return;
        }

        i_hi -= 1;
        i_lo = i_hi switch
        {
            < 0 => i_hi + 1,
            _ => i_lo
        };
    }

    public static void r8vec_sorted_split(int n, double[] a, double split, ref int i_lt,
            ref int i_gt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Given a splitting value SPLIT, the routine seeks indices
        //    I_LT and I_GT so that
        //
        //      A(I_LT) < SPLIT < A(I_GT),
        //
        //    and if there are intermediate index values between I_LT and
        //    I_GT, then those entries of A are exactly equal to SPLIT.
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
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int N, the number of entries in A.
        //
        //    Input, double A[N], a sorted array.
        //
        //    Input, double SPLIT, a value to which the entries in A are
        //    to be compared.
        //
        //    Output, int &I_LT:
        //    0 if no entries are less than SPLIT;
        //    N if all entries are less than SPLIT;
        //    otherwise, the index of the last entry in A less than SPLIT.
        //
        //    Output, int &I_GT:
        //    1 if all entries are greater than SPLIT;
        //    N+1 if no entries are greater than SPLIT;
        //    otherwise the index of the first entry in A greater than SPLIT.
        //
    {
        int i;

        switch (n)
        {
            case < 1:
                i_lt = -1;
                i_gt = -1;
                return;
        }

        if (split < a[0])
        {
            i_lt = 0;
            i_gt = 1;
            return;
        }

        if (a[n - 1] < split)
        {
            i_lt = n;
            i_gt = n + 1;
            return;
        }

        int lo = 1;
        int hi = n;

        for (;;)
        {
            if (lo + 1 == hi)
            {
                i_lt = lo;
                break;
            }

            int mid = (lo + hi) / 2;

            if (split <= a[mid - 1])
            {
                hi = mid;
            }
            else
            {
                lo = mid;
            }
        }

        for (i = i_lt + 1; i <= n; i++)
        {
            if (!(split < a[i - 1]))
            {
                continue;
            }

            i_gt = i;
            return;
        }

        i_gt = n + 1;
    }

    public static void r8vec_sorted_undex(int x_num, double[] x_val, int x_unique_num,
            double tol, ref int[] undx, ref int[] xdnu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORTED_UNDEX returns unique sorted indexes for a sorted R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The goal of this routine is to determine a vector UNDX,
        //    which points, to the unique elements of X, in sorted order,
        //    and a vector XDNU, which identifies, for each entry of X, the index of
        //    the unique sorted element of X.
        //
        //    This is all done with index vectors, so that the elements of
        //    X are never moved.
        //
        //    Assuming X is already sorted, we examine the entries of X in order,
        //    noting the unique entries, creating the entries of XDNU and
        //    UNDX as we go.
        //
        //    Once this process has been completed, the vector X could be
        //    replaced by a compressed vector XU, containing the unique entries
        //    of X in sorted order, using the formula
        //
        //      XU(I) = X(UNDX(I)).
        //
        //    We could then, if we wished, reconstruct the entire vector X, or
        //    any element of it, by index, as follows:
        //
        //      X(I) = XU(XDNU(I)).
        //
        //    We could then replace X by the combination of XU and XDNU.
        //
        //    Later, when we need the I-th entry of X, we can locate it as
        //    the XDNU(I)-th entry of XU.
        //
        //    Here is an example of a vector X, the sort and inverse sort
        //    index vectors, and the unique sort and inverse unique sort vectors
        //    and the compressed unique sorted vector.
        //
        //      I      X      XU  Undx  Xdnu
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
        //    INDX(2) = 3 means that sorted item(2) is X(3).
        //    XDNI(2) = 5 means that X(2) is sorted item(5).
        //
        //    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
        //    XDNU(8) = 2 means that X(8) is at unique sorted item(2).
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
        //    01 November 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of data values.
        //
        //    Input, double X_VAL[X_NUM], the data values.
        //
        //    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
        //    This value is only required for languages in which the size of
        //    UNDX must be known in advance.
        //
        //    Input, double TOL, a tolerance for equality.
        //
        //    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.
        //
        //    Output, int XDNU[X_NUM], the XDNU vector.
        //
    {
        //
        //  Walk through the sorted array X.
        //
        int i = 0;

        int j = 0;
        undx[j] = i;

        xdnu[i] = j;

        for (i = 1; i < x_num; i++)
        {
            if (tol < Math.Abs(x_val[i] - x_val[undx[j]]))
            {
                j += 1;
                undx[j] = i;
            }

            xdnu[i] = j;
        }
    }

    public static void r8vec_sort_bubble_a ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input/output, double A[N].
        //    On input, an unsorted array of doubles.
        //    On output, A has been sorted.
        //
    {
        int i;

        for ( i = 0; i < n-1; i++ )
        {
            int j;
            for ( j = i+1; j < n; j++ )
            {
                if ( a[j] < a[i] )
                {
                    (a[i], a[j]) = (a[j], a[i]);
                }
            }
        }
    }

    public static void r8vec_sort_bubble_d ( int n, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input/output, double A[N].
        //    On input, an unsorted array of doubles.
        //    On output, A has been sorted.
        //
    {
        int i;

        for ( i = 0; i < n-1; i++ )
        {
            int j;
            for ( j = i+1; j < n; j++ )
            {
                if ( a[i] < a[j] )
                {
                    (a[i], a[j]) = (a[j], a[i]);
                }
            }
        }
    }
        
}