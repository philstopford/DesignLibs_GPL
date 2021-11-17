using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_index_delete_all(int n, double[] x, int[] indx, double xval,
            ref int n2, ref double[] x2, ref int[] indx2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Note that the value of N is adjusted because of the deletions.
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
        //    Input, int N, the size of the current list.
        //
        //    Input, double X[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, the value to be sought.
        //
        //    Output, int &N2, the size of the current list.
        //
        //    Output, double X2[N2], the list.
        //
        //    Output, int INDX2[N2], the sort index of the list.
        //
    {
        int equal = 0;
        int get;
        int i;
        int less = 0;
        int more = 0;

        switch (n)
        {
            case < 1:
                n2 = 0;
                return;
        }

        i4vec_copy(n, indx, ref indx2);
        r8vec_copy(n, x, ref x2);
        n2 = n;

        r8vec_index_search(n2, x2, indx2, xval, ref less, ref equal, ref more);

        switch (equal)
        {
            case 0:
                return;
        }

        int equal1 = equal;

        for (;;)
        {
            if (equal1 <= 1)
            {
                break;
            }

            if (Math.Abs(x2[indx2[equal1 - 2] - 1] - xval) > double.Epsilon)
            {
                break;
            }

            equal1 -= 1;
        }

        int equal2 = equal;

        for (;;)
        {
            if (n2 <= equal2)
            {
                break;
            }

            if (Math.Abs(x2[indx2[equal2] - 1] - xval) > double.Epsilon)
            {
                break;
            }

            equal2 += 1;
        }

        //
        //  Discard certain X values.
        //
        int put = 0;

        for (get = 1; get <= n2; get++)
        {
            if (Math.Abs(x2[get - 1] - xval) > double.Epsilon)
            {
                put += 1;
                x2[put - 1] = x2[get - 1];
            }
        }

        //
        //  Adjust the INDX values.
        //
        for (equal = equal1; equal <= equal2; equal++)
        {
            for (i = 1; i <= n2; i++)
            {
                if (indx2[equal - 1] < indx2[i - 1])
                {
                    indx2[i - 1] -= 1;
                }
            }
        }

        //
        //  Discard certain INDX values.
        //
        for (i = 0; i <= n2 - equal2 - 1; i++)
        {
            indx2[equal1 + i - 1] = indx2[equal2 + i];
        }

        for (i = n2 + equal1 - equal2; i <= n2; i++)
        {
            indx2[i - 1] = 0;
        }

        //
        //  Adjust N.
        //
        n2 = put;

    }

    public static void r8vec_index_delete_dupes(int n, double[] x, int[] indx,
            ref int n2, ref double[] x2, ref int[] indx2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted list.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The output quantities N2, X2, and INDX2 are computed from the
        //    input quantities by sorting, and eliminating duplicates.
        //
        //    The output arrays should be dimensioned of size N, unless the user
        //    knows in advance what the value of N2 will be.
        //
        //    The output arrays may be identified with the input arrays.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the input list.
        //
        //    Input, double X[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Output, int &N2, the number of unique entries in X.
        //
        //    Output, double X2[N2], a copy of the list which has
        //    been sorted, and made unique.
        //
        //    Output, int INDX2[N2], the sort index of the new list.
        //
    {
        int i = 0;
        int n3 = 0;
        double[] x3 = new double[n];

        for (;;)
        {
            i += 1;

            if (n < i)
            {
                break;
            }

            switch (i)
            {
                case > 1 when x[indx[i - 1] - 1] == x3[n3 - 1]:
                    continue;
                default:
                    n3 += 1;
                    x3[n3 - 1] = x[indx[i - 1] - 1];
                    break;
            }
        }

        //
        //  Set the output data.
        //
        n2 = n3;
        r8vec_copy(n3, x3, ref x2);
        for (i = 0; i < n3; i++)
        {
            indx2[i] = i + 1;
        }
    }

    public static void r8vec_index_delete_one(int n, double[] x, int[] indx, double xval,
            ref int n2, ref double[] x2, ref int[] indx2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted list.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    If the value occurs in the list more than once, only one copy is deleted.
        //
        //    Note that the value of N is adjusted because of the deletions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 October 2000
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the current list.
        //
        //    Input, double X[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, the value to be sought.
        //
        //    Output, int &N2, the size of the current list.
        //
        //    Output, double X2[N2], the list.
        //
        //    Output, int INDX2[N2], the sort index of the list.
        //
    {
        int equal = 0;
        int less = 0;
        int more = 0;

        switch (n)
        {
            case < 1:
                n2 = 0;
                return;
        }

        n2 = n;
        i4vec_copy(n2, indx, ref indx2);
        r8vec_copy(n2, x, ref x2);

        r8vec_index_search(n2, x2, indx2, xval, ref less, ref equal, ref more);

        if (equal != 0)
        {
            int j = indx2[equal - 1];
            int i;
            for (i = j; i <= n2 - 1; i++)
            {
                x2[i - 1] = x[i];
            }

            for (i = equal; i <= n2 - 1; i++)
            {
                indx2[i - 1] = indx2[i];
            }

            for (i = 1; i <= n2 - 1; i++)
            {
                if (j < indx2[i - 1])
                {
                    indx2[i - 1] -= 1;
                }
            }

            n2 -= 1;
        }
    }

    public static void r8vec_index_insert(ref int n, double[] x, int[] indx, double xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_INSERT inserts a value in an indexed sorted list.
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
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N, the size of the current list.
        //
        //    Input, double X[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, the value to be sought.
        //
    {
        int equal = 0;
        int i;
        int less = 0;
        int more = 0;

        switch (n)
        {
            case <= 0:
                n = 1;
                x[0] = xval;
                indx[0] = 1;
                return;
        }

        r8vec_index_search(n, x, indx, xval, ref less, ref equal, ref more);

        x[n] = xval;
        for (i = n; more <= i; i--)
        {
            indx[i] = indx[i - 1];
        }

        indx[more - 1] = n + 1;
        n += 1;
    }

    public static void r8vec_index_insert_unique(ref int n, ref double[] x, ref int[] indx, double xval)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted list.
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
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &N, the size of the current list.
        //    If the input value XVAL does not already occur in X, then N is increased.
        //
        //    Input/output, double X[N], the list.
        //    If the input value XVAL does not already occur in X, then it is added
        //    to X.
        //
        //    Input/output, int INDX[N], the sort index of the list.
        //    If the input value XVAL does not already occur in X, then INDX is updated.
        //
        //    Input, double XVAL, the value which will be inserted into the X
        //    vector if it is not there already.
        //
    {
        int equal = 0;
        int less = 0;
        int more = 0;

        switch (n)
        {
            case <= 0:
                n = 1;
                x[0] = xval;
                indx[0] = 1;
                return;
        }

        //
        //  Does XVAL already occur in X?
        //
        r8vec_index_search(n, x, indx, xval, ref less, ref equal, ref more);

        switch (equal)
        {
            case 0:
            {
                x[n] = xval;
                int i;
                for (i = n; more <= i; i--)
                {
                    indx[i] = indx[i - 1];
                }

                indx[more - 1] = n + 1;
                n += 1;
                break;
            }
        }
    }

    public static void r8vec_index_order(int n, ref double[] x, int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The index vector itself is not modified.  Therefore, the pair
        //    (X,INDX) no longer represents an index sorted vector.
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
        //    Input, int N, the size of the current list.
        //
        //    Input/output, double X[N], the list.  On output, the list
        //    has been sorted.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
    {
        int i;

        double[] y = new double[n];

        for (i = 0; i < n; i++)
        {
            y[i] = x[indx[i] - 1];
        }

        for (i = 0; i < n; i++)
        {
            x[i] = y[i];
        }
    }

    public static void r8vec_index_search(int n, double[] x, int[] indx, double xval, ref int less,
            ref int equal, ref int more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.
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
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the current list.
        //
        //    Input, double X[N], the list.
        //
        //    Input, int INDX[N], the sort index of the list.
        //
        //    Input, double XVAL, the value to be sought.
        //
        //    Output, int &LESS, &EQUAL, &MORE, the indexes in INDX of the
        //    entries of X that are just less than, equal to, and just greater
        //    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
        //    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
        //    is the greatest entry of X, then MORE is N+1.
        //
    {
        switch (n)
        {
            case <= 0:
                less = 0;
                equal = 0;
                more = 0;
                return;
        }

        int lo = 1;
        int hi = n;
        double xlo = x[indx[lo - 1] - 1];
        double xhi = x[indx[hi - 1] - 1];

        if (xval < xlo)
        {
            less = 0;
            equal = 0;
            more = 1;
            return;
        }

        if (Math.Abs(xval - xlo) <= double.Epsilon)
        {
            less = 0;
            equal = 1;
            more = 2;
            return;
        }

        if (xhi < xval)
        {
            less = n;
            equal = 0;
            more = n + 1;
            return;
        }

        if (Math.Abs(xval - xhi) <= double.Epsilon)
        {
            less = n - 1;
            equal = n;
            more = n + 1;
            return;
        }

        for (;;)
        {
            if (lo + 1 == hi)
            {
                less = lo;
                equal = 0;
                more = hi;
                return;
            }

            int mid = (lo + hi) / 2;
            double xmid = x[indx[mid - 1] - 1];

            if (Math.Abs(xval - xmid) <= double.Epsilon)
            {
                equal = mid;
                less = mid - 1;
                more = mid + 1;
                return;
            }

            if (xval < xmid)
            {
                hi = mid;
            }
            else if (xmid < xval)
            {
                lo = mid;
            }
        }
    }

    public static void r8vec_index_sort_unique(int n, double[] x, ref int n2, ref double[] x2,
            ref int[] indx2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_SORT_UNIQUE creates a sort index for an R8VEC.
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
        //    16 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the current list.
        //
        //    Input, double X[N], the list.
        //
        //    Output, int &N2, the number of unique elements in X.
        //
        //    Output, double X2[N2], a list of the unique elements of X.
        //
        //    Output, int INDX2[N2], the sort index of the list.
        //
    {
        int i;

        n2 = 0;

        for (i = 0; i < n; i++)
        {
            r8vec_index_insert_unique(ref n2, ref x2, ref indx2, x[i]);
        }

        for (i = n2; i < n; i++)
        {
            x2[i] = -1;
        }

        for (i = n2; i < n; i++)
        {
            indx2[i] = -1;
        }
    }

    public static void r8vec_index_sorted_range(int n, double[] r, int[] indx, double r_lo,
            double r_hi, ref int i_lo, ref int i_hi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in the vector.
        //
        //    Input, double R[N], the index sorted vector.
        //
        //    Input, int INDX[N], the vector used to sort R.
        //    The vector R[INDX[*]] is sorted.
        //
        //    Input, double R_LO, R_HI, the limits of the range.
        //
        //    Output, int &I_LO, &I_HI, the range of indices
        //    so that I_LO <= I <= I_HI => R_LO <= R[INDX[I]] <= R_HI.  If no
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
        if (r[indx[n - 1]] < r_lo)
        {
            i_lo = n;
            i_hi = n - 1;
            return;
        }

        if (r_hi < r[indx[0]])
        {
            i_lo = 0;
            i_hi = -1;
            return;
        }

        switch (n)
        {
            //
            //  Are there are least two intervals?
            //
            case 1:
            {
                if (r_lo <= r[indx[0]] && r[indx[0]] <= r_hi)
                {
                    i_lo = 0;
                    i_hi = 0;
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
        if (r_lo <= r[indx[0]])
        {
            i_lo = 0;
        }
        else
        {
            //
            //  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
            //  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
            //  Does R_LO lie here, or below or above?
            //
            j1 = 0;
            j2 = n - 1;
            i1 = (j1 + j2 - 1) / 2;
            i2 = i1 + 1;

            for (;;)
            {
                if (r_lo < r[indx[i1]])
                {
                    j2 = i1;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else if (r[indx[i2]] < r_lo)
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
        //  Bracket R_HI.
        //
        if (r[indx[n - 1]] <= r_hi)
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
                if (r_hi < r[indx[i1]])
                {
                    j2 = i1;
                    i1 = (j1 + j2 - 1) / 2;
                    i2 = i1 + 1;
                }
                else if (r[indx[i2]] < r_hi)
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
        //    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
        //  but what we want is actually
        //    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
        //  which we can usually get simply by incrementing I_LO and decrementing I_HI.
        //
        if (r[indx[i_lo]] < r_lo)
        {
            i_lo += 1;
            if (n - 1 < i_lo)
            {
                i_hi = i_lo - 1;
            }
        }

        if (r_hi < r[indx[i_hi]])
        {
            i_hi -= 1;
            i_lo = i_hi switch
            {
                < 0 => i_hi + 1,
                _ => i_lo
            };
        }
    }

    public static void r8vec_indexed_heap_d(int n, double[] a, ref int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEXED_HEAP_D creates a descending heap from an indexed R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
        //    each referencing an entry of the data vector.
        //
        //    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
        //    we have:
        //      A[INDX[2*J+1]]   <= A[INDX[J]]
        //    and
        //      A[INDX[2*J+2]] <= A[INDX[J]]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        //
        //  Parameters:
        //
        //    Input, int N, the size of the index array.
        //
        //    Input, double A[*], the data vector.
        //
        //    Input/output, int INDX[N], the index array.
        //    Each entry of INDX must be a valid index for the array A.
        //    On output, the indices have been reordered into a descending heap.
        //
    {
        int i;
        //
        //  Only nodes N/2 - 1 down to 0 can be "parent" nodes.
        //
        for (i = n / 2 - 1; 0 <= i; i--)
        {
            //
            //  Copy the value out of the parent node.
            //  Position IFREE is now "open".
            //
            int key = indx[i];
            int ifree = i;

            for (;;)
            {
                //
                //  Positions 2*IFREE+1 and 2*IFREE+2 are the descendants of position
                //  IFREE.  (One or both may not exist because they exceed N-1.)
                //
                int m = 2 * ifree + 1;
                //
                //  Does the first position exist?
                //
                if (n - 1 < m)
                {
                    break;
                }

                //
                //  Does the second position exist?
                //
                if (m + 1 <= n - 1)
                {
                    //
                    //  If both positions exist, take the larger of the two values,
                    //  and update M if necessary.
                    //
                    if (a[indx[m]] < a[indx[m + 1]])
                    {
                        m += 1;
                    }
                }

                //
                //  If the large descendant is larger than KEY, move it up,
                //  and update IFREE, the location of the free position, and
                //  consider the descendants of THIS position.
                //
                if (a[indx[m]] <= a[key])
                {
                    break;
                }

                indx[ifree] = indx[m];
                ifree = m;
            }

            //
            //  Once there is no more shifting to do, KEY moves into the free spot IFREE.
            //
            indx[ifree] = key;
        }
    }

    public static int r8vec_indexed_heap_d_extract(ref int n, double[] a, ref int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
        //    each referencing an entry of the data vector.
        //
        //    The routine finds the maximum value in the heap, returns that value to the
        //    user, deletes that value from the heap, and restores the heap to its
        //    proper form.
        //
        //    Note that the argument N must be a variable, which will be decremented
        //    before return, and that INDX will hold one less value on output than it
        //    held on input.
        //
        //    This is one of three functions needed to model a priority queue.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Cormen, Charles Leiserson, Ronald Rivest,
        //    Introduction to Algorithms,
        //    MIT Press, 2001,
        //    ISBN: 0262032937,
        //    LC: QA76.C662.
        //
        //  Parameters:
        //
        //    Input/output, int &N, the number of items in the index vector.
        //
        //    Input, double A[*], the data vector.
        //
        //    Input/output, int INDX[N], the index vector.
        //
        //    Output, int R8VEC_INDEXED_HEAP_D_EXTRACT, the index in A of the item of
        //    maximum value, which has now been removed from the heap.
        //
    {
        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!");
                Console.WriteLine("  The heap is empty.");
                return 1;
        }

        //
        //  Get the index of the maximum value.
        //
        int indx_extract = indx[0];

        switch (n)
        {
            case 1:
                n = 0;
                return indx_extract;
        }

        //
        //  Shift the last index down.
        //
        indx[0] = indx[n - 1];
        //
        //  Restore the heap structure.
        //
        n -= 1;
        r8vec_indexed_heap_d(n, a, ref indx);

        return indx_extract;
    }

    public static void r8vec_indexed_heap_d_insert(ref int n, double[] a, ref int[] indx,
            int indx_insert)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
        //    each referencing an entry of the data vector.
        //
        //    Note that the argument N must be a variable, and will be incremented before
        //    return, and that INDX must be able to hold one more entry on output than
        //    it held on input.
        //
        //    This is one of three functions needed to model a priority queue.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Cormen, Charles Leiserson, Ronald Rivest,
        //    Introduction to Algorithms,
        //    MIT Press, 2001,
        //    ISBN: 0262032937,
        //    LC: QA76.C662.
        //
        //  Parameters:
        //
        //    Input/output, int &N, the number of items in the index vector.
        //
        //    Input, double A[*], the data vector.
        //
        //    Input/output, int INDX[N], the index vector.
        //
        //    Input, int INDX_INSERT, the index in A of the value
        //    to be inserted into the heap.
        //
    {
        n += 1;
        int i = n - 1;

        while (0 < i)
        {
            int parent = (i - 1) / 2;

            if (a[indx_insert] <= a[indx[parent]])
            {
                break;
            }

            indx[i] = indx[parent];
            i = parent;
        }

        indx[i] = indx_insert;
    }

    public static int r8vec_indexed_heap_d_max(int n, double[] a, int[] indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
        //    each referencing an entry of the data vector.
        //
        //    This is one of three functions needed to model a priority queue.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Thomas Cormen, Charles Leiserson, Ronald Rivest,
        //    Introduction to Algorithms,
        //    MIT Press, 2001,
        //    ISBN: 0262032937,
        //    LC: QA76.C662.
        //
        //  Parameters:
        //
        //    Input, int N, the number of items in the index vector.
        //
        //    Input, double A[*], the data vector.
        //
        //    Input, int INDX[N], the index vector.
        //
        //    Output, int R8VEC_INDEXED_HEAP_D_MAX, the index in A of the maximum value
        //    in the heap.
        //
    {
        int indx_max = indx[0];

        return indx_max;
    }
}