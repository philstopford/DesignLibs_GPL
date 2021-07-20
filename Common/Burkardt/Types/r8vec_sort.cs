using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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
            int j;
            double x;

            for (i = 1; i < n; i++)
            {
                x = a[i];

                j = i;

                while (1 <= j && x < a[j - 1])
                {
                    a[j] = a[j - 1];
                    j = j - 1;
                }

                a[j] = x;
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
            double aval;
            int i;
            int[] indx;
            int indxt;
            int ir;
            int j;
            int l;

            indx = new int[n];

            for (i = 1; i <= n; i++)
            {
                indx[i - 1] = i;
            }

            l = n / 2 + 1;
            ir = n;

            for (;;)
            {
                if (1 < l)
                {
                    l = l - 1;
                    indxt = indx[l - 1];
                    aval = a[indxt - 1];
                }
                else
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt - 1];
                    indx[ir - 1] = indx[0];
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        for (i = 0; i < n; i++)
                        {
                            indx[i] = indx[i] - 1;
                        }

                        return indx;
                    }

                }

                i = l;
                j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (a[indx[j - 1] - 1] < a[indx[j] - 1])
                        {
                            j = j + 1;
                        }
                    }

                    if (aval < a[indx[j - 1] - 1])
                    {
                        indx[i - 1] = indx[j - 1];
                        i = j;
                        j = j + j;
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
                    l = l - 1;
                    indxt = indx[l - 1];
                    aval = a[indxt - 1];
                }
                else
                {
                    indxt = indx[ir - 1];
                    aval = a[indxt - 1];
                    indx[ir - 1] = indx[0];
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        for (int i = 0; i < n; i++)
                        {
                            indx[i] = indx[i] - 1;
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
                            j = j + 1;
                        }
                    }

                    if (aval < a[indx[j - 1] - 1])
                    {
                        indx[i2 - 1] = indx[j - 1];
                        i2 = j;
                        j = j + j;
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
            int LEVEL_MAX = 30;

            int base_;
            int l_segment = 0;
            int level;
            int n_segment;
            int[] rsave = new int[LEVEL_MAX];
            int r_segment = 0;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R8VEC_SORT_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            }
            else if (n == 1)
            {
                return;
            }

            level = 1;
            rsave[0] = n + 1;
            base_ = 1;
            n_segment = n;

            while (0 < n_segment)
            {
                //
                //  Partition the segment.
                //
                r8vec_part_quick_a(n_segment, ref a, (base_ - 1), ref l_segment, ref r_segment);
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                if (1 < l_segment)
                {

                    if (LEVEL_MAX < level)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R8VEC_SORT_QUICK_A - Fatal error!");
                        Console.WriteLine("  Exceeding recursion maximum of " + LEVEL_MAX + "");
                        return;
                    }

                    level = level + 1;
                    n_segment = l_segment;
                    rsave[level - 1] = r_segment + base_ - 1;
                }
                //
                //  The left segment and the middle segment are sorted.
                //  Must the right segment be partitioned?
                //
                else if (r_segment < n_segment)
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
                            level = level - 1;
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
            int i;
            int indx;
            int isgn;
            int j;
            double temp;
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
            for (;;)
            {
                Helpers.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);
                //
                //  Interchange the I and J objects.
                //
                if (0 < indx)
                {
                    temp = a1[i - 1];
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
                else if (indx == 0)
                {
                    break;
                }
            }

            return;
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
            int i;
            int indx;
            int isgn;
            int j;
            double temp;
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
            for (;;)
            {
                Helpers.sort_heap_external(ref data, n, ref indx, ref i, ref j, isgn);
                //
                //  Interchange the I and J objects.
                //
                if (0 < indx)
                {
                    temp = a1[i - 1];
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
                else if (indx == 0)
                {
                    break;
                }
            }

            return;
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
            int[] indx;
            int indxt;
            int ir;
            int j;
            int l;
            double xval;
            double yval;

            if (n < 1)
            {
                return null;
            }

            indx = new int[n];

            for (i = 0; i < n; i++)
            {
                indx[i] = i;
            }

            if (n == 1)
            {
                indx[0] = indx[0];
                return indx;
            }

            l = n / 2 + 1;
            ir = n;

            for (;;)
            {
                if (1 < l)
                {
                    l = l - 1;
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
                    ir = ir - 1;

                    if (ir == 1)
                    {
                        indx[0] = indxt;
                        break;
                    }
                }

                i = l;
                j = l + l;

                while (j <= ir)
                {
                    if (j < ir)
                    {
                        if (x[indx[j - 1]] < x[indx[j]] ||
                            (x[indx[j - 1]] == x[indx[j]] && y[indx[j - 1]] < y[indx[j]]))
                        {
                            j = j + 1;
                        }
                    }

                    if (xval < x[indx[j - 1]] ||
                        (xval == x[indx[j - 1]] && yval < y[indx[j - 1]]))
                    {
                        indx[i - 1] = indx[j - 1];
                        i = j;
                        j = j + j;
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

            if (n <= 0)
            {
                return;
            }

            unique_num = 1;

            for (itest = 1; itest < n; itest++)
            {
                if (a1[itest] != a1[unique_num - 1] ||
                    a2[itest] != a2[unique_num - 1])
                {
                    a1[unique_num] = a1[itest];
                    a2[unique_num] = a2[itest];
                    unique_num = unique_num + 1;
                }
            }

            return;
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

            if (n <= 0)
            {
                unique_num = 0;
                return;
            }

            i4vec_zeros(n, ref indx);

            unique_num = 1;
            indx[0] = 1;

            for (itest = 2; itest <= n; itest++)
            {
                if (a1[itest - 2] != a1[itest - 1] || a2[itest - 2] != a2[itest - 1])
                {
                    unique_num = unique_num + 1;
                    indx[unique_num - 1] = itest;
                }
            }
        }
    }
}