using System;
using System.Linq;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r82vec_part_quick_a(int n, ref double[] a, int startIndexA, ref int l, ref int r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_PART_QUICK_A reorders an R82 vector as part of a quick sort.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 2, stored as a vector
        //    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
        //
        //    The routine reorders the entries of A.  Using A(1:2,1) as a
        //    key, all entries of A that are less than or equal to the key will
        //    precede the key, which precedes all entries that are greater than the key.
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 8
        //
        //      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
        //
        //    Output:
        //
        //      L = 2, R = 4
        //
        //      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
        //             -----------          ----------------------------------
        //             LEFT          KEY    RIGHT
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of A.
        //
        //    Input/output, double A[N*2].  On input, the array to be checked.
        //    On output, A has been reordered as described above.
        //
        //    Output, int *L, *R, the indices of A that define the three segments.
        //    Let KEY = the input value of A(1:2,1).  Then
        //    I <= L                 A(1:2,I) < KEY;
        //         L < I < R         A(1:2,I) = KEY;
        //                 R <= I    A(1:2,I) > KEY.
        //
    {
        int i;
        int j;
        double[] key = new double[2];
        switch (n)
        {
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R82VEC_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case 1:
                l = 0;
                r = 2;
                return;
        }

        key[0] = a[startIndexA + 2 * 0 + 0];
        key[1] = a[startIndexA + 2 * 0 + 1];
        int m = 1;
        //
        //  The elements of unknown size have indices between L+1 and R-1.
        //
        int ll = 1;
        int rr = n + 1;

        for (i = 2; i <= n; i++)
        {
            if (r8vec_gt(2, a, startIndexA + 2 * ll, key))
            {
                rr -= 1;
                r8vec_swap(2, ref a, ref a, startIndexA + 2 * (rr - 1), startIndexA + 2 * ll);
            }
            else if (r8vec_eq(2, a, key, startIndexA + 2 * ll))
            {
                m += 1;
                r8vec_swap(2, ref a, ref a, startIndexA + 2 * (m - 1), startIndexA + 2 * ll);
                ll += 1;
            }
            else if (r8vec_lt(2, a, key, startIndexA + 2 * ll))
            {
                ll += 1;
            }

        }

        //
        //  Now shift small elements to the left, and KEY elements to center.
        //
        for (i = 0; i < ll - m; i++)
        {
            for (j = 0; j < 2; j++)
            {
                a[startIndexA + 2 * i + j] = a[startIndexA + 2 * (i + m) + j];
            }
        }

        ll -= m;

        for (i = ll; i < ll + m; i++)
        {
            for (j = 0; j < 2; j++)
            {
                a[startIndexA + 2 * i + j] = key[j];
            }
        }

        l = ll;
        r = rr;
    }

    public static void r82vec_permute(int n, ref double[] a, int[] p)

        //****************************************************************************80*
        //
        //  Purpose:
        //
        //    R82VEC_PERMUTE permutes an R82VEC in place.
        //
        //  Discussion:
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
        //    19 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects.
        //
        //    Input/output, double A[2*N], the array to be permuted.
        //
        //    Input, int P[N], the permutation.  P(I) = J means
        //    that the I-th element of the output array should be the J-th
        //    element of the input array.  P must be a legal permutation
        //    of the integers from 1 to N, otherwise the algorithm will
        //    fail catastrophically.
        //
    {
        double[] a_temp = new double[2];
        int i;
        int istart;

        if (!perm_check(n, p))
        {
            Console.WriteLine("");
            Console.WriteLine("R82VEC_PERMUTE - Fatal error!");
            Console.WriteLine("  The input array does not represent");
            Console.WriteLine("  a proper permutation.");
            i4vec_print(n, p, "  The faulty permutation:");
            return;
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
                        a_temp[0] = a[0 + (istart - 1) * 2];
                        a_temp[1] = a[1 + (istart - 1) * 2];
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
                                Console.WriteLine("R82VEC_PERMUTE - Fatal error!");
                                Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                                Console.WriteLine("  an illegal value IGET = " + iget + ".");
                                return;
                            }

                            if (iget == istart)
                            {
                                a[0 + (iput - 1) * 2] = a_temp[0];
                                a[1 + (iput - 1) * 2] = a_temp[1];
                                break;
                            }

                            a[0 + (iput - 1) * 2] = a[0 + (iget - 1) * 2];
                            a[1 + (iput - 1) * 2] = a[1 + (iget - 1) * 2];
                        }
                    }

                    break;
                }
            }
        }

        //
        //  Restore the signs of the entries.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = -p[i];
        }
    }

    public static void r82vec_permute(int n, int[] p, int base_, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_PERMUTE permutes an R82VEC in place.
        //
        //  Discussion:
        //
        //    An R82VEC is a vector whose entries are R82's.
        //    An R82 is a vector of type double precision with two entries.
        //    An R82VEC may be stored as a 2 by N array.
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
        //    30 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects.
        //
        //    Input, int P[N], the permutation.  P(I) = J means
        //    that the I-th element of the output array should be the J-th
        //    element of the input array.
        //
        //    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
        //
        //    Input/output, double A[2*N], the array to be permuted.
        //
    {
        double[] a_temp = new double[2];
        int i;
        int istart;

        if (!perm_check2(n, p, base_))
        {
            Console.WriteLine("");
            Console.WriteLine("R82VEC_PERMUTE - Fatal error!");
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
                        a_temp[0] = a[0 + (istart - 1) * 2];
                        a_temp[1] = a[1 + (istart - 1) * 2];
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
                                Console.WriteLine("R82VEC_PERMUTE - Fatal error!");
                                Console.WriteLine("  Entry IPUT = " + iput + " of the permutation has");
                                Console.WriteLine("  an illegal value IGET = " + iget + ".");
                                return;
                            }

                            if (iget == istart)
                            {
                                a[0 + (iput - 1) * 2] = a_temp[0];
                                a[1 + (iput - 1) * 2] = a_temp[1];
                                break;
                            }

                            a[0 + (iput - 1) * 2] = a[0 + (iget - 1) * 2];
                            a[1 + (iput - 1) * 2] = a[1 + (iget - 1) * 2];
                        }
                    }

                    break;
                }
            }
        }

        //
        //  Restore the signs of the entries.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = -p[i];
        }

        //
        //  Restore the base of the entries.
        //
        for (i = 0; i < n; i++)
        {
            p[i] = p[i] - 1 + base_;
        }
    }

    public static int[] r82vec_sort_heap_index_a(int n, int base_, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
        //
        //  Discussion:
        //
        //    An R82VEC is a vector whose entries are R82's.
        //    An R82 is a vector of type double precision with two entries.
        //    An R82VEC may be stored as a 2 by N array.
        //
        //    The sorting is not actually carried out.  Rather an index array is
        //    created which defines the sorting.  This array may be used to sort
        //    or index the array, or to sort or index related arrays keyed on the
        //    original array.
        //
        //    Once the index array is computed, the sorting can be carried out
        //    "implicitly:
        //
        //      a(*,indx(*))
        //
        //    or explicitly, by the call
        //
        //      r82vec_permute ( n, indx, base, a )
        //
        //    after which a(*,*) is sorted.
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
        //    Input, int BASE, the desired indexing for the sort index:
        //    0 for 0-based indexing,
        //    1 for 1-based indexing.
        //
        //    Input, double A[2*N], an array to be index-sorted.
        //
        //    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
        //    I-th element of the sorted array is A(0:1,R8VEC_SORT_HEAP_INDEX_A(I)).
        //
    {
        double[] aval = new double[2];
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
                indx[0] += base_;
                return indx;
        }

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval[0] = a[0 + indxt * 2];
                aval[1] = a[1 + indxt * 2];
            }
            else
            {
                indxt = indx[ir - 1];
                aval[0] = a[0 + indxt * 2];
                aval[1] = a[1 + indxt * 2];
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
                    if (a[0 + indx[j - 1] * 2] < a[0 + indx[j] * 2] ||
                        Math.Abs(a[0 + indx[j - 1] * 2] - a[0 + indx[j] * 2]) <= double.Epsilon &&
                        a[1 + indx[j - 1] * 2] < a[1 + indx[j] * 2])
                    {
                        j += 1;
                    }
                }

                if (aval[0] < a[0 + indx[j - 1] * 2] ||
                    Math.Abs(aval[0] - a[0 + indx[j - 1] * 2]) <= double.Epsilon &&
                    aval[1] < a[1 + indx[j - 1] * 2])
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

        //
        //  Take care of the base.
        //
        for (i = 0; i < n; i++)
        {
            indx[i] += base_;
        }

        return indx;
    }

    public static void r82vec_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_PRINT prints an R82VEC.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 2, stored as a vector
        //    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, double A[N*2], the vector to be printed.
        //
        //    Input, string TITLE, a title to be printed first.
        //    TITLE may be blank.
        //
    {
        int i;

        if (s_len_trim(title) != 0)
        {
            Console.WriteLine("");
            Console.WriteLine(title + "");
        }

        Console.WriteLine("");
        for (i = 0; i <= n - 1; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(6) + "  "
                                                      + a[2 * i + 0].ToString().PadLeft(14) + "  "
                                                      + a[2 * i + 1].ToString().PadLeft(14) + "");
        }
    }

    public static int[] r82vec_sort_heap_index_a(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
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
        //      A(1:2,INDX(I)), I = 1 to N is sorted,
        //
        //    or explicitly, by the call
        //
        //      call R82VEC_PERMUTE ( N, A, INDX )
        //
        //    after which A(1:2,I), I = 1 to N is sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[2*N], an array to be index-sorted.
        //
        //    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
        //    I-th element of the sorted array is A(0:1,R82VEC_SORT_HEAP_INDEX_A(I-1)).
        //
    {
        double[] aval = new double[2];
        int[] indx;
        switch (n)
        {
            //
            case < 1:
                return null;
            case 1:
                indx = new int[1];
                indx[0] = 1;
                return indx;
        }

        indx = i4vec_indicator(n);

        int l = n / 2 + 1;
        int ir = n;

        for (;;)
        {
            int indxt;
            if (1 < l)
            {
                l -= 1;
                indxt = indx[l - 1];
                aval[0] = a[0 + (indxt - 1) * 2];
                aval[1] = a[1 + (indxt - 1) * 2];
            }
            else
            {
                indxt = indx[ir - 1];
                aval[0] = a[0 + (indxt - 1) * 2];
                aval[1] = a[1 + (indxt - 1) * 2];
                indx[ir - 1] = indx[0];
                ir -= 1;

                if (ir == 1)
                {
                    indx[0] = indxt;
                    break;
                }

            }

            int i = l;
            int j = l + l;

            while (j <= ir)
            {
                if (j < ir)
                {
                    if (a[0 + (indx[j - 1] - 1) * 2] < a[0 + (indx[j] - 1) * 2] ||
                        Math.Abs(a[0 + (indx[j - 1] - 1) * 2] - a[0 + (indx[j] - 1) * 2]) <= double.Epsilon &&
                        a[1 + (indx[j - 1] - 1) * 2] < a[1 + (indx[j] - 1) * 2])
                    {
                        j += 1;
                    }
                }

                if (aval[0] < a[0 + (indx[j - 1] - 1) * 2] ||
                    Math.Abs(aval[0] - a[0 + (indx[j - 1] - 1) * 2]) <= double.Epsilon &&
                    aval[1] < a[1 + (indx[j - 1] - 1) * 2])
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

    public static void r82vec_sort_quick_a(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 2, stored as a vector
        //    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input/output, double A[N*2].
        //    On input, the array to be sorted.
        //    On output, the array has been sorted.
        //
    {
        const int LEVEL_MAX = 25;

        int l_segment = 0;
        int[] rsave = new int[LEVEL_MAX];
        int r_segment = 0;

        switch (n)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R82VEC_SORT_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case 1:
                return;
        }

        int level = 1;
        rsave[level - 1] = n + 1;
        int base_ = 1;
        int n_segment = n;

        while (0 < n_segment)
        {
            //
            //  Partition the segment.
            //
            r82vec_part_quick_a(n_segment, ref a, 2 * (base_ - 1) + 0, ref l_segment, ref r_segment);
            switch (l_segment)
            {
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                case > 1 when LEVEL_MAX < level:
                    Console.WriteLine("");
                    Console.WriteLine("R82VEC_SORT_QUICK_A - Fatal error!");
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
                            if (level <= 1)
                            {
                                n_segment = 0;
                                break;
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

    public static void r82vec_uniform(int n, double[] alo, double[] ahi, ref int seed,
            ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_UNIFORM returns a random R82VEC in a given range.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 2, stored as a vector
        //    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double ALO[2], AHI[2], the range allowed for the entries.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double A[N*2], the vector of randomly chosen integers.
        //
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                a[2 * i + j] = UniformRNG.r8_uniform(alo[j], ahi[j], ref seed);
            }
        }
    }

    public static void r82vec_print_part(int n, double[] a, int max_print, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82VEC_PRINT_PART prints "part" of an R82VEC.
        //
        //  Discussion:
        //
        //    The user specifies MAX_PRINT, the maximum number of lines to print.
        //
        //    If N, the size of the vector, is no more than MAX_PRINT, then
        //    the entire vector is printed, one entry per line.
        //
        //    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
        //    followed by a line of periods suggesting an omission,
        //    and the last entry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries of the vector.
        //
        //    Input, double A[2*N], the vector to be printed.
        //
        //    Input, int MAX_PRINT, the maximum number of lines
        //    to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (max_print)
        {
            case <= 0:
                return;
        }

        switch (n)
        {
            case <= 0:
                return;
        }

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        if (n <= max_print)
        {
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + a[0 + i * 2].ToString().PadLeft(14)
                                       + "  " + a[1 + i * 2].ToString().PadLeft(14) + "");
            }
        }
        else
        {
            switch (max_print)
            {
                case >= 3:
                {
                    for (i = 0; i < max_print - 2; i++)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[0 + i * 2].ToString().PadLeft(14)
                                               + "  " + a[1 + i * 2].ToString().PadLeft(14) + "");
                    }

                    Console.WriteLine("  ........  ..............  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[0 + i * 2].ToString().PadLeft(14)
                                           + "  " + a[1 + i * 2].ToString().PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for (i = 0; i < max_print - 1; i++)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + ": " + a[0 + i * 2].ToString().PadLeft(14)
                                               + "  " + a[1 + i * 2].ToString().PadLeft(14) + "");
                    }

                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + ": " + a[0 + i * 2].ToString().PadLeft(14)
                                           + "  " + a[1 + i * 2].ToString().PadLeft(14)
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }

}