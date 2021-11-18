using System;
using System.Globalization;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r82row_max(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_MAX returns the maximum value in an R82ROW.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[2*N], the array.
        //
        //    Output, double R82ROW_MAX[2]; the largest entries in each row.
        //
    {
        const int DIM_NUM = 2;

        int i;

        switch (n)
        {
            case <= 0:
                return null;
        }

        double[] amax = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            amax[i] = a[i + 0 * DIM_NUM];
            int j;
            for (j = 1; j < n; j++)
            {
                if (amax[i] < a[0 + j * DIM_NUM])
                {
                    amax[i] = a[0 + j * DIM_NUM];
                }
            }
        }

        return amax;
    }

    public static double[] r82row_min(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_MIN returns the minimum value in an R82ROW.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the array.
        //
        //    Input, double A[2*N], the array.
        //
        //    Output, double R82ROW_MIN[2]; the smallest entries in each row.
        //
    {
        const int DIM_NUM = 2;

        int i;

        switch (n)
        {
            case <= 0:
                return null;
        }

        double[] amin = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            amin[i] = a[i + 0 * DIM_NUM];
            int j;
            for (j = 1; j < n; j++)
            {
                if (a[0 + j * DIM_NUM] < amin[i])
                {
                    amin[i] = a[0 + j * DIM_NUM];
                }
            }
        }

        return amin;
    }

    public static int r82row_order_type(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_ORDER_TYPE finds if an R82ROW is (non)strictly ascending/descending.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
        //
        //    The dictionary or lexicographic ordering is used.
        //
        //    (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).
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
        //    Input, int N, the number of entries of the array.
        //
        //    Input, double A[2*N], the array to be checked.
        //
        //    Output, int R82ROW_ORDER_TYPE, order indicator:
        //    -1, no discernable order;
        //    0, all entries are equal;
        //    1, ascending order;
        //    2, strictly ascending order;
        //    3, descending order;
        //    4, strictly descending order.
        //
    {
        int order;
        //
        //  Search for the first value not equal to A(1,1).
        //
        int i = 0;

        for (;;)
        {
            i += 1;

            if (n <= i)
            {
                order = 0;
                return order;
            }

            if (a[0 + 0 * 2] < a[0 + i * 2] || Math.Abs(a[0 + 0 * 2] - a[0 + i * 2]) <= double.Epsilon && a[1 + 0 * 2] < a[1 + i * 2])
            {
                order = i switch
                {
                    2 => 2,
                    _ => 1
                };

                break;
            }

            if (!(a[0 + i * 2] < a[0 + 0 * 2]) && (!(Math.Abs(a[0 + i * 2] - a[0 + 0 * 2]) <= double.Epsilon) ||
                                                   !(a[1 + i * 2] < a[1 + 0 * 2])))
            {
                continue;
            }

            order = i switch
            {
                2 => 4,
                _ => 3
            };

            break;
        }

        //
        //  Now we have a "direction".  Examine subsequent entries.
        //
        for (;;)
        {
            i += 1;
            if (n <= i)
            {
                break;
            }

            if (order == 1)
            {
                if (!(a[0 + i * 2] < a[0 + (i - 1) * 2]) &&
                    (!(Math.Abs(a[0 + i * 2] - a[0 + (i - 1) * 2]) <= double.Epsilon) ||
                     !(a[1 + i * 2] < a[1 + (i - 1) * 2])))
                {
                    continue;
                }

                order = -1;
                break;
            }
            else if (order == 2)
            {
                if (a[0 + i * 2] < a[0 + (i - 1) * 2] ||
                    Math.Abs(a[0 + i * 2] - a[0 + (i - 1) * 2]) <= double.Epsilon && a[1 + i * 2] < a[1 + (i - 1) * 2])
                {
                    order = -1;
                    break;
                }

                if (Math.Abs(a[0 + i * 2] - a[0 + (i - 1) * 2]) < double.Epsilon && Math.Abs(a[1 + i * 2] - a[1 + (i - 1) * 2]) <= double.Epsilon)
                {
                    order = 1;
                }
            }
            else if (order == 3)
            {
                if (!(a[0 + (i - 1) * 2] < a[0 + i * 2]) &&
                    (!(Math.Abs(a[0 + (i - 1) * 2] - a[0 + i * 2]) <= double.Epsilon) ||
                     !(a[1 + (i - 1) * 2] < a[1 + i * 2])))
                {
                    continue;
                }

                order = -1;
                break;
            }
            else
            {
                if (a[0 + (i - 1) * 2] < a[0 + i * 2] ||
                    Math.Abs(a[0 + (i - 1) * 2] - a[0 + i * 2]) <= double.Epsilon && a[1 + (i - 1) * 2] < a[1 + i * 2])
                {
                    order = -1;
                    break;
                }

                if (Math.Abs(a[0 + i * 2] - a[0 + (i - 1) * 2]) <= double.Epsilon && Math.Abs(a[1 + i * 2] - a[1 + (i - 1) * 2]) <= double.Epsilon)
                {
                    order = 3;
                }
            }
        }

        return order;
    }

    public static void r82row_part_quick_a(int n, ref double[] a, ref int l, ref int r, int index = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_PART_QUICK_A reorders an R82ROW as part of a quick sort.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
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
        //    Output, int &L, &R, the indices of A that define the three segments.
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
                Console.WriteLine("R82ROW_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            case 1:
                l = 0;
                r = 2;
                return;
        }

        key[0] = a[index + 2 * 0 + 0];
        key[1] = a[index + 2 * 0 + 1];
        int m = 1;
        //
        //  The elements of unknown size have indices between L+1 and R-1.
        //
        int ll = 1;
        int rr = n + 1;

        for (i = 2; i <= n; i++)
        {
            if (r8vec_gt(2, a, key, a1Index: index + + 2 * ll))
            {
                rr -= 1;
                r8vec_swap(2, ref a, ref a, index + + 2 * ll, index + + 2 * (rr - 1));
            }
            else if (r8vec_eq(2, a, key, index + + 2 * ll))
            {
                m += 1;
                r8vec_swap(2, ref a, ref a, index + + 2 * (m - 1), index + + 2 * ll);
                ll += 1;
            }
            else if (r8vec_lt(2, a, key, index + + 2 * ll))
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
                a[index + 2 * i + j] = a[index + 2 * (i + m) + j];
            }
        }

        ll -= m;

        for (i = ll; i < ll + m; i++)
        {
            for (j = 0; j < 2; j++)
            {
                a[index + 2 * i + j] = key[j];
            }
        }

        l = ll;
        r = rr;
    }

    public static void r82row_permute(int n, int[] p, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_PERMUTE permutes an R82ROW in place.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
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
        //    Input/output, double A[2*N], the array to be permuted.
        //
    {
        double[] a_temp = new double[2];
        int i;
        int istart;

        if (!perm0_check(n, p))
        {
            Console.WriteLine("");
            Console.WriteLine("R82ROW_PERMUTE - Fatal error!");
            Console.WriteLine("  PERM0_CHECK rejects permutation.");
            return;
        }

        //
        //  In order for the sign negation trick to work, we need to assume that the
        //  entries of P are strictly positive.  Presumably, the lowest number is 0.
        //  So temporarily add 1 to each entry to force positivity.
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
                                Console.WriteLine("R82ROW_PERMUTE - Fatal error!");
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
        //  Restore the entries.
        //
        for (i = 0; i < n; i++)
        {
            p[i] -= 1;
        }
    }

    public static void r82row_print(int n, double[] a, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_PRINT prints an R82ROW.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components of the vector.
        //
        //    Input, double A[2*N], the vector to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");
        for (j = 0; j < n; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + ": " + a[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + a[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void r82row_print_part(int n, double[] a, int max_print, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_PRINT_PART prints "part" of an R82ROW.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
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
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + a[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + a[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + ": " + a[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                               + "  " + a[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }

                    Console.WriteLine("  ........  ..............  ..............");
                    i = n - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + ": " + a[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + a[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    break;
                }
                default:
                {
                    for (i = 0; i < max_print - 1; i++)
                    {
                        Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                               + ": " + a[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                               + "  " + a[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                    }

                    i = max_print - 1;
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + ": " + a[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + a[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + "...more entries...");
                    break;
                }
            }
        }
    }

    public static int[] r82row_sort_heap_index_a(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82ROW.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
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
        //      r82row_permute ( n, indx, a )
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
        //    Input, double A[2*N], an array to be index-sorted.
        //
        //    Output, int R82ROW_SORT_HEAP_INDEX_A[N], the sort index.  The
        //    I-th element of the sorted array is A(0:1,R82ROW_SORT_HEAP_INDEX_A(I)).
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
                indx[0] = indx[0];
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

        return indx;
    }

    public static void r82row_sort_quick_a(int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R82ROW_SORT_QUICK_A ascending sorts an R82ROW using quick sort.
        //
        //  Discussion:
        //
        //    An R82ROW is a (2,N) array of R8's.
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
        //    Input/output, double A[2*N].
        //    On input, the array to be sorted.
        //    On output, the array has been sorted.
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
                Console.WriteLine("R82ROW_SORT_QUICK_A - Fatal error!");
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
            r82row_part_quick_a(n_segment, ref a, ref l_segment, ref r_segment, index: + 2 * (base_ - 1) + 0);
            switch (l_segment)
            {
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                case > 1 when LEVEL_MAX < level:
                    Console.WriteLine("");
                    Console.WriteLine("R82ROW_SORT_QUICK_A - Fatal error!");
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

}