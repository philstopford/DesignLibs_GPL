﻿using System;
using System.Linq;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r82vec_part_quick_a(int n, double[] a, ref int l, ref int r )

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
            int ll;
            int m;
            int rr;
            //
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R82VEC_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            }

            if (n == 1)
            {
                l = 0;
                r = 2;
                return;
            }

            key[0] = a[2 * 0 + 0];
            key[1] = a[2 * 0 + 1];
            m = 1;
            //
            //  The elements of unknown size have indices between L+1 and R-1.
            //
            ll = 1;
            rr = n + 1;

            for (i = 2; i <= n; i++)
            {
                if (r8vec_gt(2, a.Skip(2 * ll).ToArray(), key))
                {
                    rr = rr - 1;
                    double[] tmp = a.Skip(2 * (rr - 1)).ToArray();
                    double[] tmp2 = a.Skip(2 * ll).ToArray();
                    r8vec_swap(2, ref tmp, ref tmp2);
                    for (int t = 0; t < tmp.Length; t++)
                    {
                        a[(2 * (rr - 1)) + t] = tmp[t];
                    }
                    for (int t = 0; t < tmp2.Length; t++)
                    {
                        a[(2 * ll) + t] = tmp2[t];
                    }
                }
                else if (r8vec_eq(2, a.Skip(2 * ll).ToArray(), key))
                {
                    m = m + 1;
                    double[] tmp = a.Skip(2 * (m - 1)).ToArray();
                    double[] tmp2 = a.Skip(2 * ll).ToArray();
                    r8vec_swap(2, ref tmp, ref tmp2);
                    for (int t = 0; t < tmp.Length; t++)
                    {
                        a[(2 * (m - 1)) + t] = tmp[t];
                    }
                    for (int t = 0; t < tmp2.Length; t++)
                    {
                        a[(2 * ll) + t] = tmp2[t];
                    }
                    ll = ll + 1;
                }
                else if (r8vec_lt(2, a.Skip(2 * ll).ToArray(), key))
                {
                    ll = ll + 1;
                }

            }

            //
            //  Now shift small elements to the left, and KEY elements to center.
            //
            for (i = 0; i < ll - m; i++)
            {
                for (j = 0; j < 2; j++)
                {
                    a[2 * i + j] = a[2 * (i + m) + j];
                }
            }

            ll = ll - m;

            for (i = ll; i < ll + m; i++)
            {
                for (j = 0; j < 2; j++)
                {
                    a[2 * i + j] = key[j];
                }
            }

            l = ll;
            r = rr;

            return;
        }

        public static void r82vec_permute(int n, double[] a, int[] p )

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
            int iget;
            int iput;
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
                if (p[istart - 1] < 0)
                {
                    continue;
                }
                else if (p[istart - 1] == istart)
                {
                    p[istart - 1] = -p[istart - 1];
                    continue;
                }
                else
                {
                    a_temp[0] = a[0 + (istart - 1) * 2];
                    a_temp[1] = a[1 + (istart - 1) * 2];
                    iget = istart;
                    //
                    //  Copy the new value into the vacated entry.
                    //
                    for (;;)
                    {
                        iput = iget;
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
            }

            //
            //  Restore the signs of the entries.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = -p[i];
            }

            return;
        }

        public static void r82vec_print(int n, double[] a, string title )

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
            int i;
            int[] indx;
            int indxt;
            int ir;
            int j;
            int l;
            //
            if (n < 1)
            {
                return null;
            }

            if (n == 1)
            {
                indx = new int[1];
                indx[0] = 1;
                return indx;
            }

            indx = i4vec_indicator(n);

            l = n / 2 + 1;
            ir = n;

            for (;;)
            {
                if (1 < l)
                {
                    l = l - 1;
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
                        if (a[0 + (indx[j - 1] - 1) * 2] < a[0 + (indx[j] - 1) * 2] ||
                            (a[0 + (indx[j - 1] - 1) * 2] == a[0 + (indx[j] - 1) * 2] &&
                             a[1 + (indx[j - 1] - 1) * 2] < a[1 + (indx[j] - 1) * 2]))
                        {
                            j = j + 1;
                        }
                    }

                    if (aval[0] < a[0 + (indx[j - 1] - 1) * 2] ||
                        (aval[0] == a[0 + (indx[j - 1] - 1) * 2] &&
                         aval[1] < a[1 + (indx[j - 1] - 1) * 2]))
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
            int LEVEL_MAX = 25;

            int base_;
            int l_segment = 0;
            int level;
            int n_segment;
            int[] rsave = new int[LEVEL_MAX];
            int r_segment = 0;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R82VEC_SORT_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            }

            if (n == 1)
            {
                return;
            }

            level = 1;
            rsave[level - 1] = n + 1;
            base_ = 1;
            n_segment = n;

            while (0 < n_segment)
            {
                //
                //  Partition the segment.
                //
                r82vec_part_quick_a(n_segment, a.Skip(2 * (base_ - 1) + 0).ToArray(), ref l_segment, ref r_segment);
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                if (1 < l_segment)
                {
                    if (LEVEL_MAX < level)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R82VEC_SORT_QUICK_A - Fatal error!");
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
                        if (level <= 1)
                        {
                            n_segment = 0;
                            break;
                        }

                        base_ = rsave[level - 1];
                        n_segment = rsave[level - 2] - rsave[level - 1];
                        level = level - 1;

                        if (0 < n_segment)
                        {
                            break;
                        }

                    }

                }

            }
        }

        public static void r82vec_uniform(int n, double[] alo, double[] ahi, ref int seed,
        ref double[] a )

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
    }
}