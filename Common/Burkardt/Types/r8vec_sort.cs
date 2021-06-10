using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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
        
    }
}