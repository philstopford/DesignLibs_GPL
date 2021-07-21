using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r83row_max(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83ROW_MAX returns the maximum value in an R83ROW.
            //
            //  Discussion:
            //
            //    An R83ROW is a (3,N) array of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[3*N], the array.
            //
            //    Output, double R83ROW_MAX[3]; the largest entries in each row.
            //
        {
            int DIM_NUM = 3;

            double[] amax;
            int i;
            int j;

            if (n <= 0)
            {
                return null;
            }

            amax = new double[DIM_NUM];

            for (i = 0; i < DIM_NUM; i++)
            {
                amax[i] = a[i + 0 * DIM_NUM];
                for (j = 1; j < n; j++)
                {
                    if (amax[i] < a[i + j * DIM_NUM])
                    {
                        amax[i] = a[i + j * DIM_NUM];
                    }
                }
            }

            return amax;
        }

        public static double[] r83row_min(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83ROW_MIN returns the minimum value in an R83ROW.
            //
            //  Discussion:
            //
            //    An R83ROW is a (3,N) array of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[3*N], the array.
            //
            //    Output, double R83ROW_MIN[3]; the smallest entries in each row.
            //
        {
            int DIM_NUM = 3;

            double[] amin;
            int i;
            int j;

            if (n <= 0)
            {
                return null;
            }

            amin = new double[DIM_NUM];

            for (i = 0; i < DIM_NUM; i++)
            {
                amin[i] = a[i + 0 * DIM_NUM];
                for (j = 1; j < n; j++)
                {
                    if (a[i + j * DIM_NUM] < amin[i])
                    {
                        amin[i] = a[i + j * DIM_NUM];
                    }
                }
            }

            return amin;
        }

        public static void r83row_part_quick_a(int n, ref double[] a, ref int l, ref int r, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83ROW_PART_QUICK_A reorders an R83ROW as part of a quick sort.
            //
            //  Discussion:
            //
            //    An R83ROW is a (3,N) array of R8's.
            //
            //    The routine reorders the entries of A.  Using A(1:3,1) as a
            //    key, all entries of A that are less than or equal to the key will
            //    precede the key, which precedes all entries that are greater than the key.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries of A.
            //
            //    Input/output, double A[3*N].  On input, the array to be checked.
            //    On output, A has been reordered as described above.
            //
            //    Output, int &L, &R, the indices of A that define the three segments.
            //    Let KEY = the input value of A(1:3,1).  Then
            //    I <= L                 A(1:3,I) < KEY;
            //         L < I < R         A(1:3,I) = KEY;
            //                 R <= I    A(1:3,I) > KEY.
            //
        {
            int i;
            int j;
            double[] key = new double[3];
            int ll;
            int m;
            int rr;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("R83ROW_PART_QUICK_A - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
            }

            if (n == 1)
            {
                l = 0;
                r = 2;
                return;
            }

            key[0] = a[index + (3 * 0 + 0)];
            key[1] = a[index + (3 * 0 + 1)];
            key[2] = a[index + (3 * 0 + 2)];
            m = 1;
            //
            //  The elements of unknown size have indices between L+1 and R-1.
            //
            ll = 1;
            rr = n + 1;

            for (i = 2; i <= n; i++)
            {
                if (r8vec_gt(3, a, key, a1Index: index + (+3 * ll)))
                {
                    rr = rr - 1;
                    r8vec_swap(3, ref a, index + (+3 * (rr - 1)), ref a, index + (+3 * ll));
                }
                else if (r8vec_eq(3, a, index + (+3 * ll), key))
                {
                    m = m + 1;
                    r8vec_swap(3, ref a, index + (+3 * (m - 1)), ref a, index + (+3 * ll));
                    ll = ll + 1;
                }
                else if (r8vec_lt(3, a, index + (+3 * ll), key))
                {
                    ll = ll + 1;
                }
            }

            //
            //  Now shift small elements to the left, and KEY elements to center.
            //
            for (i = 0; i < ll - m; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    a[index + (3 * i + j)] = a[index + (3 * (i + m) + j)];
                }
            }

            ll = ll - m;

            for (i = ll; i < ll + m; i++)
            {
                for (j = 0; j < 3; j++)
                {
                    a[index + (3 * i + j)] = key[j];
                }
            }

            l = ll;
            r = rr;
        }

        public static void r83row_print_part(int n, double[] a, int max_print, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83ROW_PRINT_PART prints "part" of an R83ROW.
            //
            //  Discussion:
            //
            //    An R83ROW is a (3,N) array of R8's.
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
            //    11 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries of the vector.
            //
            //    Input, double A[3*N], the vector to be printed.
            //
            //    Input, int MAX_PRINT, the maximum number of lines
            //    to print.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;

            if (max_print <= 0)
            {
                return;
            }

            if (n <= 0)
            {
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
                                           + "  " + a[0 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[1 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[2 + i * 3] + i.ToString().PadLeft(14) + "");
                }
            }
            else if (3 <= max_print)
            {
                for (i = 0; i < max_print - 2; i++)
                {
                    Console.WriteLine("  " + (i + i).ToString().PadLeft(8)
                                           + ": " + a[0 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[1 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[2 + i * 3] + i.ToString().PadLeft(14) + "");
                }

                Console.WriteLine("  ........  ..............  ..............  ..............");
                i = n - 1;
                Console.WriteLine("  " + (i + i).ToString().PadLeft(8)
                                       + ": " + a[0 + i * 3] + i.ToString().PadLeft(14)
                                       + "  " + a[1 + i * 3] + i.ToString().PadLeft(14)
                                       + "  " + a[2 + i * 3] + i.ToString().PadLeft(14) + "");
            }
            else
            {
                for (i = 0; i < max_print - 1; i++)
                {
                    Console.WriteLine("  " + (i + i).ToString().PadLeft(8)
                                           + ": " + a[0 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[1 + i * 3] + i.ToString().PadLeft(14)
                                           + "  " + a[2 + i * 3] + i.ToString().PadLeft(14) + "");
                }

                i = max_print - 1;
                Console.WriteLine("  " + (i + i).ToString().PadLeft(8)
                                       + ": " + a[0 + i * 3] + i.ToString().PadLeft(14)
                                       + "  " + a[1 + i * 3] + i.ToString().PadLeft(14)
                                       + "  " + a[2 + i * 3] + i.ToString().PadLeft(14)
                                       + "  " + "...more entries...");
            }
        }

        public static void r83row_sort_quick_a(int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83ROW_SORT_QUICK_A ascending sorts an R83ROW using quick sort.
            //
            //  Discussion:
            //
            //    An R83ROW is a (3,N) array of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, double A[3*N].
            //    On input, the array to be sorted.
            //    On output, the array has been sorted.
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
                Console.WriteLine("R83ROW_SORT_QUICK_A - Fatal error!");
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
                r83row_part_quick_a(n_segment, ref a, ref l_segment, ref r_segment, index: +3 * (base_ - 1) + 0);
                //
                //  If the left segment has more than one element, we need to partition it.
                //
                if (1 < l_segment)
                {
                    if (LEVEL_MAX < level)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("R83ROW_SORT_QUICK_A - Fatal error!");
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

    }
}