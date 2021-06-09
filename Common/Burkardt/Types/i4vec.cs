using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {


        public static void i4vec_heap_d(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    A heap is an array A with the property that, for every index J,
            //    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
            //    2*J+1 and 2*J+2 are legal).
            //
            //  Diagram:
            //
            //                  A(0)
            //
            //            A(1)         A(2)
            //
            //      A(3)       A(4)  A(5) A(6)
            //
            //    A(7) A(8)  A(9) A(10)
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
            //    Input, int N, the size of the input array.
            //
            //    Input/output, int A[N].
            //    On input, an unsorted array.
            //    On output, the array has been reordered into a heap.
            //
        {
            int i;
            int ifree;
            int key;
            int m;
            //
            //  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
            //
            for (i = (n / 2) - 1; 0 <= i; i--)
            {
                //
                //  Copy the value out of the parent node.
                //  Position IFREE is now "open".
                //
                key = a[i];
                ifree = i;

                for (;;)
                {
                    //
                    //  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
                    //  IFREE.  (One or both may not exist because they equal or exceed N.)
                    //
                    m = 2 * ifree + 1;
                    //
                    //  Does the first position exist?
                    //
                    if (n <= m)
                    {
                        break;
                    }
                    else
                    {
                        //
                        //  Does the second position exist?
                        //
                        if (m + 1 < n)
                        {
                            //
                            //  If both positions exist, take the larger of the two values,
                            //  and update M if necessary.
                            //
                            if (a[m] < a[m + 1])
                            {
                                m = m + 1;
                            }
                        }

                        //
                        //  If the large descendant is larger than KEY, move it up,
                        //  and update IFREE, the location of the free position, and
                        //  consider the descendants of THIS position.
                        //
                        if (key < a[m])
                        {
                            a[ifree] = a[m];
                            ifree = m;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                //
                //  When you have stopped shifting items up, return the item you
                //  pulled out back to the heap.
                //
                a[ifree] = key;
            }

            return;
        }

        public static int[] i4vec_indicator(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Output, int I4VEC_INDICATOR(N), the initialized array.
            //
        {
            int[] a;
            int i;

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = i + 1;
            }

            return a;
        }

        public static void i4vec_print(int n, int[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_PRINT prints an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
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
            //    Input, int A[N], the vector to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + ": " + a[i].ToString().PadLeft(8) + "");
            }
        }

        public static void i4vec_sort_heap_a(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
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
            //    A Nijenhuis and H Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input/output, int A[N].
            //    On input, the array to be sorted;
            //    On output, the array has been sorted.
            //
        {
            int n1;
            int temp;

            if (n <= 1)
            {
                return;
            }

            //
            //  1: Put A into descending heap form.
            //
            i4vec_heap_d(n, ref a);
            //
            //  2: Sort A.
            //
            //  The largest object in the heap is in A[0].
            //  Move it to position A[N-1].
            //
            temp = a[0];
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
                i4vec_heap_d(n1, ref a);
                //
                //  Take the largest object from A[0] and move it to A[N1-1].
                //
                temp = a[0];
                a[0] = a[n1 - 1];
                a[n1 - 1] = temp;

            }

            return;
        }

        public static int i4vec_sorted_unique(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements in A.
            //
            //    Input/output, int A[N].  On input, the sorted
            //    integer array.  On output, the unique elements in A.
            //
            //    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
            //
        {
            int i;
            int unique_num;

            unique_num = 0;

            if (n <= 0)
            {
                return unique_num;
            }

            unique_num = 1;

            for (i = 1; i < n; i++)
            {
                if (a[i] != a[unique_num - 1])
                {
                    unique_num = unique_num + 1;
                    a[unique_num - 1] = a[i];
                }
            }

            return unique_num;
        }

        public static int i4vec2_compare(int n, int[] a1, int[] a2, int i, int j)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_COMPARE compares pairs of integers stored in two I4VECs.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of data items.
            //
            //    Input, int A1[N], A2[N], contain the two components of each item.
            //
            //    Input, int I, J, the items to be compared.  These values will be
            //    1-based indices for the arrays A1 and A2.
            //
            //    Output, int I4VEC2_COMPARE, the results of the comparison:
            //    -1, item I < item J,
            //     0, item I = item J,
            //    +1, item J < item I.
            //
        {
            int isgn;

            isgn = 0;

            if (a1[i - 1] < a1[j - 1])
            {
                isgn = -1;
            }
            else if (a1[i - 1] == a1[j - 1])
            {
                if (a2[i - 1] < a2[j - 1])
                {
                    isgn = -1;
                }
                else if (a2[i - 1] < a2[j - 1])
                {
                    isgn = 0;
                }
                else if (a2[j - 1] < a2[i - 1])
                {
                    isgn = +1;
                }
            }
            else if (a1[j - 1] < a1[i - 1])
            {
                isgn = +1;
            }

            return isgn;
        }

        public static void i4vec2_sort_a(int n, ref int[] a1, ref int[] a2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
            //
            //  Discussion:
            //
            //    Each item to be sorted is a pair of integers (I,J), with the I
            //    and J values stored in separate vectors A1 and A2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items of data.
            //
            //    Input/output, int A1[N], A2[N], the data to be sorted..
            //
        {
            int i;
            int indx;
            int isgn;
            int j;
            int temp;
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
                    isgn = i4vec2_compare(n, a1, a2, i, j);
                }
                else if (indx == 0)
                {
                    break;
                }

            }
        }

        public static void i4vec2_sorted_unique(int n, ref int[] a1, ref int[] a2, ref int nuniq)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC2_SORTED_UNIQUE finds unique elements in a sorted I4VEC2.
            //
            //  Discussion:
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
            //    09 July 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items.
            //
            //    Input/output, int A1[N], A2[N].
            //    On input, the array of N items.
            //    On output, an array of NUNIQ unique items.
            //
            //    Output, int *NUNIQ, the number of unique items.
            //
        {
            int itest;

            nuniq = 0;

            if (n <= 0)
            {
                return;
            }

            nuniq = 1;

            for (itest = 1; itest < n; itest++)
            {
                if (a1[itest] != a1[nuniq - 1] ||
                    a2[itest] != a2[nuniq - 1])
                {
                    a1[nuniq] = a1[itest];
                    a2[nuniq] = a2[itest];
                    nuniq = nuniq + 1;
                }
            }

            return;
        }

        public static void i4vec_dec(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DEC decrements an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, int A[N], the vector to be decremented.
            //
        {
            for (int i = 0; i < n; i++)
            {
                a[i] = a[i] - 1;
            }
        }

        public static void i4vec_inc(int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_INC increments an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, int A[N], the vector to be incremented.
            //
        {
            for (int i = 0; i < n; i++)
            {
                a[i] = a[i] + 1;
            }
        }

        public static void i4vec_data_read(string input_filename, int n, ref int[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_DATA_READ reads data from an I4VEC file.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //    The file is assumed to contain one record per line.
            //
            //    Records beginning with '#' are comments, and are ignored.
            //    Blank lines are also ignored.
            //
            //    Each line that is not ignored is assumed to contain exactly one value.
            //
            //    There are assumed to be exactly (or at least) N such records.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string INPUT_FILENAME, the name of the input file.
            //
            //    Input, int N, the number of points.  The program
            //    will stop reading data once N values have been read.
            //
            //    Output, int TABLE[N], the data.
            //
        {

            try
            {
                string[] input = File.ReadAllLines(input_filename);

                int j = 0;

                foreach (string line in input)
                {
                    if (line[0] == '#' || s_len_trim(line) == 0)
                    {
                        continue;
                    }

                    table[j] = Convert.ToInt32(line);
                    j = j + 1;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_DATA_READ - Fatal error!");
                Console.WriteLine("  Could not open the input file: \"" + input_filename + "\"");
            }

        }

        public static void i4vec_write(string output_filename, int n, int[] table)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_WRITE writes an I4VEC to a file.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 July 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, string OUTPUT_FILENAME, the output filename.
            //
            //    Input, int N, the number of points.
            //
            //    Input, int TABLE[N], the data.
            //
        {
            try
            {
                List<string> output = new List<string>();

                for (int j = 0; j < n; j++)
                {
                    output.Add(table[j].ToString());
                }

                File.WriteAllLines(output_filename, output);
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("I4VEC_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }

        public static int i4vec2_sorted_unique_count(int n, int[] a1, int[] a2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC2_SORTED_UNIQUE_COUNT counts unique elements in an I4VEC2.
        //
        //  Discussion:
        //
        //    Item I is stored as the pair A1(I), A2(I).
        //
        //    The items must have been sorted, or at least it must be the
        //    case that equal items are stored in adjacent vector locations.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of items.
        //
        //    Input, int A1[N], A2[N], the array of N items.
        //
        //    Output, int I4VEC2_SORTED_UNIQUE_COUNT, the number of unique items.
        //
        {
            int i;
            int iu;
            int unique_num;

            unique_num = 0;

            if (n <= 0)
            {
                return unique_num;
            }

            iu = 0;
            unique_num = 1;

            for (i = 1; i < n; i++)
            {
                if (a1[i] != a1[iu] ||
                    a2[i] != a2[iu])
                {
                    iu = i;
                    unique_num = unique_num + 1;
                }
            }

            return unique_num;
        }

        public static void i4vec2_sorted_uniquely(int n1, int[] a1, int[] b1, int n2, int[] a2,
        int[] b2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4VEC2_SORTED_UNIQUELY keeps the unique elements in an I4VEC2.
        //
        //  Discussion:
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
        //    15 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, the number of items.
        //
        //    Input, int A1[N1], B1[N1], the input array.
        //
        //    Input, int N2, the number of unique items.
        //
        //    Input, int A2[N2], B2[N2], the output array of unique items.
        //
        {
            int i1;
            int i2;

            i1 = 0;
            i2 = 0;

            if (n1 <= 0)
            {
                return;
            }

            a2[i2] = a1[i1];
            b2[i2] = b1[i1];

            for (i1 = 1; i1 < n1; i1++)
            {
                if (a1[i1] != a2[i2] || b1[i1] != b2[i2])
                {
                    i2 = i2 + 1;
                    a2[i2] = a1[i1];
                    b2[i2] = b1[i1];
                }
            }
        }

        public static int[] i4vec_copy_new ( int n, int[] a1 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4VEC_COPY_NEW copies an I4VEC.
            //
            //  Discussion:
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 July 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, int A1[N], the vector to be copied.
            //
            //    Output, int I4VEC_COPY_NEW[N], the copy of A1.
            //
        {
            int[] a2 = new int[n];

            for (int i = 0; i < n; i++ )
            {
                a2[i] = a1[i];
            }
            return a2;
        }

    }
}