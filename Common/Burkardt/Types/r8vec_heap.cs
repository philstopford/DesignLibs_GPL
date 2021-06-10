namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_heap_a(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_HEAP_A reorders an R8VEC into a ascending heap.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    An ascending heap is an array A with the property that, for every index J,
            //    A[J] <= A[2*J+1] and A[J] <= A[2*J+2], (as long as the indices
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
            //    17 September 2005
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
            //    Input/output, double A[N].
            //    On input, an unsorted array.
            //    On output, the array has been reordered into a heap.
            //
        {
            int i;
            int ifree;
            double key;
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
                            if (a[m + 1] < a[m])
                            {
                                m = m + 1;
                            }
                        }

                        //
                        //  If the large descendant is larger than KEY, move it up,
                        //  and update IFREE, the location of the free position, and
                        //  consider the descendants of THIS position.
                        //
                        if (a[m] <= key)
                        {
                            break;
                        }

                        a[ifree] = a[m];
                        ifree = m;
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

        public static void r8vec_heap_d(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_HEAP_D reorders an R8VEC into a descending heap.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
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
            //    Input/output, double A[N].
            //    On input, an unsorted array.
            //    On output, the array has been reordered into a heap.
            //
        {
            int i;
            int ifree;
            double key;
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
    }
}