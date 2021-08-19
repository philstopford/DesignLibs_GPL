namespace Burkardt.SortNS
{
    public class SortHeapExternalData
    {
        public int i_save = 0;
        public int j_save = 0;
        public int k = 0;
        public int k1 = 0;
        public int n1 = 0;
    }

    public static partial class Sort
    {
        public static void sort_heap_external(ref SortHeapExternalData data, int n, ref int indx, ref int i, ref int j, int isgn )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
        //
        //  Discussion:
        //
        //    The actual list is not passed to the routine.  Hence it may
        //    consist of integers, reals, numbers, names, etc.  The user,
        //    after each return from the routine, will be asked to compare or
        //    interchange two items.
        //
        //    The current version of this code mimics the FORTRAN version,
        //    so the values of I and J, in particular, are FORTRAN indices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt
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
        //    Input, int N, the length of the input list.
        //
        //    Input/output, int &INDX.
        //    The user must set INDX to 0 before the first call.
        //    On return,
        //      if INDX is greater than 0, the user must interchange
        //      items I and J and recall the routine.
        //      If INDX is less than 0, the user is to compare items I
        //      and J and return in ISGN a negative value if I is to
        //      precede J, and a positive value otherwise.
        //      If INDX is 0, the sorting is done.
        //
        //    Output, int &I, &J.  On return with INDX positive,
        //    elements I and J of the user's list should be
        //    interchanged.  On return with INDX negative, elements I
        //    and J are to be compared by the user.
        //
        //    Input, int ISGN. On return with INDX negative, the
        //    user should compare elements I and J of the list.  If
        //    item I is to precede item J, set ISGN negative,
        //    otherwise set ISGN positive.
        //
        {
            //
            //  INDX = 0: This is the first call.
            //
            if (indx == 0)
            {

                data.i_save = 0;
                data.j_save = 0;
                data.k = n / 2;
                data.k1 = data.k;
                data.n1 = n;
            }
            //
            //  INDX < 0: The user is returning the results of a comparison.
            //
            else if (indx < 0)
            {
                if (indx == -2)
                {
                    if (isgn < 0)
                    {
                        data.i_save = data.i_save + 1;
                    }

                    data.j_save = data.k1;
                    data.k1 = data.i_save;
                    indx = -1;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                if (0 < isgn)
                {
                    indx = 2;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                if (data.k <= 1)
                {
                    if (data.n1 == 1)
                    {
                        data.i_save = 0;
                        data.j_save = 0;
                        indx = 0;
                    }
                    else
                    {
                        data.i_save = data.n1;
                        data.j_save = 1;
                        data.n1 = data.n1 - 1;
                        indx = 1;
                    }

                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                data.k = data.k - 1;
                data.k1 = data.k;
            }
            //
            //  0 < INDX: the user was asked to make an interchange.
            //
            else if (indx == 1)
            {
                data.k1 = data.k;
            }

            for (;;)
            {

                data.i_save = 2 * data.k1;

                if (data.i_save == data.n1)
                {
                    data.j_save = data.k1;
                    data.k1 = data.i_save;
                    indx = -1;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }
                else if (data.i_save <= data.n1)
                {
                    data.j_save = data.i_save + 1;
                    indx = -2;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                if (data.k <= 1)
                {
                    break;
                }

                data.k = data.k - 1;
                data.k1 = data.k;
            }

            if (data.n1 == 1)
            {
                data.i_save = 0;
                data.j_save = 0;
                indx = 0;
                i = data.i_save;
                j = data.j_save;
            }
            else
            {
                data.i_save = data.n1;
                data.j_save = 1;
                data.n1 = data.n1 - 1;
                indx = 1;
                i = data.i_save;
                j = data.j_save;
            }
        }
    }
}