namespace Burkardt.SortNS;

public class SortHeapExternalData
{
    public int i_save;
    public int j_save;
    public int k;
    public int k1;
    public int n1;
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
        switch (indx)
        {
            //
            //  INDX = 0: This is the first call.
            //
            case 0:
                data.i_save = 0;
                data.j_save = 0;
                data.k = n / 2;
                data.k1 = data.k;
                data.n1 = n;
                break;
            //
            //  INDX < 0: The user is returning the results of a comparison.
            //
            case < 0 when indx == -2:
            {
                switch (isgn)
                {
                    case < 0:
                        data.i_save += 1;
                        break;
                }

                data.j_save = data.k1;
                data.k1 = data.i_save;
                indx = -1;
                i = data.i_save;
                j = data.j_save;
                return;
            }
            case < 0 when 0 < isgn:
                indx = 2;
                i = data.i_save;
                j = data.j_save;
                return;
            case < 0 when data.k <= 1:
            {
                switch (data.n1)
                {
                    case 1:
                        data.i_save = 0;
                        data.j_save = 0;
                        indx = 0;
                        break;
                    default:
                        data.i_save = data.n1;
                        data.j_save = 1;
                        data.n1 -= 1;
                        indx = 1;
                        break;
                }

                i = data.i_save;
                j = data.j_save;
                return;
            }
            case < 0:
                data.k -= 1;
                data.k1 = data.k;
                break;
            //
            //  0 < INDX: the user was asked to make an interchange.
            //
            case 1:
                data.k1 = data.k;
                break;
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

            if (data.i_save <= data.n1)
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

            data.k -= 1;
            data.k1 = data.k;
        }

        switch (data.n1)
        {
            case 1:
                data.i_save = 0;
                data.j_save = 0;
                indx = 0;
                i = data.i_save;
                j = data.j_save;
                break;
            default:
                data.i_save = data.n1;
                data.j_save = 1;
                data.n1 -= 1;
                indx = 1;
                i = data.i_save;
                j = data.j_save;
                break;
        }
    }
}