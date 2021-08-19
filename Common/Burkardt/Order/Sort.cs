namespace Burkardt.SortNS
{
    public static partial class Sort
    {
        public class SortRCData
        {
            public int i_save = 0;
            public int j_save = 0;
            public int k_save = 0;
            public int l_save = 0;
            public int n_save = 0;

        }

        public static void sort_rc(ref SortRCData data, int n, ref int indx, ref int i, ref int j, int isgn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SORT_RC externally sorts a list of items into ascending order.
            //
            //  Discussion:
            //
            //    The actual list of data is not passed to the routine.  Hence this
            //    routine may be used to sort integers, reals, numbers, names,
            //    dates, shoe sizes, and so on.  After each call, the routine asks
            //    the user to compare or interchange two items, until a special
            //    return value signals that the sorting is completed.
            //
            //    Note that this function uses internal persistent memory during the sort.
            //
            //  Example:
            //
            //    n = 100;
            //    indx = 0;
            //
            //    for ( ; ; )
            //    {
            //      sort_rc ( n, indx, i, j, isgn );
            //
            //      if ( indx < 0 )
            //      {
            //        isgn = 1;
            //        if ( a[i-1] <= a[j-1] )
            //        {
            //          isgn = -1;
            //        }
            //      }
            //      else if ( 0 < indx )
            //      {
            //        k      = a[i-1];
            //        a[i-1] = a[j-1];
            //        a[j-1] = k;
            //      }
            //      else
            //      {
            //        break;
            //      }
            //    }
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 March 2015
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
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
                data.k_save = n / 2;
                data.l_save = n / 2;
                data.n_save = n;
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

                    data.j_save = data.l_save;
                    data.l_save = data.i_save;
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

                if (data.k_save <= 1)
                {
                    if (data.n_save == 1)
                    {
                        data.i_save = 0;
                        data.j_save = 0;
                        indx = 0;
                    }
                    else
                    {
                        data.i_save = data.n_save;
                        data.j_save = 1;
                        data.n_save = data.n_save - 1;
                        indx = 1;
                    }

                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                data.k_save = data.k_save - 1;
                data.l_save = data.k_save;
            }
            //
            //  0 < INDX: the user was asked to make an interchange.
            //
            else if (indx == 1)
            {
                data.l_save = data.k_save;
            }

            for (;;)
            {

                data.i_save = 2 * data.l_save;

                if (data.i_save == data.n_save)
                {
                    data.j_save = data.l_save;
                    data.l_save = data.i_save;
                    indx = -1;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }
                else if (data.i_save <= data.n_save)
                {
                    data.j_save = data.i_save + 1;
                    indx = -2;
                    i = data.i_save;
                    j = data.j_save;
                    return;
                }

                if (data.k_save <= 1)
                {
                    break;
                }

                data.k_save = data.k_save - 1;
                data.l_save = data.k_save;
            }

            if (data.n_save == 1)
            {
                data.i_save = 0;
                data.j_save = 0;
                indx = 0;
                i = data.i_save;
                j = data.j_save;
            }
            else
            {
                data.i_save = data.n_save;
                data.j_save = 1;
                data.n_save = data.n_save - 1;
                indx = 1;
                i = data.i_save;
                j = data.j_save;
            }

        }

        public static void sort_safe_rc(int n, ref int indx, ref int i, ref int j, int isgn, ref int i_save,
                ref int j_save, ref int k_save, ref int l_save, ref int n_save)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SORT_SAFE_RC externally ascending sorts a list of items.
            //
            //  Discussion:
            //
            //    This is a version of SORT_RC which does not rely on
            //    storing certain work variables internally to the function.  This makes
            //    the function somewhat more awkward to call, but easier to program
            //    in a variety of languages, and safe to use in a parallel programming
            //    environment, or in cases where the sorting of several vectors is to
            //    be carried out at more or less the same time.
            //
            //    The actual list of data is not passed to the routine.  Hence this
            //    routine may be used to sort integers, reals, numbers, names,
            //    dates, shoe sizes, and so on.  After each call, the routine asks
            //    the user to compare or interchange two items, until a special
            //    return value signals that the sorting is completed.
            //
            //  Example:
            //
            //    n = 100;
            //    indx = 0;
            //
            //    for ( ; ; )
            //    {
            //      sort_rc ( n, indx, i, j, isgn, i_save, j_save, k_save, 
            //        l_save, n_save );
            //
            //      if ( indx < 0 )
            //      {
            //        isgn = 1;
            //        if ( a[i-1] <= a[j-1] )
            //        {
            //          isgn = -1;
            //        }
            //      }
            //      else if ( 0 < indx )
            //      {
            //        k      = a[i-1];
            //        a[i-1] = a[j-1];
            //        a[j-1] = k;
            //      }
            //      else
            //      {
            //        break;
            //      }
            //    }
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 March 2015
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
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
            //    Input/output, int &I_SAVE, &J_SAVE, &K_SAVE, &L_SAVE,
            //    &N_SAVE, workspace needed by the routine.  Before calling the function,
            //    the user should declare variables to hold these values, but should
            //    not change them, and need not ever examine them.
            //
        {
            //
            //  INDX = 0: This is the first call.
            //
            if (indx == 0)
            {
                i_save = 0;
                j_save = 0;
                k_save = n / 2;
                l_save = n / 2;
                n_save = n;
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
                        i_save = i_save + 1;
                    }

                    j_save = l_save;
                    l_save = i_save;
                    indx = -1;
                    i = i_save;
                    j = j_save;
                    return;
                }

                if (0 < isgn)
                {
                    indx = 2;
                    i = i_save;
                    j = j_save;
                    return;
                }

                if (k_save <= 1)
                {
                    if (n_save == 1)
                    {
                        i_save = 0;
                        j_save = 0;
                        indx = 0;
                    }
                    else
                    {
                        i_save = n_save;
                        j_save = 1;
                        n_save = n_save - 1;
                        indx = 1;
                    }

                    i = i_save;
                    j = j_save;
                    return;
                }

                k_save = k_save - 1;
                l_save = k_save;
            }
            //
            //  0 < INDX: the user was asked to make an interchange.
            //
            else if (indx == 1)
            {
                l_save = k_save;
            }

            for (;;)
            {
                i_save = 2 * l_save;

                if (i_save == n_save)
                {
                    j_save = l_save;
                    l_save = i_save;
                    indx = -1;
                    i = i_save;
                    j = j_save;
                    return;
                }
                else if (i_save <= n_save)
                {
                    j_save = i_save + 1;
                    indx = -2;
                    i = i_save;
                    j = j_save;
                    return;
                }

                if (k_save <= 1)
                {
                    break;
                }

                k_save = k_save - 1;
                l_save = k_save;
            }

            if (n_save == 1)
            {
                i_save = 0;
                j_save = 0;
                indx = 0;
                i = i_save;
                j = j_save;
            }
            else
            {
                i_save = n_save;
                j_save = 1;
                n_save = n_save - 1;
                indx = 1;
                i = i_save;
                j = j_save;
            }

        }
    }
}