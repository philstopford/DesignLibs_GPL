using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8vec_bracket ( int n, double[] x, double xval, ref int left,
            ref int right )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
        //
        //  Discussion:
        //
        //    If the values in the vector are thought of as defining intervals
        //    on the real line, then this routine searches for the interval
        //    nearest to or containing the given value.
        //
        //    It is always true that RIGHT = LEFT+1.
        //
        //    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
        //      XVAL   < X[0] < X[1];
        //    If X(1) <= XVAL < X[N-1], then
        //      X[LEFT-1] <= XVAL < X[RIGHT-1];
        //    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
        //      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
        //
        //    For consistency, this routine computes indices RIGHT and LEFT
        //    that are 1-based, although it would be more natural in C and
        //    C++ to use 0-based values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input, double X[N], an array that has been sorted into ascending order.
        //
        //    Input, double XVAL, a value to be bracketed.
        //
        //    Output, int *LEFT, *RIGHT, the results of the search.
        //
    {
        int i;

        for ( i = 2; i <= n - 1; i++ )
        {
            if ( xval < x[i-1] )
            {
                left = i - 1;
                right = i;
                return;
            }

        }

        left = n - 1;
        right = n;
    }
        
    public static void r8vec_bracket(int n, double[] x, int xindex, double xval, ref int left,
            ref int right)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
        //
        //  Discussion:
        //
        //    If the values in the vector are thought of as defining intervals
        //    on the real line, then this routine searches for the interval
        //    nearest to or containing the given value.
        //
        //    It is always true that RIGHT = LEFT+1.
        //
        //    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
        //      XVAL   < X[0] < X[1];
        //    If X(1) <= XVAL < X[N-1], then
        //      X[LEFT-1] <= XVAL < X[RIGHT-1];
        //    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
        //      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
        //
        //    For consistency, this routine computes indices RIGHT and LEFT
        //    that are 1-based, although it would be more natural in C and
        //    C++ to use 0-based values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of input array.
        //
        //    Input, double X[N], an array that has been sorted into ascending order.
        //
        //    Input, double XVAL, a value to be bracketed.
        //
        //    Output, int *LEFT, *RIGHT, the results of the search.
        //
    {
        int i;

        for (i = 2; i <= n - 1; i++)
        {
            if (xval < x[i - 1 + xindex])
            {
                left = i - 1;
                right = i;
                return;
            }

        }

        left = n - 1;
        right = n;
    }

    public static void r8vec_bracket2(int n, double[] x, double xval, int start, ref int left,
            ref int right)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET2 searches a sorted array for successive brackets of a value.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    If the values in the vector are thought of as defining intervals
        //    on the real line, then this routine searches for the interval
        //    containing the given value.
        //
        //    R8VEC_BRACKET2 is a variation on R8VEC_BRACKET.  It seeks to reduce
        //    the search time by allowing the user to suggest an interval that
        //    probably contains the value.  The routine will look in that interval
        //    and the intervals to the immediate left and right.  If this does
        //    not locate the point, a binary search will be carried out on
        //    appropriate subportion of the sorted array.
        //
        //    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
        //    and X(LEFT) <= XVAL <= X(RIGHT).
        //
        //    Special cases:
        //      Value is less than all data values:
        //    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
        //      Value is greater than all data values:
        //    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
        //      Value is equal to a data value:
        //    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of the input array.
        //
        //    Input, double X[N], an array that has been sorted into
        //    ascending order.
        //
        //    Input, double XVAL, a value to be bracketed by entries of X.
        //
        //    Input, int START, between 1 and N, specifies that XVAL
        //    is likely to be in the interval:
        //      [ X(START), X(START+1) ]
        //    or, if not in that interval, then either
        //      [ X(START+1), X(START+2) ]
        //    or
        //      [ X(START-1), X(START) ].
        //
        //    Output, int &LEFT, &RIGHT, the results of the search.
        //
    {
        int high;
        int low;
        switch (n)
        {
            //
            //  Check.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_BRACKET2 - Fatal error!");
                Console.WriteLine("  N < 1.");
                return;
        }

        if (start < 1 || n < start)
        {
            start = (n + 1) / 2;
        }

        //
        //  XVAL = X(START)?
        //
        if (Math.Abs(x[start - 1] - xval) <= double.Epsilon)
        {
            left = start;
            right = start;
            return;
        }
        //
        //  X(START) < XVAL?
        //

        if (x[start - 1] < xval)
        {
            //
            //  X(START) = X(N) < XVAL < oo?
            //
            if (n < start + 1)
            {
                left = start;
                right = -1;
                return;
            }
            //
            //  XVAL = X(START+1)?
            //

            if (Math.Abs(xval - x[start]) <= double.Epsilon)
            {
                left = start + 1;
                right = start + 1;
                return;
            }
            //
            //  X(START) < XVAL < X(START+1)?
            //
            if (xval < x[start])
            {
                left = start;
                right = start + 1;
                return;
            }
            //
            //  X(START+1) = X(N) < XVAL < oo?
            //
            if (n < start + 2)
            {
                left = start + 1;
                right = -1;
                return;
            }
            //
            //  XVAL = X(START+2)?
            //
            if (Math.Abs(xval - x[start + 1]) <= double.Epsilon)
            {
                left = start + 2;
                right = start + 2;
                return;
            }
            //
            //  X(START+1) < XVAL < X(START+2)?
            //
            if (xval < x[start + 1])
            {
                left = start + 1;
                right = start + 2;
                return;
            }
            //
            //  Binary search for XVAL in [ X(START+2), X(N) ],
            //  where XVAL is guaranteed to be greater than X(START+2).
            //
            low = start + 2;
            high = n;

            r8vec_bracket(high + 1 - low, x, +low - 1, xval, ref left, ref right);

            left = left + low - 1;
            right = right + low - 1;
        }
        else
        {
            switch (start)
            {
                //
                //  -oo < XVAL < X(START) = X(1).
                //
                case 1:
                    left = -1;
                    right = start;
                    return;
                //
                default:
                {
                    if (Math.Abs(xval - x[start - 2]) <= double.Epsilon)
                    {
                        left = start - 1;
                        right = start - 1;
                    }
                    //
                    //  X(START-1) < XVAL < X(START)?
                    //
                    else if (x[start - 2] <= xval)
                    {
                        left = start - 1;
                        right = start;
                    }
                    //
                    //  Binary search for XVAL in [ X(1), X(START-1) ],
                    //  where XVAL is guaranteed to be less than X(START-1).
                    //
                    else
                    {
                        low = 1;
                        high = start - 1;
                        r8vec_bracket(high + 1 - low, x, 0, xval, ref left, ref right);
                    }

                    break;
                }
            }
        }
    }

    public static void r8vec_bracket3(int n, double[] t, double tval, ref int left)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET3 finds the interval containing or nearest a given value.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The routine always returns the index LEFT of the sorted array
        //    T with the property that either
        //    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
        //    *  T < T[LEFT] = T[0], or
        //    *  T > T[LEFT+1] = T[N-1].
        //
        //    The routine is useful for interpolation problems, where
        //    the abscissa must be located within an interval of data
        //    abscissas for interpolation, or the "nearest" interval
        //    to the (extreme) abscissa must be found so that extrapolation
        //    can be carried out.
        //
        //    This version of the function has been revised so that the value of
        //    LEFT that is returned uses the 0-based indexing natural to C++.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, length of the input array.
        //
        //    Input, double T[N], an array that has been sorted into ascending order.
        //
        //    Input, double TVAL, a value to be bracketed by entries of T.
        //
        //    Input/output, int *LEFT.
        //    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
        //    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
        //    is searched first, followed by the appropriate interval to the left
        //    or right.  After that, a binary search is used.
        //    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
        //    is the closest to TVAL; it either contains TVAL, or else TVAL
        //    lies outside the interval [ T[0], T[N-1] ].
        //
    {
        int high;
        int low;
        int mid;
        switch (n)
        {
            //  
            //  Check the input data.
            //
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_BRACKET3 - Fatal error//");
                Console.WriteLine("  N must be at least 2.");
                return;
        }

        //
        //  If *LEFT is not between 0 and N-2, set it to the middle value.
        //
        if (left < 0 || n - 2 < left)
        {
            left = (n - 1) / 2;
        }

        //
        //  CASE 1: TVAL < T[*LEFT]:
        //  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
        //
        if (tval < t[left])
        {
            switch (left)
            {
                case 0:
                    return;
                case 1:
                    left = 0;
                    return;
                default:
                {
                    if (t[left - 1] <= tval)
                    {
                        left -= 1;
                        return;
                    }

                    if (tval <= t[1])
                    {
                        left = 0;
                        return;
                    }

                    break;
                }
            }

            // 
            //  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
            //
            low = 1;
            high = left - 2;

            for (;;)
            {
                if (low == high)
                {
                    left = low;
                    return;
                }

                mid = (low + high + 1) / 2;

                if (t[mid] <= tval)
                {
                    low = mid;
                }
                else
                {
                    high = mid - 1;
                }
            }
        }
        // 
        //  CASE 2: T[*LEFT+1] < TVAL:
        //  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
        //

        if (t[left + 1] < tval)
        {
            if (left == n - 2)
            {
                return;
            }

            if (left == n - 3)
            {
                left += 1;
                return;
            }
            if (tval <= t[left + 2])
            {
                left += 1;
                return;
            }
            if (t[n - 2] <= tval)
            {
                left = n - 2;
                return;
            }

            // 
            //  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
            //
            low = left + 2;
            high = n - 3;

            for (;;)
            {

                if (low == high)
                {
                    left = low;
                    return;
                }

                mid = (low + high + 1) / 2;

                if (t[mid] <= tval)
                {
                    low = mid;
                }
                else
                {
                    high = mid - 1;
                }
            }
        }
        //
        //  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
        //  T is just where the user said it might be.
        //
    }


    public static void r8vec_bracket4(int nt, double[] t, int ns, double[] s, ref int[] left)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET4 finds the interval containing or nearest a given value.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    The routine always returns the index LEFT of the sorted array
        //    T with the property that either
        //    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
        //    *  T < T[LEFT] = T[0], or
        //    *  T > T[LEFT+1] = T[NT-1].
        //
        //    The routine is useful for interpolation problems, where
        //    the abscissa must be located within an interval of data
        //    abscissas for interpolation, or the "nearest" interval
        //    to the (extreme) abscissa must be found so that extrapolation
        //    can be carried out.
        //
        //    This version of the function has been revised so that the value of
        //    LEFT that is returned uses the 0-based indexing natural to C++.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, length of the input array.
        //
        //    Input, double T[NT], an array that has been sorted
        //    into ascending order.
        //
        //    Input, int NS, the number of points to be bracketed.
        //
        //    Input, double S[NS], values to be bracketed by entries of T.
        //
        //    Output, int LEFT[NS].
        //    LEFT[I] is set so that the interval [ T[LEFT[I]], T[LEFT[I]+1] ]
        //    is the closest to S[I]; it either contains S[I], or else S[I]
        //    lies outside the interval [ T[0], T[NT-1] ].
        //
    {
        switch (nt)
        {
            //  
            //  Check the input data.
            //
            case < 2:
                Console.WriteLine("");
                Console.WriteLine("R8VEC_BRACKET4 - Fatal error!");
                Console.WriteLine("  NT must be at least 2.");
                return;
        }

        for (int i = 0; i < ns; i++)
        {
            left[i] = (nt - 1) / 2;
            //
            //  CASE 1: S[I] < T[LEFT]:
            //  Search for S[I] in (T[I],T[I+1]), for I = 0 to LEFT-1.
            //
            int high;
            int low;
            int mid;
            if (s[i] < t[left[i]])
            {
                switch (left[i])
                {
                    case 0:
                        continue;
                    case 1:
                        left[i] = 0;
                        continue;
                    default:
                    {
                        if (t[left[i] - 1] <= s[i])
                        {
                            left[i] -= 1;
                            continue;
                        }

                        if (s[i] <= t[1])
                        {
                            left[i] = 0;
                            continue;
                        }

                        break;
                    }
                }

                // 
                //  ...Binary search for S[I] in (T[I],T[I+1]), for I = 1 to *LEFT-2.
                //
                low = 1;
                high = left[i] - 2;

                for (;;)
                {
                    if (low == high)
                    {
                        left[i] = low;
                        break;
                    }

                    mid = (low + high + 1) / 2;

                    if (t[mid] <= s[i])
                    {
                        low = mid;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                }
            }
            // 
            //  CASE 2: T[LEFT+1] < S[I]:
            //  Search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+1 to NT-2.
            //
            else if (t[left[i] + 1] < s[i])
            {
                if (left[i] == nt - 2)
                {
                    continue;
                }

                if (left[i] == nt - 3)
                {
                    left[i] += 1;
                    continue;
                }
                if (s[i] <= t[left[i] + 2])
                {
                    left[i] += 1;
                    continue;
                }
                if (t[nt - 2] <= s[i])
                {
                    left[i] = nt - 2;
                    continue;
                }

                // 
                //  ...Binary search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+2 to NT-3.
                //
                low = left[i] + 2;
                high = nt - 3;

                for (;;)
                {

                    if (low == high)
                    {
                        left[i] = low;
                        break;
                    }

                    mid = (low + high + 1) / 2;

                    if (t[mid] <= s[i])
                    {
                        low = mid;
                    }
                    else
                    {
                        high = mid - 1;
                    }
                }
            }
            //
            //  CASE 3: T[LEFT] <= S[I] <= T[LEFT+1]:
            //
        }
    }


    public static int r8vec_bracket5(int nd, double[] xd, double xi)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.
        //
        //  Discussion:
        //
        //    We assume XD is sorted.
        //
        //    If XI is contained in the interval [XD(1),XD(N)], then the returned 
        //    value B indicates that XI is contained in [ XD(B), XD(B+1) ].
        //
        //    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.
        //
        //    This code implements a version of binary search which is perhaps more
        //    understandable than the usual ones.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data values.
        //
        //    Input, double XD[N], the sorted data.
        //
        //    Input, double XD, the query value.
        //
        //    Output, int R8VEC_BRACKET5, the bracket information.
        //
    {
        int b;

        if (xi < xd[0] || xd[nd - 1] < xi)
        {
            b = -1;
        }
        else
        {
            int l = 0;
            int r = nd - 1;

            while (l + 1 < r)
            {
                int m = (l + r) / 2;
                if (xi < xd[m])
                {
                    r = m;
                }
                else
                {
                    l = m;
                }
            }

            b = l;
        }

        return b;
    }

    public static int[] r8vec_bracket6(int nd, double[] xd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_BRACKET6 brackets data between successive entries of a sorted R8VEC.
        //
        //  Discussion:
        //
        //    We assume XD is sorted.
        //
        //    If XI(I) is contained in the interval [XD(1),XD(N)], then the value of
        //    B(I) indicates that XI(I) is contained in [ XD(B(I)), XD(B(I)+1) ].
        //
        //    If XI(I) is not contained in the interval [XD(1),XD(N)], then B(I) = -1.
        //
        //    This code implements a version of binary search which is perhaps more
        //    understandable than the usual ones.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data values.
        //
        //    Input, double XD[N], the sorted data.
        //
        //    Input, int NI, the number of inquiry values.
        //
        //    Input, double XD[NI], the query values.
        //
        //    Output, int R8VEC_BRACKET6[NI], the bracket information.
        //
    {
        int i;

        int[] b = new int[ni];

        for (i = 0; i < ni; i++)
        {
            if (xi[i] < xd[0] || xd[nd - 1] < xi[i])
            {
                b[i] = -1;
            }
            else
            {
                int l = 0;
                int r = nd - 1;

                while (l + 1 < r)
                {
                    int m = (l + r) / 2;
                    if (xi[i] < xd[m])
                    {
                        r = m;
                    }
                    else
                    {
                        l = m;
                    }
                }

                b[i] = l;
            }
        }

        return b;
    }

}