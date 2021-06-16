using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public class VecGrayData
        {
            public int[] active;
            public int[] dir;
        }
        
        public static void vec_next_gray(ref VecGrayData data, int n, int[] base_, ref int[] a, ref bool done, ref int change )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_NEXT_GRAY computes the elements of a product space.
        //
        //  Discussion:
        //
        //    The elements are produced one at a time.
        //
        //    This routine handles the case where the number of degrees of freedom may
        //    differ from one component to the next.
        //
        //    A method similar to the Gray code is used, so that successive
        //    elements returned by this routine differ by only a single element.
        //
        //    The routine uses internal static memory.
        //
        //  Example:
        //
        //    N = 2, BASE = ( 2, 3 ), DONE = TRUE
        //
        //     A    DONE  CHANGE
        //    ---  -----  ------
        //    0 0  FALSE    1
        //    0 1  FALSE    2
        //    0 2  FALSE    2
        //    1 2  FALSE    1
        //    1 1  FALSE    2
        //    1 0  FALSE    2
        //    1 0   TRUE   -1  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472.
        //
        //  Parameters:
        //
        //    Input, int N, the number of components.
        //
        //    Input, int BASE[N], contains the number of degrees of
        //    freedom of each component.  The output values of A will
        //    satisfy 0 <= A[I] < BASE[I].
        //
        //    Input/output, int A[N].  On the first call, the input value
        //    of A doesn't matter.  Thereafter, it should be the same as
        //    its output value from the previous call.  On output, if DONE
        //    is FALSE, then A contains the next element of the space.
        //
        //    Input/output, bool *DONE.  On the first call, the user must
        //    set DONE to TRUE.  This signals the program to initialize data.
        //    On every return, if DONE is FALSE, the program has computed
        //    another entry, which is contained in A.  If DONE is TRUE,
        //    then there are no more entries, and the program should not be
        //    called for any more.
        //
        //    Output, int *CHANGE, is set to the index of the element whose
        //    value was changed.  On return from the first call, CHANGE
        //    is 1, even though all the elements have been "changed".  On
        //    return with DONE equal to TRUE, CHANGE is -1.
        //
        {
            int i;
            //
            //  The user is calling for the first time.
            //
            if (done)
            {
                done = false;
                for (i = 0; i < n; i++)
                {
                    a[i] = 0;
                }

                data.active = new int[n];
                data.dir = new int[n];

                for (i = 0; i < n; i++)
                {
                    data.dir[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    data.active[i] = 1;
                }

                for (i = 0; i < n; i++)
                {
                    if (base_[i] < 1)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("VEC_NEXT_GRAY - Warning");
                        Console.WriteLine("  For index I = " + i + "");
                        Console.WriteLine("  the nonpositive value of BASE[I] = " + base_[i] + "");
                        Console.WriteLine("  which was reset to 1!");
                        base_[i] = 1;
                        data.active[i] = 0;
                    }
                    else if (base_[i] == 1)
                    {
                        data.active[i] = 0;
                    }

                }

                change = 0;
                return;

            }

            //
            //  Find the maximum active index.
            //
            change = -1;

            for (i = 0; i < n; i++)
            {
                if (data.active[i] != 0)
                {
                    change = i;
                }
            }

            //
            //  If there are NO active indices, we have generated all vectors.
            //
            if (change == -1)
            {
                data.dir = null;
                data.active = null;
                done = true;
                return;
            }

            //
            //  Increment the element with maximum active index.
            //
            a[change] = a[change] + data.dir[change];
            //
            //  If we attained a minimum or maximum value, reverse the direction
            //  vector, and deactivate the index.
            //
            if (a[change] == 0 || a[change] == base_[change] - 1)
            {
                data.dir[change] = -data.dir[change];
                data.active[change] = 0;
            }

            //
            //  Activate all subsequent indices.
            //
            for (i = change + 1; i < n; i++)
            {
                if (1 < base_[i])
                {
                    data.active[i] = 1;
                }
            }
        }

        public static void vector_constrained_next4(int n, double[] alpha, int[] x_min,
        int[] x_max, ref int[] x, double q, ref bool more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VECTOR_CONSTRAINED_NEXT4 returns the "next" constrained vector.
        //
        //  Discussion:
        //
        //    This routine is similar to VECTOR_CONSTRAINED2 and VECTOR_CONSTRAINED3.
        //
        //    We consider all vectors X of dimension N whose components
        //    satisfy X_MIN(1:N) <= X(1:N) <= X_MAX(1:N).
        //
        //    We are only interested in the subset of these vectors which
        //    satisfy the following constraint:
        //
        //      sum ( 1 <= I <= N ) ( ALPHA(I) * X(I) ) <= Q
        //
        //  Example:
        //
        //    N = 3
        //    ALPHA    4.0  3.0  5.0
        //    Q       20.0
        //    X_MIN:   1   1   1
        //    X_MAX:   5   6   4
        //
        //    P = 120
        //
        //    #  X(1)  X(2)  X(3)      Total
        //
        //    1    1     1     1       12.0
        //    2    2     1     1       20.0
        //    3    1     2     1       15.0
        //    4    2     2     1       19.0
        //    5    1     3     1       18.0
        //    6    1     1     2       17.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of components in the vector.
        //
        //    Input, double ALPHA[N], the coefficient vector.
        //
        //    Input, int X_MIN[N], X_MAX[N], the minimum and maximum
        //    values allowed in each component.
        //
        //    Input/output, integer X[N].  On first call (with MORE = FALSE),
        //    the input value of X is not important.  On subsequent calls, the
        //    input value of X should be the output value from the previous call.
        //    On output, (with MORE = TRUE), the value of X will be the "next"
        //    vector in the reverse lexicographical list of vectors that satisfy
        //    the condition.  However, on output with MORE = FALSE, the vector
        //    X is meaningless, because there are no more vectors in the list.
        //
        //    Input, double Q, the limit on the sum.
        //
        //    Input/output, bool *MORE.  On input, if the user has set MORE
        //    FALSE, the user is requesting the initiation of a new sequence
        //    of values.  If MORE is TRUE, then the user is requesting "more"
        //    values in the current sequence.  On output, if MORE is TRUE,
        //    then another value was found and returned in X, but if MORE is
        //    FALSE, then there are no more values in the sequence, and X is
        //    NOT the next value.
        //
        {
            int i;
            int j;
            double total;

            if (!(more))
            {
                for (j = 0; j < n; j++)
                {
                    x[j] = x_min[j];
                }

                total = 0.0;
                for (j = 0; j < n; j++)
                {
                    total = total + alpha[j] * (double) x[j];
                }

                if (q < total)
                {
                    more = false;
                }
                else
                {
                    more = true;
                }

                return;
            }
            else
            {
                i = 0;

                for (;;)
                {
                    if (x[i] < x_max[i])
                    {
                        x[i] = x[i] + 1;

                        total = 0;
                        for (j = 0; j < n; j++)
                        {
                            total = total + alpha[j] * (double) x[j];
                        }

                        if (total <= q)
                        {
                            break;
                        }
                    }

                    x[i] = x_min[i];

                    i = i + 1;

                    if (n <= i)
                    {
                        more = false;
                        break;
                    }
                }
            }
        }
    }
}