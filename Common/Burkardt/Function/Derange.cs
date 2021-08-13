using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Derange
    {
        public static int derange_enum(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM returns the number of derangements of N objects.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation which leaves no object
            //    unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //    D(N) is the number of ways of placing N non-attacking rooks on
            //    an N by N chessboard with one diagonal deleted.
            //
            //    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
            //
            //    The number of permutations with exactly K items in the right
            //    place is COMB(N,K) * D(N-K).
            //
            //    The formula:
            //
            //      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
            //
            //    based on the inclusion/exclusion law.
            //
            //  Recursion:
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(2) = 1
            //      D(N) = (N-1) * ( D(N-1) + D(N-2) )
            //
            //    or
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(N) = N * D(N-1) + (-1)**N
            //
            //  First values:
            //
            //     N         D(N)
            //     0           1
            //     1           0
            //     2           1
            //     3           2
            //     4           9
            //     5          44
            //     6         265
            //     7        1854
            //     8       14833
            //     9      133496
            //    10     1334961
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Output, int DERANGE_ENUM, the number of derangements of N objects.
            //
        {
            int i;
            int value;
            int value1;
            int value2;

            if (n < 0)
            {
                value = 0;
            }
            else if (n == 0)
            {
                value = 1;
            }
            else if (n == 1)
            {
                value = 0;
            }
            else if (n == 2)
            {
                value = 1;
            }
            else
            {
                value1 = 0;
                value = 1;

                for (i = 3; i <= n; i++)
                {
                    value2 = value1;
                    value1 = value;
                    value = (i - 1) * (value1 + value2);
                }
            }

            return value;
        }

        public static void derange_enum2(int n, ref int[] d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM2 returns the number of derangements of 0 through N objects.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation which leaves no object
            //    unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //    D(N) is the number of ways of placing N non-attacking rooks on
            //    an N by N chessboard with one diagonal deleted.
            //
            //    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
            //
            //    The number of permutations with exactly K items in the right
            //    place is COMB(N,K) * D(N-K).
            //
            //    The formula is:
            //
            //      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
            //
            //    based on the inclusion/exclusion law.
            //
            //  Recursion:
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(2) = 1
            //      D(N) = (N-1) * ( D(N-1) + D(N-2) )
            //
            //    or
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(N) = N * D(N-1) + (-1)^N
            //
            //  Example:
            //
            //     N         D(N)
            //     0           1
            //     1           0
            //     2           1
            //     3           2
            //     4           9
            //     5          44
            //     6         265
            //     7        1854
            //     8       14833
            //     9      133496
            //    10     1334961
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the maximum number of objects to be permuted.
            //
            //    Output, int D[N+1]; D(I) is the number of derangements of
            //    I objects.
            //
        {
            int i;

            d[0] = 1;
            d[1] = 0;

            for (i = 2; i <= n; i++)
            {
                d[i] = (i - 1) * (d[i - 1] + d[i - 2]);
            }

        }

        public static int derange_enum3(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM3 returns the number of derangements of 0 through N objects.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation which leaves no object
            //    unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //    D(N) is the number of ways of placing N non-attacking rooks on
            //    an N by N chessboard with one diagonal deleted.
            //
            //    Limit ( N -> Infinity ) D(N)/N! = 1 / e.
            //
            //    The number of permutations with exactly K items in the right
            //    place is COMB(N,K) * D(N-K).
            //
            //    The formula is:
            //
            //      D(N) = N! * ( 1 - 1/1! + 1/2! - 1/3! ... 1/N! )
            //
            //    based on the inclusion/exclusion law.
            //
            //    D(N) = nint ( N! / E )
            //
            //  Recursion:
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(2) = 1
            //      D(N) = (N-1) * ( D(N-1) + D(N-2) )
            //
            //    or
            //
            //      D(0) = 1
            //      D(1) = 0
            //      D(N) = N * D(N-1) + (-1)^N
            //
            //  Example:
            //
            //     N         D(N)
            //     0           1
            //     1           0
            //     2           1
            //     3           2
            //     4           9
            //     5          44
            //     6         265
            //     7        1854
            //     8       14833
            //     9      133496
            //    10     1334961
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the maximum number of objects to be permuted.
            //
            //    Output, int DERANGE_ENUM3, the number of derangements of N objects.
            //
        {
            double E = 2.718281828459045;

            int value;

            if (n < 0)
            {
                value = -1;
            }
            else if (n == 0)
            {
                value = 1;
            }
            else if (n == 1)
            {
                value = 0;
            }
            else
            {
                value = (int)(0.5 + (typeMethods.r8_factorial(n) / E));
            }

            return value;
        }

        public static void derange0_back_candidate(int n, int[] a, int k, ref int nstack, ref int[] stack,
                int[] ncan)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_BACK_CANDIDATE finds possible K-th entries of a derangement.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation of (0,...,N-1) which leaves 
            //    no object unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the derangement.
            //
            //    Input, int A[N].  The first K-1 entries of A
            //    record the currently set values of the derangement.
            //
            //    Input, int K, the entry of the derangement for which candidates
            //    are to be found.
            //
            //    Input/output, int &NSTACK, the length of the stack.
            //
            //    Input/output, int STACK[(N*(N+1))/2].  On output, we have added
            //    the candidates for entry K to the end of the stack.
            //
            //    Input/output, int NCAN[N], the number of candidates for each level.
            //
        {
            int ican;
            int[] ifree;
            int nfree;
            //
            //  Consider all the integers from 1 through N that have not been used yet.
            //
            nfree = n - k + 1;
            ifree = new int[n];

            Permutation.perm0_free(k - 1, a, nfree, ref ifree);
            //
            //  Everything but K is a legitimate candidate for the K-th entry.
            //
            ncan[k - 1] = 0;

            for (ican = 0; ican < nfree; ican++)
            {
                if (ifree[ican] != k - 1)
                {
                    ncan[k - 1] = ncan[k - 1] + 1;
                    stack[nstack] = ifree[ican];
                    nstack = nstack + 1;
                }
            }

        }

        public class DerangeBackData
        {
            public int indx = -1;
            public int k = -1;
            public int[] ncan = null;
            public int[] stack = null;
            public int stack_max = -1;
            public int stack_num = -1;
        }

        public static void derange0_back_next(ref DerangeBackData data, int n, ref int[] a, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_BACK_NEXT returns the next derangement of N items.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation of (0,...,N-1) which leaves 
            //    no object unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //    This routine uses backtracking.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items to be deranged.  N should be 2 or more.
            //
            //    Input/output, int A[N].
            //    On the first call, the input value of A is not important.
            //    On return with MORE = TRUE, A contains the next derangement.
            //    On subsequent input, A should not be changed.
            //
            //    Input/output, bool &MORE.
            //    On first call, set MORE to FALSE and do not alter it after.
            //    On return, MORE is TRUE if another derangement is being returned in A,
            //    and FALSE if no more derangements could be found.
            //
        {
            int i;

            if (!(more))
            {
                if (n < 2)
                {
                    more = false;
                    return;
                }

                data.indx = 0;
                data.k = 0;
                data.stack_max = (n * (n + 1)) / 2;
                data.stack_num = 0;


                data.stack = new int[data.stack_max];
                for (i = 0; i < data.stack_max; i++)
                {
                    data.stack[i] = 0;
                }

                data.ncan = new int[n];
                for (i = 0; i < n; i++)
                {
                    data.ncan[i] = 0;
                }

                more = true;
            }

            for (;;)
            {
                typeMethods.i4vec_backtrack(n, data.stack_max, data.stack, ref a, ref data.indx, ref data.k, ref data.stack_num,
                    ref data.ncan);

                if (data.indx == 1)
                {
                    break;
                }
                else if (data.indx == 2)
                {
                    derange0_back_candidate(n, a, data.k, ref data.stack_num, ref data.stack, data.ncan);
                }
                else
                {
                    more = false;
                    data.ncan = null;
                    data.stack = null;
                    break;
                }
            }
        }

        public static bool derange0_check(int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_CHECK is TRUE if a permutation is a derangement.
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation which leaves no object
            //    unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects permuted.
            //
            //    Input, int A[N], a permutation.
            //
            //    Output, bool DERANGE0_CHECK is TRUE if there was an error.
            //
        {
            bool check;
            int i;
            int j;
            //
            //  Values must be between 0 and N-1.
            //
            for (i = 0; i < n; i++)
            {
                if (a[i] < 0 || n - 1 < a[i])
                {
                    check = false;
                    return check;
                }
            }

            //
            //  Every value must be represented.
            //
            for (j = 0; j < n; j++)
            {
                check = false;
                for (i = 0; i < n; i++)
                {
                    if (a[i] == j)
                    {
                        check = true;
                        break;
                    }
                }

                if (!check)
                {
                    return check;
                }
            }

            //
            //  Values must be deranged.
            //
            for (i = 0; i < n; i++)
            {
                if (a[i] == i)
                {
                    check = false;
                    return check;
                }
            }

            check = true;

            return check;
        }

        public static void derange0_weed_next(int n, int[] a, ref bool more, ref int maxder, ref int numder)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_WEED_NEXT computes derangements of (0,...,N-1).
            //
            //  Discussion:
            //
            //    A derangement of N objects is a permutation which leaves no object
            //    unchanged.
            //
            //    A derangement of N objects is a permutation with no fixed
            //    points.  If we symbolize the permutation operation by "P",
            //    then for a derangment, P(I) is never equal to I.
            //
            //    The number of derangements of N objects is sometimes called
            //    the subfactorial function, or the derangement number D(N).
            //
            //    This routine simply generates all permutations, one at a time,
            //    and weeds out those that are not derangements.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects being permuted.
            //
            //    Input/output, int A[N].
            //    On first call, the input contents of A are unimportant.  But
            //    on the second and later calls, the input value of A should be
            //    the output value returned on the previous call.
            //    On output, A contains the next derangement.
            //
            //    Input/output, bool &MORE.
            //    Set MORE = FALSE before the first call.
            //    MORE will be reset to TRUE and a derangement will be returned.
            //    Each new call produces a new derangement until MORE is returned FALSE.
            //
            //    Input/output, int &MAXDER, &NUMDER, two parameters
            //    used by the program for bookkeeping.  The user should declare these
            //    variables, and pass the output values from one call to the next,
            //    but should not alter them.
            //
        {
            bool deranged;
            //
            //  Initialization on call with MORE = FALSE.
            //
            if (!more)
            {
                maxder = derange_enum(n);
                numder = 0;
            }

            //
            //  Watch out for cases where there are no derangements.
            //
            if (maxder == 0)
            {
                more = false;
                return;
            }

            //
            //  Get the next permutation.
            //
            for (;;)
            {
                Permutation.perm0_lex_next(n, a, ref more);
                //
                //  See if it is a derangment.
                //
                deranged = derange0_check(n, a);

                if (deranged)
                {
                    break;
                }
            }

            numder = numder + 1;

            if (maxder <= numder)
            {
                more = false;
            }

        }
    }
}