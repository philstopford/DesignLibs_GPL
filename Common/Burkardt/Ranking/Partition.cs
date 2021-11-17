using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static int[] partition_greedy(int n, ref int[] a)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PARTITION_GREEDY attacks the partition problem with a greedy algorithm.
        // 
        //  Discussion:
        // 
        //    Given a collection of N not necessarily distinct positive integers A(I),
        //    it is desired to partition the values into two groups, whose sums are
        //    as close as possible.
        // 
        //  Algorithm:
        // 
        //    Begin with sets 1 and 2 empty.
        // 
        //    Process the data in descending order of magnitude.
        // 
        //    The next item A(I) is added to set 1 or set 2, whichever has the
        //    smallest current sum.
        // 
        //    Stop as soon as all items have been allocated.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    24 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Brian Hayes,
        //    The Easiest Hard Problem,
        //    American Scientist,
        //    Volume 90, Number 2, March-April 2002, pages 113-117.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of values.  N must be positive.
        // 
        //    Input/output, int A[N], a collection of positive values.
        //    On output, A has been sorted into descending order.
        // 
        //    Output, int PARTITION_GREEDY[N]; an entry is 0 if A[I] is part of
        //    set 0, and 1 if it is assigned to set 1.
        // 
    {
        int i;
        int[] indx;
        int j;
        int[] sums = new int[2];

        sums[0] = 0;
        sums[1] = 0;

        typeMethods.i4vec_sort_insert_d(n, ref a);

        indx = new int[n];

        for (i = 0; i < n; i++)
        {
            if (sums[0] < sums[1])
            {
                j = 0;
            }
            else
            {
                j = 1;
            }

            indx[i] = j;
            sums[j] += a[i];
        }

        return indx;
    }

    public static int partn_enum(int n, int nmax)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PARTN_ENUM enumerates the partitions of N with maximum element NMAX.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    24 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the integer to be partitioned.
        //    Normally N must be positive, but for this routine any
        //    N is allowed.
        // 
        //    Input, int NMAX, the maximum element in the partition.
        //    Normally, 1 <= NMAX <= N is required,
        //    but for this routine any value of NMAX is allowed.
        // 
        //    Output, int PARTN_ENUM is the number of partitions of N
        //    with maximum element NMAX.
        // 
    {
        int[] p;
        int value;

        switch (n)
        {
            case <= 0:
                value = 0;
                break;
            default:
            {
                if (nmax <= 0 || n < nmax)
                {
                    value = 0;
                }
                else
                {
                    p = npart_table(n, nmax);

                    value = p[n + nmax * (n + 1)];
                }

                break;
            }
        }

        return value;
    }

    public static bool partn_sf_check(int n, int nmax, int npart, int[] a)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PARTN_SF_CHECK checks an SF partition of an integer with largest entry NMAX.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NMAX, the value of the largest entry.
        //    1 <= NMAX <= N.
        // 
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input, int A[NPART], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.  The entries must be in DESCENDING order.
        // 
        //    Output, bool PARTN_SF_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        int asum;
        bool check;
        int i;

        check = true;

        switch (n)
        {
            case < 1:
                check = false;
                return check;
        }

        if (nmax < 1 || n < nmax)
        {
            check = false;
            return check;
        }

        if (npart < 1 || n < npart)
        {
            check = false;
            return check;
        }

        // 
        //  Entry 1 must be NMAX.
        // 
        if (a[0] != nmax)
        {
            check = false;
            return check;
        }

        // 
        //  Every entry must lie between 1 and N.
        // 
        for (i = 0; i < npart; i++)
        {
            if (a[i] < 1 || n < a[i])
            {
                check = false;
                return check;
            }
        }

        // 
        //  The entries must be in descending order.
        // 
        for (i = 1; i < npart; i++)
        {
            if (a[i - 1] < a[i])
            {
                check = false;
                return check;
            }
        }

        // 
        //  The entries must add up to N.
        // 
        asum = typeMethods.i4vec_sum(npart, a);

        if (asum != n)
        {
            check = false;
            return check;
        }

        return check;
    }

    public static void partn_successor(int n, int nmax, ref int npart, ref int[] a, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PARTN_SUCCESSOR computes partitions whose largest part is NMAX.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    28 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NMAX, the maximum size of any part of the
        //    partition.  1 <= NMAX <= N.
        // 
        //    Input/output, int NPART, the number of parts of the
        //    partition.  1 <= NPART <= N.
        // 
        //    Input/output, int A[N], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.
        // 
        //    Input/output, int &RANK, the rank.
        //    If RANK = -1 on input, then the routine understands that this is
        //    the first call, and that the user wishes the routine to supply
        //    the first element in the ordering, which has RANK = 0.
        //    In general, the input value of RANK is increased by 1 for output,
        //    unless the very last element of the ordering was input, in which
        //    case the output value of RANK is 0.
        // 
    {
        bool check;
        int i;
        int index;
        int temp;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                npart = n + 1 - nmax;
                a[0] = nmax;
                for (i = 1; i < npart; i++)
                {
                    a[i] = 1;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        check = partn_sf_check(n, nmax, npart, a);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PARTN_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        switch (npart)
        {
            // 
            //  If there are at least two parts, and the next to last is not NMAX,
            //  then rob the last part and pay the next to the last part.
            //  Then, if the next to last part is too big, swap it leftwards.
            // 
            case > 1:
            {
                if (a[npart - 2] < nmax)
                {
                    a[npart - 1] -= 1;
                    a[npart - 2] += 1;
                    index = npart - 1;

                    for (;;)
                    {
                        if (index <= 1)
                        {
                            break;
                        }

                        if (a[index - 1] <= a[index - 2])
                        {
                            break;
                        }

                        temp = a[index - 2];
                        a[index - 2] = a[index - 1];
                        a[index - 1] = temp;

                        index -= 1;
                    }

                    // 
                    //  Sum the tail.
                    // 
                    temp = 0;
                    for (i = index; i < npart; i++)
                    {
                        temp += a[i];
                    }

                    // 
                    //  Spread the sum as 1''s.
                    // 
                    npart = index + temp;
                    for (i = index; i < npart; i++)
                    {
                        a[i] = 1;
                    }

                    rank += 1;
                }

                break;
            }
            // 
            default:
            {
                npart = n + 1 - nmax;
                a[0] = nmax;
                for (i = 1; i < npart; i++)
                {
                    a[i] = 1;
                }

                rank = 0;
                return;
            }
        }
    }

    public static bool setpart_check(int m, int nsub, int[] s, int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SETPART_CHECK checks a set partition.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    31 December 2015
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int M, the number of elements of the set.
        //    M must be positive.
        // 
        //    Input, int NSUB, the number of nonempty subsets into
        //    which the set is partitioned.  1 <= NSUB <= M.
        //
        //    Input, int S[M], contains the integers from 1 to M,
        //    grouped into subsets as described by INDEX.
        // 
        //    Input, int INDEX[NSUB], lists the location in S of the
        //    last element of each subset.  Thus, the elements of subset 1
        //    are S(1) through S(INDEX(1)), the elements of subset 2
        //    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
        // 
        //    Output, bool SETPART_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        bool check;
        int i;
        int imin;
        int j;

        check = true;
        switch (m)
        {
            //
            //  Check M.
            //
            case < 1:
                check = false;
                return check;
        }

        switch (nsub)
        {
            //
            //  Check NSUB.
            //
            case < 1:
                check = false;
                return check;
        }

        // 
        //  Check INDEX.
        // 
        imin = 0;
        for (i = 0; i < nsub; i++)
        {
            if (index[i] <= imin || m < index[i])
            {
                check = false;
                return check;
            }

            imin = index[i];
        }

        // 
        //  Check the elements of S.
        // 
        for (i = 0; i < m; i++)
        {
            if (s[i] <= 0 || m < s[i])
            {
                check = false;
                return check;
            }

            for (j = 0; j < i; j++)
            {
                if (s[j] == s[i])
                {
                    check = false;
                    return check;
                }
            }
        }

        return check;
    }

    public static int setpart_enum(int m)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SETPART_ENUM enumerates the partitions of a set of M elements.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    27 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int M, the number of elements in the set.
        //    M must be positive.  However, for the enumeration routine only,
        //    it is legal to call with any value of M.
        // 
        //    Output, int SETPART_ENUM, the number of partitions of the set.
        // 
    {
        int[] b;
        int i;
        int j;
        int value;

        switch (m)
        {
            case < 0:
                value = 0;
                break;
            case 0:
                value = 1;
                break;
            default:
            {
                b = new int[m + 1];
                b[0] = 1;
                for (j = 1; j <= m; j++)
                {
                    b[j] = 0;
                    for (i = 0; i < j; i++)
                    {
                        b[j] += typeMethods.i4_choose(j - 1, i) * b[i];
                    }
                }

                value = b[m];
                break;
            }
        }

        return value;
    }

    public static int[] setpart_to_rgf(int m, int nsub, int[] s, int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SETPART_TO_RGF converts a set partition to a restricted growth function.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int M, the number of elements of the set.
        //    M must be positive.
        // 
        //    Input, int NSUB, the number of nonempty subsets into
        //    which the set is partitioned.  1 <= NSUB <= M.
        //
        //    Input, int S(M), contains the integers from 1 to M,
        //    grouped into subsets as described by INDEX.
        // 
        //    Input, int INDEX(NSUB), lists the location in S of the
        //    last element of each subset.  Thus, the elements of subset 1
        //    are S(1) through S(INDEX(1)), the elements of subset 2
        //    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
        // 
        //    Output, int SETPART_TO_RGF[M], the restricted growth function from
        //    M to NSUB.
        // 
    {
        bool check;
        int[] f;
        int i;
        int k;
        int khi;
        int klo;
        // 
        //  Check.
        // 
        check = setpart_check(m, nsub, s, index);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("SETPART_TO_RGF - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return null;
        }

        f = new int[m];

        khi = 0;
        for (i = 1; i <= nsub; i++)
        {
            klo = khi + 1;
            khi = index[i - 1];
            for (k = klo; k <= khi; k++)
            {
                f[s[k - 1] - 1] = i;
            }
        }

        return f;
    }

}