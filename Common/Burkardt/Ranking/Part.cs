using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static int part_enum(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_ENUM enumerates the number of partitions of N.
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
        //    Output, int PART_ENUM is the number of partitions of N.
        // 
    {
        int value;

        switch (n)
        {
            case < 0:
                value = 0;
                break;
            default:
                int[] p = part_table(n);

                value = p[n];
                break;
        }

        return value;
    }

    public static bool part_rsf_check(int n, int npart, int[] a)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_RSF_CHECK checks a reverse standard form partition of an integer.
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
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input, int A[NPART], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.  The entries must be in ASCENDING order.
        // 
        //    Output, bool PART_RSF_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        int i;

        switch (n)
        {
            case < 1:
                return false;
        }

        if (npart < 1 || n < npart)
        {
            return false;
        }

        // 
        //  Every entry must lie between 1 and N.
        // 
        for (i = 0; i < npart; i++)
        {
            if (a[i] < 1 || n < a[i])
            {
                return false;
            }
        }

        // 
        //  The entries must be in ascending order.
        // 
        for (i = 1; i < npart; i++)
        {
            if (a[i] < a[i - 1])
            {
                return false;
            }
        }

        // 
        //  The entries must add up to N.
        // 
        return typeMethods.i4vec_sum(npart, a) == n;
    }

    public static bool part_sf_check(int n, int npart, int[] a)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_SF_CHECK checks a standard form partition of an integer.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input, int A[NPART], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.  The entries must be in DESCENDING order.
        //
        //    Output, bool PART_SF_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
    {
        int i;

        switch (n)
        {
            case < 1:
                return false;
        }

        if (npart < 1 || n < npart)
        {
            return false;
        }

        // 
        //  Every entry must lie between 1 and N.
        // 
        for (i = 0; i < npart; i++)
        {
            if (a[i] < 1 || n < a[i])
            {
                return false;
            }
        }

        // 
        //  The entries must be in descending order.
        // 
        for (i = 1; i < npart; i++)
        {
            if (a[i - 1] < a[i])
            {
                return false;
            }
        }

        // 
        //  The entries must add up to N.
        // 
        return typeMethods.i4vec_sum(npart, a) == n;
    }

    public static int[] part_sf_conjugate(int n, int npart, int[] a, ref int npart2 )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_SF_CONJUGATE computes the conjugate of a partition.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input, int A[N], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.
        // 
        //    Output, int &NPART2, the number of parts of the conjugate
        //    partition.
        // 
        //    Output, int PART_SF_CONJUGATE[N], contains the conjugate partition.
        // 
    {
        int i;
        // 
        //  Check.
        // 
        bool check = part_sf_check(n, npart, a);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PART_SF_CHECK - Fatal error!");
                Console.WriteLine("  The partition is illegal.");
                return null;
        }

        npart2 = a[0];

        int[] b = new int[n];

        for (i = 0; i < npart2; i++)
        {
            b[i] = 0;
        }

        for (i = 0; i < npart; i++)
        {
            int j;
            for (j = 0; j < a[i]; j++)
            {
                b[j] += 1;
            }
        }

        return b;
    }

    public static int part_sf_majorize(int n, int nparta, int[] a, int npartb, int[] b )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_SF_MAJORIZE determines if partition A majorizes partition B.
        // 
        //  Discussion:
        // 
        //    The partitions must be in standard form.
        // 
        //    If A, with NPARTA parts, and B, with NPARTB parts, are both partitions
        //    of the same positive integer N, then we say that A majorizes B if,
        //    for every index K from 1 to N, it is true that
        // 
        //      sum ( 1 <= I <= K ) B(I) <= sum ( 1 <= I <= K ) A(I)
        // 
        //    where entries of A beyond index NPARTA, and of B beyond BPARTB
        //    are assumed to be 0.  We say that A strictly majorizes B if
        //    A majorizes B, and for at least one index K the inequality is strict.
        // 
        //    For any two partitions of N, it is possible that A majorizes B,
        //    B majorizes A, both partitions majorize each other (in which case
        //    they are equal), or that neither majorizes the other.
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
        //    Jack vanLint, Richard Wilson,
        //    A Course in Combinatorics,
        //    Cambridge, 1992,
        //    ISBN: 0-521-42260-4,
        //    LC: QA164.L56.
        // 
        //  Parameters:
        // 
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NPARTA, the number of parts in partition A.
        //    1 <= NPARTA <= N.
        // 
        //    Input, int A[NPARTA], contains partition A in standard
        //    form.  A(1) through A(NPARTA) contain nonzero integers which sum to N.
        // 
        //    Input, int NPARTB, the number of parts in partition B.
        //    1 <= NPARTB <= N.
        // 
        //    Input, int B[NPARTB], contains partition B in standard
        //    form.  B(1) through B(NPARTB) contain nonzero integers which sum to N.
        // 
        //    Output, int PART_SF_MAJORIZE, the result of the comparison.
        //    -2, A and B are incomparable, but would have been -1.
        //    -1, A < B, (A is strictly majorized by B),
        //     0, A = B, (A and B are identical),
        //    +1, A > B, (A strictly majorizes B),
        //    +2, A and B are incomparable, but would have been +1.
        // 
    {
        int i;
        // 
        //  Check.
        // 
        bool check = part_sf_check(n, nparta, a);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PART_SF_MAJORIZE - Fatal error!");
                Console.WriteLine("  The partition is illegal.");
                return 1;
        }

        check = part_sf_check(n, npartb, b);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PART_SF_MAJORIZE - Fatal error!");
                Console.WriteLine("  The partition is illegal.");
                return 1;
        }

        int result = 0;
        int suma = 0;
        int sumb = 0;

        for (i = 0; i < Math.Min(nparta, npartb); i++)
        {
            if (i < nparta)
            {
                suma += a[i];
            }

            if (i < npartb)
            {
                sumb += b[i];
            }

            switch (result)
            {
                case -1 when sumb < suma:
                    result = -2;
                    return result;
                case 0 when suma < sumb:
                    result = -1;
                    break;
                case 0:
                {
                    if (sumb < suma)
                    {
                        result = +1;
                    }

                    break;
                }
                case +1 when suma < sumb:
                    result = +2;
                    return result;
            }
        }

        return result;
    }

    public static void part_successor(int n, ref int npart, ref int[] a, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_SUCCESSOR computes the lexicographic partition successor.
        // 
        //  Discussion:
        // 
        //    PART_SUCCESSOR is "inspired by" the GenPartitions algorithm,
        //    but instead of relying on recursion, generates the partitions
        //    one at a time.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input/output, int &NPART, the number of parts of the
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
        int i;
        int j;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 1;
                }

                npart = n;
                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        bool check = part_sf_check(n, npart, a);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PART_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The partition is illegal.");
                return;
        }

        // 
        //  If possible, increment the first intermediate position that
        //  is less than its left hand neighbor, and has at least one
        //  right hand neighbor.
        // 
        int ihi = npart - 1;

        for (i = ihi; 2 <= i; i--)
        {
            if (a[i - 1] >= a[i - 2])
            {
                continue;
            }

            int asum = -1;
            for (j = i + 1; j <= npart; j++)
            {
                asum += a[j - 1];
            }

            a[i - 1] += 1;
            for (j = i + 1; j <= npart; j++)
            {
                a[j - 1] = 0;
            }

            npart = i + asum;
            for (j = i + 1; j <= npart; j++)
            {
                a[j - 1] = 1;
            }

            rank += 1;
            return;
        }

        switch (npart)
        {
            // 
            //  A) there are two or more parts
            //  Increment the first, replace the rest by 1''s.
            // 
            case >= 2:
            {
                a[0] += 1;
                for (j = 2; j <= npart; j++)
                {
                    a[j - 1] = 0;
                }

                npart = n - a[0] + 1;
                for (j = 2; j <= npart; j++)
                {
                    a[j - 1] = 1;
                }

                rank += 1;
                break;
            }
            // 
            //  B) there is only one part.
            //  We have reached the last item.
            //  Return the first one.
            // 
            case 1:
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 1;
                }

                npart = n;
                rank = 0;
                break;
            }
        }
    }

    public static int[] part_table(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PART_TABLE tabulates the number of partitions of N.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Output, int P[N+1], P(I) is the number of partitions of I.
        // 
    {
        int i;

        int[] p = new int[n + 1];

        p[0] = 1;

        switch (n)
        {
            case <= 0:
                return p;
        }

        p[1] = 1;

        for (i = 2; i <= n; i++)
        {
            int sign = 1;
            int psum = 0;
            int w = 1;
            int j = 1;
            int wprime = w + j;

            while (w < n)
            {
                switch (i - w)
                {
                    case >= 0 when sign == 1:
                        psum += p[i - w];
                        break;
                    case >= 0:
                        psum -= p[i - w];
                        break;
                }

                if (wprime <= i)
                {
                    switch (sign)
                    {
                        case 1:
                            psum += p[i - wprime];
                            break;
                        default:
                            psum -= p[i - wprime];
                            break;
                    }
                }

                w = w + 3 * j + 1;
                j += 1;
                wprime = w + j;
                sign = -sign;
            }

            p[i] = psum;
        }

        return p;
    }
}