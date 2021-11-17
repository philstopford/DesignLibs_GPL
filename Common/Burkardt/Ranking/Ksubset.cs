using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool ksubset_colex_check(int k, int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_COLEX_CHECK checks a K subset in colex form.
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
        //    Input, int K, the number of elements each K subset must
        //    have. 0 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    0 <= N.
        // 
        //    Input, int T[K], describes a K subset.  T(I) is the I-th
        //    element of the K subset.  The elements must be listed in
        //    DESCENDING order.
        // 
        //    Output, bool KSUBSET_COLEX_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        bool check;
        int i;
        int tmax;

        check = true;

        switch (n)
        {
            case < 0:
                check = false;
                return check;
        }

        if (k < 0 || n < k)
        {
            check = false;
            return check;
        }

        tmax = n + 1;

        for (i = 0; i < k; i++)
        {
            if (t[i] <= 0 || tmax <= t[i])
            {
                check = false;
                return check;
            }

            tmax = t[i];
        }

        return check;
        ;
    }

    public static int ksubset_colex_rank(int k, int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_COLEX_RANK computes the colex rank of a K subset.
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
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input, int T[K[, describes a K subset.  T(I) is the I-th
        //    element of the K subset.  The elements must be listed in DESCENDING order.
        // 
        //    Output, int KSUBSET_COLEX_RANK, the rank of the subset.
        // 
    {
        bool check;
        int i;
        int rank;
        // 
        //  Check.
        // 
        check = ksubset_colex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_COLEX_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return 1;
        }

        rank = 0;

        for (i = 0; i < k; i++)
        {
            rank += typeMethods.i4_choose(t[i] - 1, k - i);
        }

        return rank;
    }

    public static void ksubset_colex_successor(int k, int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_COLEX_SUCCESSOR computes the K subset colex successor.
        // 
        //  Discussion:
        // 
        //    In the original code, there is a last element with no successor.
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
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input/output, int T[K], describes a K subset.  T(I) is the
        //    I-th element.  The elements must be listed in DESCENDING order.
        //    On input, T describes a K subset.
        //    On output, T describes the next K subset in the ordering.
        //    If the input T was the last in the ordering, then the output T
        //    will be the first.
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
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                for (i = 1; i <= k; i++)
                {
                    t[i - 1] = k + 1 - i;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        check = ksubset_colex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_COLEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        for (i = k - 1; 1 <= i; i--)
        {
            if (t[k - i] + 1 < t[k - i - 1])
            {
                t[k - i] += 1;
                rank += 1;
                return;
            }
        }

        if (t[0] < n)
        {
            t[0] += 1;
            for (i = 1; i <= k - 1; i++)
            {
                t[k - i] = i;
            }

            rank += 1;
            return;
        }

        // 
        //  The last K subset was input.
        //  Return the first one.
        // 
        for (i = 1; i <= k; i++)
        {
            t[i - 1] = k + 1 - i;
        }

        rank = 0;
    }

    public static int[] ksubset_colex_unrank(int rank, int k, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_COLEX_UNRANK computes the K subset of given colex rank.
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
        //    Input, int RANK, the rank of the K subset.
        // 
        //    Input, int K, the number of elements each K subset must
        //    have.  0 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Output, int KSUBSET_COLEX_UNRANK[K], describes the K subset of the given
        //    rank.  T(I) is the I-th element.  The elements must be listed in
        //    DESCENDING order.
        // 
    {
        int i;
        int nksub;
        int rank_copy;
        int[] t;
        int x;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_COLEX_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
        }

        switch (k)
        {
            case 0:
                t = new int[k];
                return t;
        }

        if (k < 0 || n < k)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_COLEX_UNRANK - Fatal error!");
            Console.WriteLine("  Input K is illegal.");
            return null;
        }

        nksub = ksubset_enum(k, n);

        if (rank < 0 || nksub < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_COLEX_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        // 
        rank_copy = rank;

        x = n;

        t = new int[k];

        for (i = 1; i <= k; i++)
        {
            while (rank_copy < typeMethods.i4_choose(x, k + 1 - i))
            {
                x -= 1;
            }

            t[i - 1] = x + 1;
            rank_copy -= typeMethods.i4_choose(x, k + 1 - i);
        }

        return t;
    }

    public static int ksubset_enum(int k, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_ENUM enumerates the K element subsets of an N set.
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
        //  Parameters:
        // 
        //    Input, int K, the number of elements each K subset must
        //    have. 0 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    0 <= N.
        // 
        //    Output, int KSUBSET_ENUM, the number of distinct elements.
        // 
    {
        int value = typeMethods.i4_choose(n, k);

        return value;
    }

    public static bool ksubset_lex_check(int k, int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_LEX_CHECK checks a K subset in lex form.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 November 2015
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
        //    Input, int K, the number of elements each K subset must
        //    have. 0 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    0 <= N.
        // 
        //    Input, int T[K], describes a K subset.  T(I) is the I-th
        //    element of the K subset.  The elements must be listed in
        //    ASCENDING order.
        // 
        //    Output, bool KSUBSET_LEX_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, data is not legal.
        // 
    {
        bool check;
        int i;
        int tmin;

        check = true;

        switch (n)
        {
            case < 0:
                check = false;
                return check;
        }

        if (k < 0 || n < k)
        {
            check = false;
            return check;
        }

        tmin = 0;

        for (i = 0; i < k; i++)
        {
            if (t[i] <= tmin || n < t[i])
            {
                check = false;
                return check;
            }

            tmin = t[i];
        }

        return check;
    }

    public static int ksubset_lex_rank(int k, int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_LEX_RANK computes the lexicographic rank of a K subset.
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
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input, int T[K], describes a K subset.  T(I) is the I-th
        //    element.  The elements must be listed in ascending order.
        // 
        //    Output, int KSUBSET_LEX_RANK, the rank of the K subset.
        // 
    {
        bool check;
        int i;
        int j;
        int rank;
        int tim1;
        // 
        //  Check.
        // 
        check = ksubset_lex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_LEX_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return 0;
        }

        rank = 0;

        for (i = 1; i <= k; i++)
        {
            tim1 = i switch
            {
                1 => 0,
                _ => t[i - 2]
            };

            if (tim1 + 1 <= t[i - 1] - 1)
            {
                for (j = tim1 + 1; j <= t[i - 1] - 1; j++)
                {
                    rank += typeMethods.i4_choose(n - j, k - i);
                }
            }
        }

        return rank;
    }

    public static void ksubset_lex_successor(int k, int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_LEX_SUCCESSOR computes the K subset lexicographic successor.
        // 
        //  Discussion:
        // 
        //    In the original code, there is a last element with no successor.
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
        //    Input, int K, the number of elements each K subset must
        //    have. 1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input/output, int T[K], describes a K subset.  T(I) is
        //    the I-th element.  The elements must be listed in ascending order.
        //    On input, T describes a K subset.
        //    On output, T describes the next K subset in the ordering.
        //    If the input T was the last in the ordering, then the output T
        //    will be the first.
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
        int isave;
        int j;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
                typeMethods.i4vec_indicator1(k, ref t);
                rank = 0;
                return;
        }

        // 
        //  Check.
        // 
        check = ksubset_lex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_LEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        isave = 0;

        for (i = k; 1 <= i; i--)
        {
            if (t[i - 1] != n - k + i)
            {
                isave = i;
                break;
            }
        }

        switch (isave)
        {
            // 
            //  The last K subset was input.
            //  Return the first one.
            // 
            case 0:
                typeMethods.i4vec_indicator1(k, ref t);
                rank = 0;
                break;
            default:
            {
                for (j = k; isave <= j; j--)
                {
                    t[j - 1] = t[isave - 1] + 1 + j - isave;
                }

                rank += 1;
                break;
            }
        }
    }

    public static int[] ksubset_lex_unrank(int rank, int k, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_LEX_UNRANK computes the K subset of given lexicographic rank.
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
        //    Input, int RANK, the rank of the K subset.
        // 
        //    Input, int K, the number of elements each K subset must
        //    have.  0 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Output, int KSUBSET_LEX_RANK[K], describes the K subset of the given
        //    rank.  T(I) is the I-th element.  The elements must be listed in
        //    ascending order.
        // 
    {
        int i;
        int nksub;
        int rank_copy;
        int[] t;
        int x;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
        }

        switch (k)
        {
            case 0:
                t = new int[k];
                return t;
        }

        if (k < 0 || n < k)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_LEX_UNRANK - Fatal error!");
            Console.WriteLine("  Input K is illegal.");
            return null;
        }

        nksub = ksubset_enum(k, n);

        if (rank < 0 || nksub < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_LEX_UNRANK - Fatal error!");
            Console.WriteLine("  Input rank is illegal.");
            return null;
        }

        t = new int[k];

        rank_copy = rank;

        x = 1;

        for (i = 1; i <= k; i++)
        {
            while (typeMethods.i4_choose(n - x, k - i) <= rank_copy)
            {
                rank_copy -= typeMethods.i4_choose(n - x, k - i);
                x += 1;
            }

            t[i - 1] = x;
            x += 1;
        }

        return t;
    }

    public static int ksubset_revdoor_rank(int k, int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_RANK computes the revolving door rank of a K subset.
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
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input, int T[K], describes a K subset.  T(I) is the I-th
        //    element.  The elements must be listed in ascending order.
        // 
        //    Output, int KSUBSET_REVDOOR_RANK, the rank of the K subset.
        // 
    {
        bool check;
        int i;
        int rank;
        int s;
        // 
        //  Check.
        // 
        check = ksubset_lex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_REVDOOR_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return 0;
        }

        rank = (k % 2) switch
        {
            0 => 0,
            _ => -1
        };

        s = 1;

        for (i = k; 1 <= i; i--)
        {
            rank += s * typeMethods.i4_choose(t[i - 1], i);
            s = -s;
        }

        return rank;
    }

    public static void ksubset_revdoor_successor(int k, int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_SUCCESSOR computes the K subset revolving door successor.
        // 
        //  Discussion:
        // 
        //    After numerous attempts to implement the algorithm published in
        //    Kreher and Stinson, the Nijenhuis and Wilf version was implemented
        //    instead.  The K and S algorithm is supposedly based on the N and W one.
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
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms for Computers and Calculators,
        //    Second Edition,
        //    Academic Press, 1978,
        //    ISBN: 0-12-519260-6,
        //    LC: QA164.N54.
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input/output, int T[K], describes a K subset.  T(I) is the
        //    I-th element.  The elements must be listed in ascending order.
        //    On input, T describes a K subset.
        //    On output, T describes the next K subset in the ordering.
        //    If the input T was the last in the ordering, then the output T
        //    will be the first.
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
        int j;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
                typeMethods.i4vec_indicator1(k, ref t);
                rank = 0;
                return;
        }

        // 
        //  Check.
        // 
        check = ksubset_lex_check(k, n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_REVDOOR_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        j = 0;

        for (;;)
        {
            if (0 < j || k % 2 == 0)
            {
                j += 1;

                if (k < j)
                {
                    t[k - 1] = k;
                    rank = 0;
                    return;
                }

                if (t[j - 1] != j)
                {
                    t[j - 1] -= 1;

                    if (j != 1)
                    {
                        t[j - 2] = j - 1;
                    }

                    rank += 1;
                    return;
                }
            }

            j += 1;

            if (j < k)
            {
                if (t[j - 1] != t[j] - 1)
                {
                    break;
                }
            }
            else
            {
                if (t[j - 1] != n)
                {
                    break;
                }
            }
        }

        t[j - 1] += 1;

        if (j != 1)
        {
            t[j - 2] = t[j - 1] - 1;
        }

        rank += 1;
    }

    public static int[] ksubset_revdoor_unrank(int rank, int k, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_UNRANK computes the K subset of given revolving door rank.
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
        //    Input, int RANK, the rank of the K subset.
        // 
        //    Input, int K, the number of elements each K subset must
        //    have.  1 <= K <= N.
        // 
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Output, int KSUBSET_REVDOOR_UNRANK[K], describes the K subset of the given
        //    rank.  T(I) is the I-th element.  The elements must be listed in
        //    ascending order.
        // 
    {
        int i;
        int nksub;
        int rank_copy;
        int[] t;
        int x;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KSUBSET_REVDOOR_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
        }

        if (k < 1 || n < k)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_REVDOOR_UNRANK - Fatal error!");
            Console.WriteLine("  Input K is illegal.");
            return null;
        }

        nksub = ksubset_enum(k, n);

        if (rank < 0 || nksub < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("KSUBSET_REVDOOR_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        rank_copy = rank;

        t = new int[k];

        x = n;

        for (i = k; 1 <= i; i--)
        {
            while (rank_copy < typeMethods.i4_choose(x, i))
            {
                x -= 1;
            }

            t[i - 1] = x + 1;
            rank_copy = typeMethods.i4_choose(x + 1, i) - rank_copy - 1;
        }

        return t;
    }
}