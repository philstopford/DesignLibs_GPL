using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool cycle_check(int n, int ncycle, int[] t, int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    CYCLE_CHECK checks a permutation in cycle form.
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
        //    Input, int N, the number of items permuted.
        //    N must be positive.
        // 
        //    Input, int NCYCLE, the number of cycles.
        //    1 <= NCYCLE <= N.
        // 
        //    Input, int T[N], INDEX[NCYCLE], describes the permutation
        //    as a collection of NCYCLE cycles.  The first cycle is
        //    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
        //
        //    Output, bool CYCLE_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        //
    {
        int i;
        int iseek;

        switch (n)
        {
            // 
            //  N must be at least 1.
            // 
            case < 1:
                return false;
        }

        // 
        //  1 <= NCYCLE <= N.
        // 
        if (ncycle < 1 || n < ncycle)
        {
            return false;
        }

        // 
        //  1 <= INDEX(I) <= N.
        //
        for (i = 0; i < ncycle; i++)
        {
            if (index[i] < 1 || n < index[i])
            {
                return false;
            }
        }

        //
        //  The INDEX values sum to N.
        //
        if (typeMethods.i4vec_sum(ncycle, index) != n)
        {
            return false;
        }

        // 
        //  1 <= T(I) <= N.
        // 
        for (i = 0; i < n; i++)
        {
            if (t[i] < 1 || n < t[i])
            {
                return false;
            }
        }

        // 
        //  Verify that every value from 1 to N occurs in T.
        // 
        for (iseek = 1; iseek <= n; iseek++)
        {
            int ifind = -1;

            for (i = 0; i < n; i++)
            {
                if (t[i] != iseek)
                {
                    continue;
                }

                ifind = i + 1;
                break;
            }

            switch (ifind)
            {
                case -1:
                    return false;
            }
        }

        return true;
    }

    public static int[] cycle_to_perm(int n, int ncycle, int[] t, int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    CYCLE_TO_PERM converts a permutation from cycle to array form.
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
        //    Input, int N, the number of items permuted.
        //    N must be positive.
        // 
        //    Input, int NCYCLE, the number of cycles.
        //    1 <= NCYCLE <= N.
        // 
        //    Input, int T[N], INDEX[NCYCLE], describes the permutation
        //    as a collection of NCYCLE cycles.  The first cycle is
        //    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
        // 
        //    Output, int CYCLE_TO_PERM[N], describes the permutation using a
        //    single array.  For each index I, I -> P(I).
        // 
    {
        int i;
        // 
        //  Check.
        // 
        bool check = cycle_check(n, ncycle, t, index);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("CYCLE_TO_PERM - Fatal error!");
                Console.WriteLine("  The cycle is not legal.");
                return null;
        }

        int[] p = new int[n];

        int jhi = 0;

        for (i = 1; i <= ncycle; i++)
        {
            int jlo = jhi + 1;
            jhi += index[i - 1];

            int j;
            for (j = jlo; j <= jhi; j++)
            {
                if (j < jhi)
                {
                    p[t[j - 1] - 1] = t[j];
                }
                else
                {
                    p[t[j - 1] - 1] = t[jlo - 1];
                }
            }
        }

        return p;
    }
}