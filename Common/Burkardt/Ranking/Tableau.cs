using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool tableau_check(int n, int[] tab)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    TABLEAU_CHECK checks a 2 by N tableau.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 December 2015
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
        //    Input, int N, the number of columns in the tableau.
        //    N must be positive.
        // 
        //    Input, int TAB[2*N], a 2 by N tableau.
        //
        //    Output, bool TABLEAU_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
    {
        bool check;
        int i;
        int j;

        check = true;

        switch (n)
        {
            case < 1:
                check = false;
                return check;
        }

        // 
        //  The entries must be between 1 and 2*N.
        // 
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (tab[i + j * 2] < 1 || 2 * n < tab[i + j * 2])
                {
                    check = false;
                    return check;
                }
            }
        }

        // 
        //  The entries must be increasing to the right.
        //
        for (i = 0; i < 2; i++)
        {
            for (j = 1; j < n; j++)
            {
                if (tab[i + j * 2] <= tab[i + (j - 1) * 2])
                {
                    check = false;
                    return check;
                }
            }
        }

        // 
        //  The entries must be increasing down.
        // 
        i = 1;
        for (j = 0; j < n; j++)
        {
            if (tab[i + j * 2] <= tab[i - 1 + j * 2])
            {
                check = false;
                return check;
            }
        }

        return check;
    }

    public static int tableau_enum(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    TABLEAU_ENUM enumerates the 2 by N standard tableaus.
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
        //    Input, int N, the number of columns in the tableau.
        //    N must be nonnegative.
        // 
        //    Output, int TABLEAU_ENUM, the number of 2 by N standard tableaus.
        // 
    {
        int value = typeMethods.i4_choose(2 * n, n) / (n + 1);

        return value;
    }

    public static int[] tableau_to_bal_seq(int n, int[] tab)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    TABLEAU_TO_BAL_SEQ converts a 2 by N tableau to a balanced sequence.
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
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be positive.
        // 
        //    Input, int TAB[2*N], a 2 by N tableau.
        // 
        //    Output, int TABLEAU_TO_BAL_SEQ[2*N], a balanced sequence.
        // 
    {
        bool check;
        int i;
        int j;
        int[] t;
        // 
        //  Check.
        // 
        check = tableau_check(n, tab);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("TABLEAU_TO_BAL_SEQ - Fatal error!");
                Console.WriteLine("  The tableau is illegal.");
                return null;
        }

        t = new int[2 * n];

        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < n; j++)
            {
                t[tab[i + j * 2] - 1] = i;
            }
        }

        return t;
    }
}