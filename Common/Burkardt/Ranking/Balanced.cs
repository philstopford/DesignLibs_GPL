using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool bal_seq_check(int n, int[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_CHECK checks a balanced sequence.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 November 2015
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
        //    Input, int T[2*N], a balanced sequence.
        // 
        //    Output, bool BAL_SEQ_CHECK.
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

        int one_count = 0;
        int zero_count = 0;

        for (i = 0; i < 2 * n; i++)
        {
            switch (t[i])
            {
                case 0:
                    zero_count += 1;
                    break;
                case 1:
                    one_count += 1;
                    break;
                default:
                    return false;
            }

            if (zero_count < one_count)
            {
                return false;
            }
        }

        if (one_count != zero_count)
        {
            return false;
        }

        return true;
    }

    public static int bal_seq_enum(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BAL_SEQ_ENUM enumerates the balanced sequences.
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
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be nonnegative.
        // 
        //    Output, int BAL_SEQ_ENUM, the number of balanced sequences.
        // 
    {
        int value = typeMethods.i4_choose(2 * n, n) / (n + 1);

        return value;
    }

    public static int bal_seq_rank(int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BAL_SEQ_RANK ranks a balanced sequence.
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
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be positive.
        // 
        //    Input, int T[2*N], a balanced sequence.
        // 
        //    Output, int BAL_SEQ_RANK, the rank of the balanced sequence.
        // 
    {
        int x;
        // 
        //  Check.
        // 
        bool check = bal_seq_check(n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("BAL_SEQ_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return 1;
        }

        int y = 0;
        int rank = 0;

        for (x = 1; x <= 2 * n - 1; x++)
        {
            switch (t[x - 1])
            {
                case 0:
                    y += 1;
                    break;
                default:
                    int mxy = mountain(n, x, y + 1);
                    rank += mxy;
                    y -= 1;
                    break;
            }
        }

        return rank;
    }

    public static void bal_seq_successor(int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BAL_SEQ_SUCCESSOR computes the lexical balanced sequence successor.
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
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be positive.
        // 
        //    Input/output, int T[2*N], on input, a balanced sequence,
        //    and on output, its lexical successor.
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
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                for (i = n; i < 2 * n; i++)
                {
                    t[i] = 1;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        bool check = bal_seq_check(n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("BAL_SEQ_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        // 
        //  After the I-th 0 there is a 'slot' with the capacity to
        //  hold between 0 and I ones.
        // 
        //  The first element of the sequence has all the 1''s cowering
        //  behind the N-th 0.
        // 
        //  We seek to move a 1 to the left, and to do it lexically,
        //  we will move a 1 to the rightmost slot that is under capacity.
        // 
        //  Find the slot.
        // 
        int slot = 0;
        int slot_index = 0;
        int slot_ones = 0;

        int open = 0;
        int open_index = 0;

        for (i = 1; i <= 2 * n; i++)
        {
            switch (t[i - 1])
            {
                case 0:
                {
                    switch (slot)
                    {
                        case > 0:
                        {
                            if (slot_ones < slot)
                            {
                                open = slot;
                                open_index = slot_index;
                            }

                            break;
                        }
                    }

                    slot += 1;
                    slot_index = i;
                    break;
                }
                default:
                    slot_ones += 1;
                    break;
            }
        }

        // 
        //  If OPEN is not 0, then preserve the string up to the OPEN-th 0,
        //  preserve the 1''s that follow, but then write a 1, then
        //  all the remaining 0's and all the remaining 1's.
        // 
        if (open != 0)
        {
            int j = open_index + 1;

            while (t[j - 1] == 1)
            {
                j += 1;
            }

            t[j - 1] = 1;

            for (i = open + 1; i <= n; i++)
            {
                j += 1;
                t[j - 1] = 0;
            }

            for (i = j + 1; i <= 2 * n; i++)
            {
                t[i - 1] = 1;
            }
        }
        // 
        //  If OPEN is 0, the last element was input.
        //  Return the first one.
        // 
        else
        {
            for (i = 0; i < n; i++)
            {
                t[i] = 0;
            }

            for (i = n; i < 2 * n; i++)
            {
                t[i] = 1;
            }

            rank = 0;
            return;
        }

        rank += 1;
    }

    public static int[] bal_seq_to_tableau(int n, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BAL_SEQ_TO_TABLEAU converts a balanced sequence to a 2 by N tableau.
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
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be positive.
        // 
        //    Input, int T[2*N], a balanced sequence.
        // 
        //    Output, int BAL_SEQ_TO_TABLEAU[2*N], a 2 by N tableau.
        // 
    {
        int[] c = new int[2];
        int i;
        // 
        //  Check.
        // 
        bool check = bal_seq_check(n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("BAL_SEQ_TO_TABLEAU - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return null;
        }

        int[] tab = new int[2 * n];

        c[0] = 0;
        c[1] = 0;

        for (i = 1; i <= 2 * n; i++)
        {
            int r = t[i - 1] + 1;
            c[r - 1] += 1;
            tab[r - 1 + (c[r - 1] - 1) * 2] = i;
        }

        return tab;
    }

    public static int[] bal_seq_unrank(int rank, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BAL_SEQ_UNRANK unranks a balanced sequence.
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
        //    Input, int RANK, the rank of the balanced sequence.
        // 
        //    Input, int N, the number of 0's (and 1's) in the sequence.
        //    N must be positive.
        // 
        //    Output, int BAL_SEQ_UNRANK[2*N], a balanced sequence.
        // 
    {
        int x;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("BAL_SEQ_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
        }

        int nseq = bal_seq_enum(n);

        if (rank < 0 || nseq < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("BAL_SEQ_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        int[] t = new int[2 * n];

        int y = 0;
        int low = 0;

        for (x = 0; x < 2 * n; x++)
        {
            int m = mountain(n, x + 1, y + 1);

            if (rank <= low + m - 1)
            {
                y += 1;
                t[x] = 0;
            }
            else
            {
                low += m;
                y -= 1;
                t[x] = 1;
            }
        }

        return t;
    }
}