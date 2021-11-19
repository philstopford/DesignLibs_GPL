using System;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool pruefer_check(int n, int[] p)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_CHECK checks a Pruefer code.
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
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input, int P[N-2], the Pruefer code for the tree.
        // 
        //    Output, bool PRUEFER_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        int i;

        switch (n)
        {
            case < 3:
                return false;
        }

        for (i = 0; i < n - 2; i++)
        {
            if (p[i] < 1 || n < p[i])
            {
                return false;
            }
        }

        return true;
    }

    public static int pruefer_enum(int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_ENUM enumerates the Pruefer codes on N-2 digits.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    30 November 2015
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of digits in the code, plus 2.
        // 
        //    Output, int NCODE, the number of distinct elements.
        // 
    {
        int value = n switch
        {
            < 2 => 0,
            2 => 1,
            _ => (int) Math.Pow(n, n - 2)
        };

        return value;
    }

    public static int pruefer_rank(int n, int[] p)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_RANK ranks a Pruefer code.
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
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input, int P[N-2], the Pruefer code for the tree.
        // 
        //    Output, int PRUEFER_RANK, the rank of the Pruefer code.
        // 
    {
        int i;
        // 
        //  Check.
        // 
        bool check = pruefer_check(n, p);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PRUEFER_RANK - Fatal error!");
                Console.WriteLine("  Input array is illegal.");
                return 1;
        }

        int rank = 0;
        int k = 1;
        for (i = n - 3; 0 <= i; i--)
        {
            rank += k * (p[i] - 1);
            k *= n;
        }

        return rank;
    }

    public static void pruefer_successor(int n, ref int[] p, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_SUCCESSOR computes the lexical Pruefer sequence successor.
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
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input/output, int P(N-2), on input, the Pruefer code
        //    for a tree, and on output, its lexical successor.
        // 
        //    Input/output, int RANK, the rank.
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
                for (i = 0; i < n - 2; i++)
                {
                    p[i] = 1;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        bool check = pruefer_check(n, p);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PRUEFER_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
        }

        int j = n - 2;

        for (;;)
        {
            if (p[j - 1] != n)
            {
                break;
            }

            j -= 1;

            if (j <= 0)
            {
                break;
            }
        }

        if (j != 0)
        {
            p[j - 1] += 1;
            for (i = j + 1; i <= n - 2; i++)
            {
                p[i - 1] = 1;
            }

            rank += 1;
        }
        else
        {
            for (i = 0; i < n - 2; i++)
            {
                p[i] = 1;
            }

            rank = 0;
        }
    }

    public static void pruefer_to_tree(int n, int[] p, ref int[] t )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_TO_TREE converts a Pruefer code to a tree.
        // 
        //  Discussion:
        // 
        //    The original code attempts to tack on an extra entry to P.
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
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input, int P[N-2], the Pruefer code for the tree.
        // 
        //    Output, int T[2*(N-1)], describes the edges of the tree
        //    as pairs of nodes.
        // 
    {
        int i;
        int j;
        // 
        //  Check.
        // 
        bool check = pruefer_check(n, p);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("PRUEFER_TO_TREE - Fatal error!");
                Console.WriteLine("  The input array is illegal!");
                return;
        }

        // 
        //  Initialize the tree to 0.
        // 
        for (j = 0; j < n - 1; j++)
        {
            for (i = 0; i < 2; i++)
            {
                t[i + j * 2] = 0;
            }
        }

        int[] d = new int[n];

        for (i = 0; i < n; i++)
        {
            d[i] = 1;
        }

        for (i = 0; i < n - 2; i++)
        {
            d[p[i] - 1] += 1;
        }

        for (i = 1; i <= n - 1; i++)
        {
            int x = n;
            while (d[x - 1] != 1)
            {
                x -= 1;
            }

            int y = i == n - 1 ? 1 : p[i - 1];

            d[x - 1] -= 1;
            d[y - 1] -= 1;

            t[0 + (i - 1) * 2] = x;
            t[1 + (i - 1) * 2] = y;
        }
    }

    public static int[] pruefer_to_tree_new(int n, int[] p)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_TO_TREE_NEW converts a Pruefer code to a tree.
        // 
        //  Discussion:
        // 
        //    The original code attempts to tack on an extra entry to P.
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
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Input, int P[N-2], the Pruefer code for the tree.
        // 
        //    Output, int PRUEFER_TO_TREE_NEW[2*(N-1)], describes the edges of the tree
        //    as pairs of nodes.
        // 
    {
        int[] t = new int[2 * (n - 1)];

        pruefer_to_tree(n, p, ref t);

        return t;
    }

    public static int[] pruefer_unrank(int rank, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PRUEFER_UNRANK unranks a Pruefer code.
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
        //    Input, int RANK, the rank of the Pruefer code.
        // 
        //    Input, int N, the number of nodes in the tree.
        //    N must be at least 3.
        // 
        //    Output, int P[N-2], the Pruefer code for the tree.
        // 
    {
        int i;
        int[] p;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("PRUEFER_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            case < 3:
                return null;
        }

        int ncode = pruefer_enum(n);

        if (rank < 0 || ncode < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("PRUEFER_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        int rank_copy = rank;
        p = new int[n - 2];

        for (i = n - 3; 0 <= i; i--)
        {
            p[i] = rank_copy % n + 1;
            rank_copy = (rank_copy - p[i] + 1) / n;
        }

        return p;
    }
}