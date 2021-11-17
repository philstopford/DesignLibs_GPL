using System;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool edge_check(int n_node, int n_edge, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    EDGE_CHECK checks a graph stored by edges.
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
        //  Parameters:
        // 
        //    Input, int N_NODE, the number of nodes in the graph.
        //    N_NODE must be positive.
        // 
        //    Input, int N_EDGE, the number of edges in the graph.
        //    N_EDGE must be positive.
        // 
        //    Input, int T(2,N_EDGE), describes the edges of the tree
        //    as pairs of nodes.
        // 
        //    Output, bool EDGE_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
        // 
    {
        bool check;
        int i;
        int j;
        int j2;

        check = true;

        switch (n_node)
        {
            case < 1:
                check = false;
                return check;
        }

        switch (n_edge)
        {
            case < 1:
                check = false;
                return check;
        }

        // 
        //  Every edge must join two legal nodes.
        // 
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < n_edge; j++)
            {
                if (t[i + j * 2] < 1 || n_node < t[i + j * 2])
                {
                    check = false;
                    return check;
                }
            }
        }

        // 
        //  Every edge must join distinct nodes.
        // 
        for (j = 0; j < n_edge; j++)
        {
            if (t[0 + j * 2] == t[1 + j * 2])
            {
                check = false;
                return check;
            }
        }

        // 
        //  Every edge must be distinct.
        // 
        for (j = 0; j < n_edge - 1; j++)
        {
            for (j2 = j + 1; j2 < n_edge; j2++)
            {
                if (t[0 + j * 2] == t[0 + j2 * 2] && t[1 + j * 2] == t[1 + j2 * 2])
                {
                    check = false;
                    return check;
                }

                if (t[0 + j * 2] == t[1 + j2 * 2] && t[1 + j * 2] == t[0 + j2 * 2])
                {
                    check = false;
                    return check;
                }
            }
        }

        return check;
    }

    public static int[] edge_degree(int n_node, int n_edge, int[] t)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    EDGE_DEGREE returns the degree of the nodes of a graph stored by edges.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    29 July 2011
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
        //    Input, int N_NODE, the number of nodes in the graph.
        //    N_NODE must be positive.
        // 
        //    Input, int N_EDGE, the number of edges in the graph.
        //    N_EDGE must be positive.
        // 
        //    Input, int T[2*N_EDGE], describes the edges of the tree
        //    as pairs of nodes.
        // 
        //    Output, int EDGE_DEGREE[N_NODE], the degree of each node.
        // 
    {
        bool check;
        int[] d;
        int i;
        int j;
        // 
        //  Check.
        // 
        check = edge_check(n_node, n_edge, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("EDGE_DEGREE - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return null;
        }

        // 
        //  Compute the degree of each node.
        //
        d = new int[n_node];

        for (i = 0; i < n_node; i++)
        {
            d[i] = 0;
        }

        for (j = 0; j < n_edge; j++)
        {
            d[t[0 + j * 2] - 1] += 1;
            d[t[1 + j * 2] - 1] += 1;
        }

        return d;
    }

    public static int edge_enum(int n_node)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    EDGE_ENUM enumerates the maximum number of edges in a graph on N_NODE nodes.
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
        //    Input, int N_NODE, the number of nodes in the graph.
        //    N_NODE must be positive.
        // 
        //    Output, int EDGE_ENUM, the maximum number of edges in a graph
        //    on N_NODE nodes.
        // 
    {
        int value = n_node * (n_node - 1) / 2;

        return value;
    }
}