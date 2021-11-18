using System;
using Burkardt.Types;

namespace Burkardt.Treepack;

public static class Graph
{
    public static int graph_adj_edge_count(int[] adj, int nnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ADJ_EDGE_COUNT counts the number of edges in a graph.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ADJ[NNODE*NNODE], the adjacency information.
        //    ADJ(I,J) is 1 if there is an edge from node I to node J.
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Output, int GRAPH_ADJ_EDGE_COUNT, the number of edges in the graph.
        //
    {
        int i;
        int j;
        int nedge;

        nedge = 0;

        for (i = 0; i < nnode; i++)
        {
            for (j = 0; j < nnode; j++)
            {
                if (i == j)
                {
                    nedge += 2 * adj[i + j * nnode];
                }
                else
                {
                    nedge += adj[i + j * nnode];
                }
            }
        }

        nedge /= 2;

        return nedge;
    }

    public static int graph_adj_is_node_connected(int[] adj, int nnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ADJ_IS_NODE_CONNECTED determines if a graph is nodewise connected.
        //
        //  Discussion:
        //
        //    A graph is nodewise connected if, from every node, there is a path
        //    to any other node.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the 
        //    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Output, int GRAPH_ADJ_IS_NODE_CONNECTED.
        //    0, the graph is not nodewise connected.
        //    1, the graph is nodewise connected.
        //
    {
        int[] found;
        int i;
        int ihi;
        int ii;
        int ilo;
        int j;
        int jhi;
        int jlo;
        int[] list;
        int result;
        //
        //  FOUND(I) is 1 if node I has been reached.
        //  LIST(I) contains a list of the nodes as they are reached.
        //
        list = new int[nnode];
        found = new int[nnode];

        for (i = 0; i < nnode; i++)
        {
            list[i] = 0;
            found[i] = 0;
        }

        //
        //  Start at node 1.
        //
        found[1 - 1] = 1;
        list[1 - 1] = 1;
        ilo = 1;
        ihi = 1;
        //
        //  From the batch of nodes found last time, LIST(ILO:IHI),
        //  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
        //
        for (;;)
        {
            jlo = ihi + 1;
            jhi = ihi;

            for (ii = ilo; ii <= ihi; ii++)
            {
                i = list[ii - 1];

                for (j = 1; j <= nnode; j++)
                {
                    if (adj[i - 1 + (j - 1) * nnode] != 0 || adj[j - 1 + (i - 1) * nnode] != 0)
                    {
                        switch (found[j - 1])
                        {
                            case 0:
                                jhi += 1;
                                list[jhi - 1] = j;
                                found[j - 1] = 1;
                                break;
                        }
                    }
                }
            }

            //
            //  If no neighbors were found, exit.
            //
            if (jhi < jlo)
            {
                break;
            }

            //
            //  If neighbors were found, then go back and find THEIR neighbors.
            //
            ilo = jlo;
            ihi = jhi;
        }

        //
        //  No more neighbors were found.  Have we reached all nodes?
        //
        if (ihi == nnode)
        {
            result = 1;
        }
        else
        {
            result = 0;
        }

        return result;
    }

    public static int graph_adj_is_tree(int[] adj, int nnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ADJ_IS_TREE determines whether a graph is a tree.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the 
        //    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Output, int GRAPH_ADJ_IS_TREE.
        //    0, the graph is not a tree.
        //    1, the graph is a tree.
        //
    {
        int nedge;
        int result;

        switch (nnode)
        {
            case <= 1:
                result = 1;
                return result;
        }

        //
        //  Every node must be connected to every other node.
        //
        result = graph_adj_is_node_connected(adj, nnode);

        switch (result)
        {
            case 0:
                return result;
        }

        //
        //  There must be exactly NNODE-1 edges.
        //
        nedge = graph_adj_edge_count(adj, nnode);

        if (nedge == nnode - 1)
        {
            result = 1;
        }
        else
        {
            result = 0;
        }

        return result;
    }

    public static int[] graph_arc_degree(int nnode, int nedge, int[] inode, int[] jnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_DEGREE determines the degree of the nodes of a graph.
        //
        //  Discussion:
        //
        //    The degree of a node is the number of edges that include the node.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE], the pairs of nodes
        //    that form the edges.
        //
        //    Output, int GRAPH_ARC_DEGREE[NNODE], the degree of each node,
        //    that is, the number of edges that include the node.
        //
    {
        int[] degree;
        int i;
        int n;

        degree = new int[nnode];

        for (i = 0; i < nnode; i++)
        {
            degree[i] = 0;
        }

        for (i = 0; i < nedge; i++)
        {
            n = inode[i];

            switch (n)
            {
                case >= 1 when n <= nnode:
                    degree[n - 1] += 1;
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("GRAPH_ARC_DEGREE - Fatal error!");
                    Console.WriteLine("  Out-of-range node value = " + n + "");
                    return null;
            }

            n = jnode[i];
            switch (n)
            {
                case >= 1 when n <= nnode:
                    degree[n - 1] += 1;
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("GRAPH_ARC_DEGREE - Fatal error!");
                    Console.WriteLine("  Out-of-range node value = " + n + "");
                    return null;
            }
        }

        return degree;
    }

    public static int graph_arc_is_tree(int nedge, int[] inode, int[] jnode, int nnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_IS_TREE determines whether a graph is a tree.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
        //    JNODE(I) are the start and end nodes of the I-th edge of the graph G.
        //
        //    Input, int NEDGE, the number of edges in the graph G.
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Output, int GRAPH_ARC_IS_TREE.
        //    0, the graph is not a tree.
        //    1, the graph is a tree.
        //
    {
        int[] adj;
        int result;

        adj = graph_arc_to_graph_adj(nedge, inode, jnode);

        result = graph_adj_is_tree(adj, nnode);

        return result;
    }

    public static int graph_arc_node_count(int nedge, int[] inode, int[] jnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_NODE_COUNT counts the number of nodes in a graph.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
        //    JNODE(I) are the start and end nodes of the I-th edge.
        //
        //    Output, int GRAPH_ARC_NODE_COUNT, the number of distinct nodes.
        //
    {
        int i;
        int[] knode;
        int nnode;
        //
        //  Copy all the node labels into KNODE,
        //  sort KNODE,
        //  count the unique entries.  
        //
        //  That's NNODE.
        //
        knode = new int[2 * nedge];

        for (i = 0; i < nedge; i++)
        {
            knode[i] = inode[i];
            knode[i + nedge] = jnode[i];
        }

        typeMethods.i4vec_sort_heap_a(2 * nedge, ref knode);

        nnode = typeMethods.i4vec_sorted_unique_count(2 * nedge, knode);

        return nnode;
    }

    public static int graph_arc_node_max(int nedge, int[] inode, int[] jnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_NODE_MAX determines the maximum node label.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
        //    JNODE(I) are the start and end nodes of the I-th edge.
        //
        //    Output, int GRAPH_ARC_NODE_MAX, the maximum node index.
        //
    {
        int i;
        int node_max;

        node_max = 0;
        for (i = 0; i < nedge; i++)
        {
            if (node_max < inode[i])
            {
                node_max = inode[i];
            }

            if (node_max < jnode[i])
            {
                node_max = jnode[i];
            }
        }

        return node_max;
    }

    public static void graph_arc_print(int nedge, int[] inode, int[] jnode, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_PRINT prints out a graph from an edge list.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE], the beginning 
        //    and end nodes of the edges.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        Console.WriteLine("");

        for (i = 0; i < nedge; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + inode[i].ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + jnode[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }
    }

    public static int[] graph_arc_to_graph_adj(int nedge, int[] inode, int[] jnode)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRAPH_ARC_TO_GRAPH_ADJ converts an arc list graph to an adjacency graph.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE], the edge array for 
        //    an undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
        //
        //    Output, int GRAPH_ARC_TO_GRAPH_ADJ[NNODE*NNODE], the adjacency information.
        //
    {
        int[] adj;
        int i;
        int j;
        int k;
        int nnode;
        //
        //  Determine the number of nodes.
        //
        nnode = graph_arc_node_count(nedge, inode, jnode);

        adj = new int[nnode * nnode];

        for (j = 0; j < nnode; j++)
        {
            for (i = 0; i < nnode; i++)
            {
                adj[i + j * nnode] = 0;
            }
        }

        for (k = 0; k < nedge; k++)
        {
            i = inode[k] - 1;
            j = jnode[k] - 1;
            adj[i + j * nnode] = 1;
            adj[j + i * nnode] = 1;
        }

        return adj;
    }
}