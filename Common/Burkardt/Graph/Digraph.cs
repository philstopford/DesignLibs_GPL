using System;
using Burkardt.Types;

namespace Burkardt.Graph;

public static class Digraph
{
    public static void digraph_arc_euler(int nnode, int nedge, int[] inode, int[] jnode,
            ref bool success, ref int[] trail)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGRAPH_ARC_EULER returns an Euler circuit in a digraph.
        //
        //  Discussion:
        //
        //    An Euler circuit of a digraph is a path which starts and ends at
        //    the same node and uses each directed edge exactly once.  A digraph is
        //    eulerian if it has an Euler circuit.  The problem is to decide whether
        //    a given digraph is eulerian and to find an Euler circuit if the
        //    answer is affirmative.
        //
        //  Method:
        //
        //    A digraph has an Euler circuit if and only if the number of incoming
        //    edges is equal to the number of outgoing edges at each node.
        //
        //    This characterization gives a straightforward procedure to decide whether
        //    a digraph is eulerian.  Furthermore, an Euler circuit in an eulerian
        //    digraph G of NEDGE edges can be determined by the following method:
        //
        //      STEP 1: Choose any node U as the starting node, and traverse any edge
        //      ( U, V ) incident to node U, and than traverse any unused edge incident
        //      to node U.  Repeat this process of traversing unused edges until the
        //      starting node U is reached.  Let P be the resulting walk consisting of
        //      all used edges.  If all edges of G are in P, than stop.
        //
        //      STEP 2: Choose any unused edge ( X,  Y) in G such that X is
        //      in P and Y is not in P.  Use node X as the starting node and
        //      find another walk Q using all unused edges as in step 1.
        //
        //      STEP 3: Walk P and walk Q share a common node X, they can be merged
        //      to form a walk R by starting at any node S of P and to traverse P
        //      until node X is reached; than, detour and traverse all edges of Q
        //      until node X is reached and continue to traverse the edges of P until
        //      the starting node S is reached.  Set P = R.
        //
        //      STEP 4: Repeat steps 2 and 3 until all edges are used.
        //
        //    The running time of the algorithm is O ( NEDGE ).
        //
        //    The digraph is assumed to be connected.
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
        //    Original FORTRAN77 version by Hang Tong Lau.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hang Tong Lau,
        //    Algorithms on Graphs,
        //    Tab Books, 1989.
        //
        //  Parameters:
        //
        //    Input, int NNODE, the number of nodes.
        //
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE(NEDGE); the I-th edge starts at node
        //    INODE(I) and ends at node JNODE(I).
        //
        //    Output, bool &SUCCESS, is TRUE if an Euler circuit was found.
        //
        //    Output, int TRAIL[NEDGE].  TRAIL[I] is the edge number of the I-th
        //    edge in the Euler circuit.
        //
    {
        bool[] candid;
        int[] endnod;
        int i;
        int istak;
        int j;
        int k;
        int l;
        int len;
        int lensol;
        int lenstk;
        int[] stack;
        //
        //  Check if the digraph is eulerian.
        //
        for (i = 0; i < nedge; i++)
        {
            trail[i] = 0;
        }

        endnod = new int[nedge];

        for (i = 0; i < nedge; i++)
        {
            endnod[i] = 0;
        }

        for (i = 1; i <= nedge; i++)
        {
            j = inode[i - 1];
            trail[j - 1] += 1;
            j = jnode[i - 1];
            endnod[j - 1] += 1;
        }

        for (i = 1; i <= nnode; i++)
        {
            if (trail[i - 1] != endnod[i - 1])
            {
                success = false;
                return;
            }
        }

        //
        //  The digraph is eulerian; find an Euler circuit.
        //
        success = true;
        lensol = 1;
        lenstk = 0;

        candid = new bool[nedge];
        stack = new int[2 * nedge];
        //
        //  Find the next edge.
        //
        for (;;)
        {
            switch (lensol)
            {
                case 1:
                    endnod[0] = inode[0];
                    stack[0] = 1;
                    stack[1] = 1;
                    lenstk = 2;
                    break;
                default:
                {
                    l = lensol - 1;

                    if (lensol != 2)
                    {
                        endnod[l - 1] = inode[trail[l - 1] - 1] + jnode[trail[l - 1] - 1] - endnod[l - 2];
                    }

                    k = endnod[l - 1];

                    for (i = 1; i <= nedge; i++)
                    {
                        candid[i - 1] = k == jnode[i - 1];
                    }

                    for (i = 1; i <= l; i++)
                    {
                        candid[trail[i - 1] - 1] = false;
                    }

                    len = lenstk;

                    for (i = 1; i <= nedge; i++)
                    {
                        switch (candid[i - 1])
                        {
                            case true:
                                len += 1;
                                stack[len - 1] = i;
                                break;
                        }
                    }

                    stack[len] = len - lenstk;
                    lenstk = len + 1;
                    break;
                }
            }

            for (;;)
            {
                istak = stack[lenstk - 1];
                lenstk -= 1;

                if (istak != 0)
                {
                    break;
                }

                lensol -= 1;

                switch (lensol)
                {
                    case 0:
                        typeMethods.i4vec_reverse(nedge, ref trail);
                        return;
                }
            }

            trail[lensol - 1] = stack[lenstk - 1];
            stack[lenstk - 1] = istak - 1;

            if (lensol == nedge)
            {
                break;
            }

            lensol += 1;
        }

        typeMethods.i4vec_reverse(nedge, ref trail);
    }

    public static void digraph_arc_print(int nedge, int[] inode, int[] jnode, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGRAPH_ARC_PRINT prints out a digraph from an edge list.
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
        //    Input, int NEDGE, the number of edges.
        //
        //    Input, int INODE[NEDGE], JNODE[NEDGE], the beginning and end
        //    nodes of the edges.
        //
        //    Input, string TITLE, a title.
        //
    {
        int i;

        switch (title.Length)
        {
            case > 0:
                Console.WriteLine("");
                Console.WriteLine(title + "");
                break;
        }

        for (i = 0; i < nedge; i++)
        {
            Console.WriteLine((i + 1).ToString().PadLeft(6) + "  "
                                                            + inode[i].ToString().PadLeft(6) + "  "
                                                            + jnode[i].ToString().PadLeft(6) + "");
        }

    }
}