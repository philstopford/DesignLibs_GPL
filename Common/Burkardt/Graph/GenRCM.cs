﻿using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.Graph;

public static class GenRCM
{
    public static void degree(int root, int adj_num, int[] adj_row, int[] adj, int[] mask,
            ref int[] deg, ref int iccsze, ref int[] ls, int node_num, int lsIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEGREE computes the degrees of the nodes in the connected component.
        //
        //  Discussion:
        //
        //    The connected component is specified by MASK and ROOT.
        //    Nodes for which MASK is zero are ignored.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Alan George, Joseph Liu.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Input, int ROOT, the node that defines the connected component.
        //
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int ADJ[ADJ_NUM], the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
        //    Input, int MASK[NODE_NUM], is nonzero for those nodes which are
        //    to be considered.
        //
        //    Output, int DEG[NODE_NUM], contains, for each  node in the connected
        //    component, its degree.
        //
        //    Output, int *ICCSIZE, the number of nodes in the connected component.
        //
        //    Output, int LS[NODE_NUM], stores in entries 1 through ICCSIZE the nodes
        //    in the connected component, starting with ROOT, and proceeding 
        //    by levels.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
    {
        int i;
        int node;
        //
        //  The sign of ADJ_ROW(I) is used to indicate if node I has been considered.
        //
        ls[lsIndex + 0] = root;
        adj_row[root - 1] = -adj_row[root - 1];
        int lvlend = 0;
        iccsze = 1;
        //
        //  LBEGIN is the pointer to the beginning of the current level, and
        //  LVLEND points to the end of this level.
        //
        for (;;)
        {
            int lbegin = lvlend + 1;
            lvlend = iccsze;
            //
            //  Find the degrees of nodes in the current level,
            //  and at the same time, generate the next level.
            //
            for (i = lbegin; i <= lvlend; i++)
            {
                node = ls[(lsIndex + (i - 1) + ls.Length) % ls.Length];
                int jstrt = -adj_row[(node - 1 + adj_row.Length) % adj_row.Length];
                int jstop = Math.Abs(adj_row[node % adj_row.Length]) - 1;
                int ideg = 0;

                int j;
                for (j = jstrt; j <= jstop; j++)
                {
                    int nbr = adj[(j - 1 + adj.Length) % adj.Length];

                    if (mask[(nbr - 1 + mask.Length) % mask.Length] != 0)
                    {
                        ideg += 1;

                        switch (adj_row[(nbr - 1 + adj_row.Length) % adj_row.Length])
                        {
                            case >= 0:
                                adj_row[(nbr - 1 + adj_row.Length) % adj_row.Length] = -adj_row[(nbr - 1 + adj_row.Length) % adj_row.Length];
                                iccsze += 1;
                                ls[(lsIndex + (iccsze - 1) + ls.Length) % ls.Length] = nbr;
                                break;
                        }
                    }
                }

                deg[(node - 1 + deg.Length) % deg.Length] = ideg;
            }

            //
            //  Compute the current level width.
            //
            int lvsize = iccsze - lvlend;
            //
            //  If the current level width is nonzero, generate another level.
            //
            if (lvsize == 0)
            {
                break;
            }
        }

        //
        //  Reset ADJ_ROW to its correct sign and return.
        //
        for (i = 0; i < iccsze; i++)
        {
            node = ls[(lsIndex + i + ls.Length) % ls.Length] - 1;
            adj_row[(node + adj_row.Length) % adj_row.Length] = -adj_row[(node + adj_row.Length) % adj_row.Length];
        }
    }

    public static int[] genrcm(int node_num, int adj_num, int[] adj_row, int[] adj)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GENRCM finds the reverse Cuthill-Mckee ordering for a general graph.
        //
        //  Discussion:
        //
        //    For each connected component in the graph, the routine obtains
        //    an ordering by calling RCM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 May 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Alan George, Joseph Liu.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int  ADJ[ADJ_NUM], the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
        //    Output, int GENRCM[NODE_NUM], the RCM ordering.
        //
        //  Local Parameters:
        //
        //    Local, int  LEVEL_ROW[NODE_NUM+1], the index vector for a level
        //    structure.  The level structure is stored in the currently unused
        //    spaces in the permutation vector PERM.
        //
        //    Local, int MASK[NODE_NUM], marks variables that have been numbered.
        //
    {
        int i;
        int iccsze = 0;
        int level_num = 0;
        //
        //  Assuming the input dat is 0 based, add 1 to ADJ_ROW and ADJ,
        //  because GENRCM uses 1-based indexing!
        //
        for (i = 0; i < node_num + 1; i++)
        {
            adj_row[i] += 1;
        }

        for (i = 0; i < adj_num; i++)
        {
            adj[i] += 1;
        }

        int[] perm = new int[node_num];
        int[] level_row = new int[node_num + 1];
        int[] mask = new int[node_num];

        for (i = 0; i < node_num; i++)
        {
            mask[i] = 1;
        }

        int num = 1;

        for (i = 0; i < node_num; i++)
        {
            //
            //  For each masked connected component...
            //
            if (mask[i] != 0)
            {
                int root = i + 1;
                //
                //  Find a pseudo-peripheral node ROOT.  The level structure found by
                //  ROOT_FIND is stored starting at PERM(NUM).
                //
                root_find(ref root, adj_num, adj_row, adj, mask, ref level_num,
                    ref level_row, ref perm, node_num, levelIndex: + num - 1);
                //
                //  RCM orders the component using ROOT as the starting node.
                //
                rcm(root, adj_num, adj_row, adj, ref mask, ref perm, ref iccsze,
                    node_num, permIndex: + num - 1);

                num += iccsze;
            }

            //
            //  We can stop once every node is in one of the connected components.
            //
            if (node_num < num)
            {
                break;
            }
        }

        //
        //  PERM is computed as a 1-based vector.
        //  Rewrite it as a 0-based vector.
        //
        for (i = 0; i < node_num; i++)
        {
            perm[i] -= 1;
        }

        //
        //  Subtract 1 from ADJ_ROW and ADJ because GENRCM used 1-based indexing!
        //
        for (i = 0; i < node_num + 1; i++)
        {
            adj_row[i] -= 1;
        }

        for (i = 0; i < adj_num; i++)
        {
            adj[i] -= 1;
        }

        return perm;
    }
        
    public static void level_set ( int root, int adj_num, int[] adj_row, int[] adj, ref int[] mask, 
            ref int level_num, ref int[] level_row, ref int[] level, int node_num, int levelIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_SET generates the connected level structure rooted at a given node.
        //
        //  Discussion:
        //
        //    Only nodes for which MASK is nonzero will be considered.
        //
        //    The root node chosen by the user is assigned level 1, and masked.
        //    All (unmasked) nodes reachable from a node in level 1 are
        //    assigned level 2 and masked.  The process continues until there
        //    are no unmasked nodes adjacent to any node in the current level.
        //    The number of levels may vary between 2 and NODE_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Alan George, Joseph Liu.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Input, int ROOT, the node at which the level structure
        //    is to be rooted.
        //
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int ADJ[ADJ_NUM], the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
        //    Input/output, int MASK[NODE_NUM].  On input, only nodes with nonzero
        //    MASK are to be processed.  On output, those nodes which were included
        //    in the level set have MASK set to 1.
        //
        //    Output, int *LEVEL_NUM, the number of levels in the level
        //    structure.  ROOT is in level 1.  The neighbors of ROOT
        //    are in level 2, and so on.
        //
        //    Output, int LEVEL_ROW[NODE_NUM+1], LEVEL[NODE_NUM], the rooted 
        //    level structure.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
    {
        int i;

        mask[root - 1] = 0;
        level[levelIndex + 0] = root;
        level_num = 0;
        int lvlend = 0;
        int iccsze = 1;
        //
        //  LBEGIN is the pointer to the beginning of the current level, and
        //  LVLEND points to the end of this level.
        //
        for (;;)
        {
            int lbegin = lvlend + 1;
            lvlend = iccsze;
            level_num += 1;
            level_row[level_num - 1] = lbegin;
            //
            //  Generate the next level by finding all the masked neighbors of nodes
            //  in the current level.
            //
            for (i = lbegin; i <= lvlend; i++)
            {
                int node = level[levelIndex + (i - 1)];
                int jstrt = adj_row[node - 1];
                int jstop = adj_row[node % adj_row.Length] - 1;

                int j;
                for (j = jstrt; j <= jstop; j++)
                {
                    int nbr = adj[(j - 1 + adj.Length) % adj.Length];

                    if (mask[(nbr - 1 + mask.Length) % mask.Length] == 0)
                    {
                        continue;
                    }

                    iccsze += 1;
                    level[levelIndex + (iccsze - 1)] = nbr;
                    mask[(nbr - 1 + mask.Length) % mask.Length] = 0;
                }
            }

            //
            //  Compute the current level width (the number of nodes encountered.)
            //  If it is positive, generate the next level.
            //
            int lvsize = iccsze - lvlend;

            if (lvsize <= 0)
            {
                break;
            }
        }

        level_row[level_num] = lvlend + 1;
        //
        //  Reset MASK to 1 for the nodes in the level structure.
        //
        for (i = 0; i < iccsze; i++)
        {
            mask[(level[(levelIndex + i) % level.Length] - 1 + mask.Length) % mask.Length] = 1;
        }
    }

    public static void level_set_print(int node_num, int level_num, int[] level_row,
            int[] level)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEVEL_SET_PRINT prints level set information.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int LEVEL_NUM, the number of levels.
        //
        //    Input, int LEVEL_ROW[LEVEL_NUM+1], organizes the entries of LEVEL.
        //    The entries for level I are in entries LEVEL_ROW(I)
        //    through LEVEL_ROW(I+1)-1.
        //
        //    Input, integer LEVEL[NODE_NUM], is simply a list of the nodes in an
        //    order induced by the levels.
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("LEVEL_SET_PRINT");
        Console.WriteLine("  Show the level set structure of a rooted graph.");
        Console.WriteLine("  The number of nodes is  " + node_num + "");
        Console.WriteLine("  The number of levels is " + level_num + "");
        Console.WriteLine("");
        Console.WriteLine("  Level Min Max      Nonzeros");
        Console.WriteLine("");

        for (i = 0; i < level_num; i++)
        {
            int jmin = level_row[i];
            int jmax = level_row[i + 1] - 1;

            if (jmax < jmin)
            {
                Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + jmin.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + jmax.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
            }
            else
            {
                int jlo;
                for (jlo = jmin; jlo <= jmax; jlo += 5)
                {
                    int jhi = Math.Min(jlo + 4, jmax);

                    string cout = "";
                    int j;
                    if (jlo == jmin)
                    {
                        cout = "  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                    + "  " + jmin.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                    + "  " + jmax.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                    + "   ";
                        for (j = jlo; j <= jhi; j++)
                        {
                            cout += level[j - 1].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                        }

                        Console.WriteLine(cout);
                    }
                    else
                    {
                        cout = "                     ";
                        for (j = jlo; j <= jhi; j++)
                        {
                            cout += level[j - 1].ToString(CultureInfo.InvariantCulture).PadLeft(8);
                        }

                        Console.WriteLine(cout);
                    }
                }
            }
        }
    }

    public static void rcm(int root, int adj_num, int[] adj_row, int[] adj, ref int[] mask,
            ref int[] perm, ref int iccsze, int node_num, int permIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RCM renumbers a connected component by the reverse Cuthill McKee algorithm.
        //
        //  Discussion:
        //
        //    The connected component is specified by a node ROOT and a mask.
        //    The numbering starts at the root node.
        //
        //    An outline of the algorithm is as follows:
        //
        //    X(1) = ROOT.
        //
        //    for ( I = 1 to N-1)
        //      Find all unlabeled neighbors of X(I),
        //      assign them the next available labels, in order of increasing degree.
        //
        //    When done, reverse the ordering.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Alan George, Joseph Liu.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //  Parameters:
        //
        //    Input, int ROOT, the node that defines the connected component.
        //    It is used as the starting point for the RCM ordering.
        //
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW(NODE_NUM+1).  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int ADJ(ADJ_NUM), the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
        //    Input/output, int MASK(NODE_NUM), a mask for the nodes.  Only 
        //    those nodes with nonzero input mask values are considered by the 
        //    routine.  The nodes numbered by RCM will have their mask values 
        //    set to zero.
        //
        //    Output, int PERM(NODE_NUM), the RCM ordering.
        //
        //    Output, int ICCSZE, the size of the connected component
        //    that has been numbered.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //  Local Parameters:
        //
        //    Workspace, int DEG[NODE_NUM], a temporary vector used to hold 
        //    the degree of the nodes in the section graph specified by mask and root.
        //
    {
        //
        //  Find the degrees of the nodes in the component specified by MASK and ROOT.
        //
        int[] deg = new int[node_num];

        degree(root, adj_num, adj_row, adj, mask, ref deg, ref iccsze, ref perm, node_num, permIndex);

        mask[root - 1] = 0;

        switch (iccsze)
        {
            case <= 1:
                return;
        }

        int lvlend = 0;
        int lnbr = 1;
        //
        //  LBEGIN and LVLEND point to the beginning and
        //  the end of the current level respectively.
        //
        while (lvlend < lnbr)
        {
            int lbegin = lvlend + 1;
            lvlend = lnbr;

            int i;
            for (i = lbegin; i <= lvlend; i++)
            {
                //
                //  For each node in the current level...
                //
                int node = perm[permIndex + (i - 1)];
                int jstrt = adj_row[(node - 1 + adj_row.Length) % adj_row.Length];
                int jstop = adj_row[(node + adj_row.Length) % adj_row.Length] - 1;
                //
                //  Find the unnumbered neighbors of NODE.
                //
                //  FNBR and LNBR point to the first and last neighbors
                //  of the current node in PERM.
                //
                int fnbr = lnbr + 1;

                int nbr;
                int j;
                for (j = jstrt; j <= jstop; j++)
                {
                    nbr = adj[(j - 1 + adj.Length) % adj.Length];

                    if (mask[(nbr - 1 + mask.Length) % mask.Length] == 0)
                    {
                        continue;
                    }

                    lnbr += 1;
                    mask[(nbr - 1 + mask.Length) % mask.Length] = 0;
                    perm[permIndex + (lnbr - 1)] = nbr;
                }

                //
                //  If no neighbors, skip to next node in this level.
                //
                if (lnbr <= fnbr)
                {
                    continue;
                }

                //
                //  Sort the neighbors of NODE in increasing order by degree.
                //  Linear insertion is used.
                //
                int k = fnbr;

                while (k < lnbr)
                {
                    int l = k;
                    k += 1;
                    nbr = perm[permIndex + (k - 1)];

                    while (fnbr < l)
                    {
                        int lperm = perm[permIndex + (l - 1)];

                        if (deg[(lperm - 1 + deg.Length) % deg.Length] <= deg[(nbr - 1 + deg.Length) % deg.Length])
                        {
                            break;
                        }

                        perm[l] = lperm;
                        l -= 1;
                    }

                    perm[permIndex + l] = nbr;
                }
            }
        }

        //
        //  We now have the Cuthill-McKee ordering.  Reverse it.
        //
        typeMethods.i4vec_reverse(iccsze, ref perm, permIndex);
    }
        
    public static void root_find ( ref int root, int adj_num, int[] adj_row, int[] adj, int[] mask, 
            ref int level_num, ref int[] level_row, ref int[] level, int node_num, int levelIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROOT_FIND finds a pseudo-peripheral node.
        //
        //  Discussion:
        //
        //    The diameter of a graph is the maximum distance (number of edges)
        //    between any two nodes of the graph.
        //
        //    The eccentricity of a node is the maximum distance between that
        //    node and any other node of the graph.
        //
        //    A peripheral node is a node whose eccentricity equals the
        //    diameter of the graph.
        //
        //    A pseudo-peripheral node is an approximation to a peripheral node;
        //    it may be a peripheral node, but all we know is that we tried our
        //    best.
        //
        //    The routine is given a graph, and seeks pseudo-peripheral nodes,
        //    using a modified version of the scheme of Gibbs, Poole and
        //    Stockmeyer.  It determines such a node for the section subgraph
        //    specified by MASK and ROOT.
        //
        //    The routine also determines the level structure associated with
        //    the given pseudo-peripheral node; that is, how far each node
        //    is from the pseudo-peripheral node.  The level structure is
        //    returned as a list of nodes LS, and pointers to the beginning
        //    of the list of nodes that are at a distance of 0, 1, 2, ...,
        //    NODE_NUM-1 from the pseudo-peripheral node.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Alan George, Joseph Liu.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Alan George, Joseph Liu,
        //    Computer Solution of Large Sparse Positive Definite Systems,
        //    Prentice Hall, 1981.
        //
        //    Norman Gibbs, William Poole, Paul Stockmeyer,
        //    An Algorithm for Reducing the Bandwidth and Profile of a Sparse Matrix,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 13, pages 236-250, 1976.
        //
        //    Norman Gibbs,
        //    Algorithm 509: A Hybrid Profile Reduction Algorithm,
        //    ACM Transactions on Mathematical Software,
        //    Volume 2, pages 378-387, 1976.
        //
        //  Parameters:
        //
        //    Input/output, int *ROOT.  On input, ROOT is a node in the
        //    the component of the graph for which a pseudo-peripheral node is
        //    sought.  On output, ROOT is the pseudo-peripheral node obtained.
        //
        //    Input, int ADJ_NUM, the number of adjacency entries.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about row I is stored
        //    in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ.
        //
        //    Input, int ADJ[ADJ_NUM], the adjacency structure.
        //    For each row, it contains the column indices of the nonzero entries.
        //
        //    Input, int MASK[NODE_NUM], specifies a section subgraph.  Nodes 
        //    for which MASK is zero are ignored by FNROOT.
        //
        //    Output, int *LEVEL_NUM, is the number of levels in the level structure
        //    rooted at the node ROOT.
        //
        //    Output, int LEVEL_ROW(NODE_NUM+1), LEVEL(NODE_NUM), the 
        //    level structure array pair containing the level structure found.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
    {
        int level_num2 = 0;
        //
        //  Determine the level structure rooted at ROOT.
        //
        level_set(root, adj_num, adj_row, adj, ref mask, ref level_num,
            ref level_row, ref level, node_num, levelIndex);
        //
        //  Count the number of nodes in this level structure.
        //
        int iccsze = level_row[level_num] - 1;
        switch (level_num)
        {
            //
            //  Extreme case:
            //    A complete graph has a level set of only a single level.
            //    Every node is equally good (or bad).
            //
            case 1:
                return;
        }

        //
        //  Extreme case:
        //    A "line graph" 0--0--0--0--0 has every node in its only level.
        //    By chance, we've stumbled on the ideal root.
        //
        if (level_num == iccsze)
        {
            return;
        }

        //
        //  Pick any node from the last level that has minimum degree
        //  as the starting point to generate a new level set.
        //
        for (;;)
        {
            int mindeg = iccsze;

            int jstrt = level_row[level_num - 1];
            root = level[levelIndex + (jstrt - 1)];

            if (jstrt < iccsze)
            {
                int j;
                for (j = jstrt; j <= iccsze; j++)
                {
                    int node = level[j - 1];
                    int ndeg = 0;
                    int kstrt = adj_row[node - 1];
                    int kstop = adj_row[node] - 1;

                    int k;
                    for (k = kstrt; k <= kstop; k++)
                    {
                        int nabor = adj[(k - 1 + adj.Length) % adj.Length];
                        switch (mask[(nabor - 1 + mask.Length) % mask.Length])
                        {
                            case > 0:
                                ndeg += 1;
                                break;
                        }
                    }

                    if (ndeg >= mindeg)
                    {
                        continue;
                    }

                    root = node;
                    mindeg = ndeg;
                }
            }

            //
            //  Generate the rooted level structure associated with this node.
            //
            level_set(root, adj_num, adj_row, adj, ref mask, ref level_num2,
                ref level_row, ref level, node_num, levelIndex);
            //
            //  If the number of levels did not increase, accept the new ROOT.
            //
            if (level_num2 <= level_num)
            {
                break;
            }

            level_num = level_num2;
            //
            //  In the unlikely case that ROOT is one endpoint of a line graph,
            //  we can exit now.
            //
            if (iccsze <= level_num)
            {
                break;
            }
        }
    }
}