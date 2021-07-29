using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Treepack
{
    public class TreeNextData
    {
        public VecNextData vData { get; set; }

        public TreeNextData()
        {
            vData = new VecNextData();
        }
    }
    public static class Tree
    {
        public static void tree_arc_center(int nnode, int[] inode, int[] jnode, ref int[] center,
                ref int eccent, ref int parity)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ARC_CENTER computes the center, eccentricity, and parity of a tree.
            //
            //  Discussion:
            //
            //    A tree is an undirected graph of N nodes, which uses N-1 edges,
            //    and is connected.  
            //
            //    A graph with N-1 edges is not guaranteed to be a tree, and so this
            //    routine must first check that condition before proceeding.
            //
            //    The edge distance between two nodes I and J is the minimum number of
            //    edges that must be traversed in a path from I and J.
            //
            //    The eccentricity of a node I is the maximum edge distance between
            //    node I and the other nodes J in the graph.
            //
            //    The radius of a graph is the minimum eccentricity over all nodes
            //    in the graph.
            //
            //    The diameter of a graph is the maximum eccentricity over all nodes
            //    in the graph.
            //
            //    The center of a graph is the set of nodes whose eccentricity is 
            //    equal to the radius, that is, the set of nodes of minimum eccentricity.
            //
            //    For a tree, the center is either a single node, or a pair of
            //    neighbor nodes.
            //
            //    The parity of the tree is 1 if the center is a single node, or 2 if
            //    the center is 2 nodes.
            //
            //    The center of a tree can be found by removing all "leaves", that is,
            //    nodes of degree 1.  This step is repeated until only 1 or 2 nodes
            //    are left.
            //
            //    Thanks to Alexander Sax for pointing out that a previous version of the
            //    code was failing when the tree had an odd parity, that is, a single
            //    center node, 15 April 2013.
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
            //    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges of
            //    the tree.  Edge I connects nodes INODE(I) and JNODE(I).
            //
            //    Output, int CENTER[2].  CENTER(1) is the index of the
            //    first node in the center.  CENTER(2) is 0 if there is only one node
            //    in the center, or else the index of the second node.
            //
            //    Output, ref int ECCENT, the eccentricity of the nodes in 
            //    the center, and the radius of the the tree.
            //
            //    Output, ref int PARITY, the parity of the tree, which is
            //    normally 1 or 2.
            //
        {
            int[] degree;
            int i;
            int iedge;
            int ileaf;
            int j;
            int[] list;
            int nedge;
            int nleaf;
            int nnode2;
            int result;

            eccent = 0;
            center[0] = 0;
            center[1] = 0;
            parity = 0;

            if (nnode <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ARC_CENTER - Fatal error!");
                Console.WriteLine("  NNODE <= 0.");
                return;
            }
            else if (nnode == 1)
            {
                eccent = 0;
                center[0] = 1;
                center[1] = 0;
                parity = 1;
                return;
            }
            else if (nnode == 2)
            {
                eccent = 1;
                center[0] = 1;
                center[1] = 2;
                parity = 2;
                return;
            }

            //
            //  Is this graph really a tree?
            //
            nedge = nnode - 1;
            result = Graph.graph_arc_is_tree(nedge, inode, jnode, nnode);

            if (result == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ARC_CENTER - Fatal error!");
                Console.WriteLine("  This graph is NOT a tree.");
                return;
            }

            //
            //  Compute the degrees.
            //
            degree = Graph.graph_arc_degree(nnode, nedge, inode, jnode);
            //
            //  Defoliate the tree.
            //
            nnode2 = nnode;
            list = new int[nnode];

            for (;;)
            {
                eccent = eccent + 1;
                //
                //  Find and mark the leaves.
                //
                nleaf = 0;

                for (i = 1; i <= nnode; i++)
                {
                    if (degree[i - 1] == 1)
                    {
                        nleaf = nleaf + 1;
                        list[nleaf - 1] = i;
                    }
                }

                //
                //  Delete the leaves.
                //
                for (ileaf = 1; ileaf <= nleaf; ileaf++)
                {
                    i = list[ileaf - 1];

                    iedge = 0;
                    j = 0;

                    for (;;)
                    {
                        iedge = iedge + 1;

                        if (nedge < iedge)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("TREE_ARC_CENTER - Fatal error!");
                            Console.WriteLine("  Data or algorithm failure.");
                            return;
                        }

                        if (inode[iedge - 1] == i)
                        {
                            j = jnode[iedge - 1];
                            inode[iedge - 1] = -inode[iedge - 1];
                            jnode[iedge - 1] = -jnode[iedge - 1];
                        }
                        else if (jnode[iedge - 1] == i)
                        {
                            j = inode[iedge - 1];
                            inode[iedge - 1] = -inode[iedge - 1];
                            jnode[iedge - 1] = -jnode[iedge - 1];
                        }

                        if (j != 0)
                        {
                            break;
                        }
                    }

                    degree[i - 1] = -1;
                    nnode2 = nnode2 - 1;
                    degree[j - 1] = degree[j - 1] - 1;
                    //
                    //  If the other node has degree 0, we must have just finished
                    //  stripping all leaves from the tree, leaving a single node.
                    //  Don't kill it here.  It is our odd center.
                    //
                    //     if ( degree(j) == 0 )
                    //     {
                    //       nnode2 = nnode2 - 1
                    //     }
                    //
                }

                //
                //  Find the remaining nodes.
                //
                nnode2 = 0;

                for (i = 1; i <= nnode; i++)
                {
                    if (0 <= degree[i - 1])
                    {
                        nnode2 = nnode2 + 1;
                        list[nnode2 - 1] = i;
                    }
                }

                //
                //  If at least 3, more pruning is required.
                //
                if (nnode2 < 3)
                {
                    break;
                }
            }

            //
            //  If only one or two nodes left, we are done.
            //
            parity = nnode2;

            for (i = 0; i < nnode2; i++)
            {
                center[i] = list[i];
            }

            for (i = 0; i < nedge; i++)
            {
                inode[i] = Math.Abs(inode[i]);
                jnode[i] = Math.Abs(jnode[i]);
            }
        }

        public static void tree_arc_diam(int nnode, int[] inode, int[] jnode, ref int diam,
                ref int[] label, ref int n1, ref int n2)

            //****************************************************************************80
            /*
            Purpose:
            
            TREE_ARC_DIAM computes the "diameter" of a tree.
            
            Discussion:
            
            A tree is an undirected graph of N nodes, which uses N-1 edges,
            and is connected.  
            
            A graph with N-1 edges is not guaranteed to be a tree, and so this
            routine must first check that condition before proceeding.
            
            The diameter of a graph is the length of the longest possible
            path that never repeats an edge.
            
            Licensing:
            
            This code is distributed under the GNU LGPL license. 
            
            Modified:
            
            20 August 2019
            
            Author:
            
            John Burkardt
            
            Parameters:
            
            Input, int NNODE, the number of nodes.
            
            Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges 
            of the tree.  Edge I connects nodes INODE(I) and JNODE(I).
            
            Output, ref int DIAM, the length of the longest path 
            in the tree.
            
            Output, int LABEL[NNODE], marks the path between 
            nodes N1 and N2.  Node I is in this path if LABEL(I) is 1.
            
            Output, ref int N1, &N2, the indices of two nodes in the 
            tree which are separated by DIAM edges.
            */
        {
            int[] degree;
            int i;
            int invals;
            int j;
            int k;
            int kstep;
            int nabe;
            int nedge;
            int result;

            for (i = 0; i < nnode; i++)
            {
                label[i] = i + 1;
            }

            if (nnode <= 0)
            {
                diam = 0;
                Console.WriteLine("");
                Console.WriteLine("TREE_ARC_DIAM - Fatal error!");
                Console.WriteLine("  NNODE <= 0.");
                return;
            }

            if (nnode == 1)
            {
                diam = 0;
                n1 = 1;
                n2 = 1;
                return;
            }

            //
            //  Is this graph really a tree?
            //
            nedge = nnode - 1;
            result = Graph.graph_arc_is_tree(nedge, inode, jnode, nnode);

            if (result == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ARC_DIAM - Fatal error!");
                Console.WriteLine("  This graph is NOT a tree.");
                return;
            }

            //
            //  On step K:
            //
            //    Identify the terminal and interior nodes.
            //
            //    If there are no interior nodes left, 
            //
            //      then there are just two nodes left at all.  The diameter is 2*K-1, 
            //      and a maximal path extends between the nodes whose labels are 
            //      contained in the two remaining terminal nodes.
            //
            //    Else
            //
            //      The label of each terminal node is passed to its interior neighbor.
            //      If more than one label arrives, take any one.
            //
            //      The terminal nodes are removed.
            //
            kstep = 0;
            degree = new int[nnode];

            for (;;)
            {
                kstep = kstep + 1;
                //
                //  Compute the degree of each node.
                //
                for (j = 1; j <= nnode; j++)
                {
                    degree[j - 1] = 0;
                }

                for (j = 1; j <= nedge; j++)
                {
                    k = inode[j - 1];
                    if (0 < k)
                    {
                        degree[k - 1] = degree[k - 1] + 1;
                    }

                    k = jnode[j - 1];
                    if (0 < k)
                    {
                        degree[k - 1] = degree[k - 1] + 1;
                    }
                }

                //
                //  Count the number of interior nodes.
                //
                invals = 0;
                for (i = 1; i <= nnode; i++)
                {
                    if (1 < degree[i - 1])
                    {
                        invals = invals + 1;
                    }
                }

                //
                //  If there are 1 or 0 interior nodes, it's time to stop.
                //
                if (invals == 1)
                {
                    diam = 2 * kstep;
                    break;
                }
                else if (invals == 0)
                {
                    diam = 2 * kstep - 1;
                    break;
                }

                //
                //  If there are at least two interior nodes, then chop off the 
                //  terminal nodes and pass their labels inward.
                //
                for (k = 1; k <= nnode; k++)
                {
                    if (degree[k - 1] == 1)
                    {
                        for (j = 1; j <= nedge; j++)
                        {
                            if (inode[j - 1] == k)
                            {
                                nabe = jnode[j - 1];
                                label[nabe - 1] = label[k - 1];
                                inode[j - 1] = -inode[j - 1];
                                jnode[j - 1] = -jnode[j - 1];
                            }
                            else if (jnode[j - 1] == k)
                            {
                                nabe = inode[j - 1];
                                label[nabe - 1] = label[k - 1];
                                inode[j - 1] = -inode[j - 1];
                                jnode[j - 1] = -jnode[j - 1];
                            }
                        }
                    }
                }
            }

            //
            //  Now get the labels from two of the remaining terminal nodes.
            //  The nodes represented by these labels will be a diameter apart.
            //
            n1 = 0;
            n2 = 0;

            for (i = 1; i <= nnode; i++)
            {
                if (degree[i - 1] == 1)
                {
                    if (n1 == 0)
                    {
                        n1 = label[i - 1];
                    }
                    else if (n2 == 0)
                    {
                        n2 = label[i - 1];
                    }
                }
            }

            //
            //  Set the labels of the interior node (if any) and nodes marked
            //  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
            //
            if (invals == 1)
            {
                for (i = 1; i <= nnode; i++)
                {
                    if (1 < degree[i - 1])
                    {
                        label[i - 1] = 1;
                    }
                }
            }

            for (i = 1; i <= nnode; i++)
            {
                if (label[i - 1] == n1 || label[i - 1] == n2)
                {
                    label[i - 1] = 1;
                }
                else
                {
                    label[i - 1] = 0;
                }
            }

            //
            //  Clean up the arrays.
            //
            for (j = 1; j <= nedge; j++)
            {
                inode[j - 1] = Math.Abs(inode[j - 1]);
                k = inode[j - 1];
                jnode[j - 1] = Math.Abs(jnode[j - 1]);
                k = jnode[j - 1];
            }
        }

        public static void tree_arc_random(int nnode, ref int seed, ref int[] code, ref int[] inode,
                ref int[] jnode)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
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
            //    Input/output, ref int SEED, a seed for the random number 
            //    generator.
            //
            //    Output, int CODE[NNODE-2], the Pruefer code for the 
            //    labeled tree.
            //
            //    Output, int INODE[NNODE-1], JNODE[NNODE-1], the edge 
            //    array for the tree.
            //
        {
            int i;

            if (nnode <= 2)
            {
                return;
            }

            Vec.vec_random(nnode - 2, nnode, ref seed, code);

            for (i = 0; i < nnode - 2; i++)
            {
                code[i] = code[i] + 1;
            }

            Pruefer.pruefer_to_tree_arc(nnode, code, ref inode, ref jnode);

        }

        public static int[] tree_arc_to_pruefer(int nnode, int[] inode, int[] jnode)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ARC_TO_PRUEFER is given a labeled tree, and computes its Pruefer code.
            //
            //  Discussion:
            //
            //    A tree is an undirected graph of N nodes, which uses N-1 edges,
            //    and is connected.  
            //
            //    A graph with N-1 edges is not guaranteed to be a tree, and so this
            //    routine must first check that condition before proceeding.
            //
            //    The Pruefer code is a correspondence between all labeled trees of
            //    N nodes, and all list of N-2 integers between 1 and N (with repetition
            //    allowed).  The number of labeled trees on N nodes is therefore N^(N-2).
            //
            //    The Pruefer code is constructed from the tree as follows:
            //
            //    A terminal node on the tree is defined as a node with only one neighbor.
            //
            //    Consider the set of all terminal nodes on the tree.  Take the one
            //    with the highest label, I.  Record the label of its neighbor, J.
            //    Delete node I and the edge between node I and J.
            //
            //    J is the first entry in the Pruefer code for the tree.   Repeat
            //    the operation a total of N-2 times to get the complete code.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2013
            //
            //  Reference:
            //
            //    Dennis Stanton, Dennis White,
            //    Constructive Combinatorics,
            //    Springer Verlage, New York, 1986.
            //
            //  Parameters:
            //
            //    Input, int NNODE, the number of nodes.
            //
            //    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edge array 
            //    of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
            //
            //    Output, int TREE_ARC_TO_PRUEFER[NNODE-2], the Pruefer code of the tree.
            //
        {
            int[] code;
            int[] degree;
            int i;
            int i2;
            int iterm;
            int j;
            int jsave = 0;
            int nedge;
            int result;
            //
            //  Is this graph really a tree?
            //
            nedge = nnode - 1;
            result = Graph.graph_arc_is_tree(nedge, inode, jnode, nnode);

            if (result == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ARC_TO_PRUEFER - Fatal error!");
                Console.WriteLine("  This graph is NOT a tree.");
                return null;
            }

            code = new int[nnode - 2];
            //
            //  Compute the degree of each node.
            //
            nedge = nnode - 1;
            degree = Graph.graph_arc_degree(nnode, nedge, inode, jnode);
            //
            //  Compute the next term of the Pruefer code.
            //
            for (i = 1; i <= nnode - 2; i++)
            {
                //
                //  Find the terminal node with the highest label.
                //
                iterm = 0;

                for (j = 1; j <= nnode; j++)
                {
                    if (degree[j - 1] == 1)
                    {
                        iterm = j;
                    }
                }

                //
                //  Find the edge that includes this node, and note the
                //  index of the other node.
                //
                i2 = -1;

                for (j = 1; j <= nnode - 1; j++)
                {
                    jsave = j;

                    if (inode[j - 1] == iterm)
                    {
                        i2 = 2;
                        break;
                    }
                    else if (jnode[j - 1] == iterm)
                    {
                        i2 = 1;
                        break;
                    }
                }

                if (i2 == -1)
                {

                    Console.WriteLine("");
                    Console.WriteLine("TREE_ARC_TO_PRUEFER - Fatal error!");
                    Console.WriteLine("  Algorithm breakdown!");
                    return null;
                }

                //
                //  Delete the edge from the tree.
                //
                degree[inode[jsave - 1] - 1] = degree[inode[jsave - 1] - 1] - 1;
                degree[jnode[jsave - 1] - 1] = degree[jnode[jsave - 1] - 1] - 1;
                //
                //  Add the neighbor of the node to the Pruefer code.
                //
                if (i2 == 1)
                {
                    code[i - 1] = inode[jsave - 1];
                }
                else
                {
                    code[i - 1] = jnode[jsave - 1];
                }

                //
                //  Negate the nodes in the edge list to mark them as used.
                //
                inode[jsave - 1] = -inode[jsave - 1];
                jnode[jsave - 1] = -jnode[jsave - 1];
            }

            //
            //  Before returning, restore the original form of the edge list.
            //
            for (i = 1; i <= nnode - 1; i++)
            {
                inode[i - 1] = Math.Abs(inode[i - 1]);
                jnode[i - 1] = Math.Abs(jnode[i - 1]);
            }

            return code;
        }

        public static int tree_enum(int nnode)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ENUM enumerates the labeled trees on NNODE nodes.
            //
            //  Discussion:
            //
            //    The formula is due to Cauchy.
            //
            //  Example:
            //
            //    NNODE      NTREE
            //
            //    0              1
            //    1              1
            //    2              1
            //    3              3
            //    4             16
            //    5            125
            //    6           1296
            //    7          16807
            //    8         262144
            //    9        4782969
            //   10      100000000
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
            //    Input, int NNODE, the number of nodes in each tree.
            //    NNODE must normally be at least 3, but for this routine,
            //    any value of NNODE is allowed.  Values of NNODE greater than 10
            //    will probably overflow.
            //
            //    Output, int TREE_ENUM, the number of distinct labeled trees.
            //
        {
            int ntree;

            if (nnode < 0)
            {
                ntree = 0;
            }
            else if (nnode == 0)
            {
                ntree = 1;
            }
            else if (nnode == 1)
            {
                ntree = 1;
            }
            else if (nnode == 2)
            {
                ntree = 1;
            }
            else
            {
                ntree = (int) Math.Pow(nnode, nnode - 2);
            }

            return ntree;
        }

        public static void tree_parent_next(ref TreeNextData data, int nnode, ref int[] code, ref int[] itree, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_PARENT_NEXT generates, one at a time, all labeled trees.
            //
            //  Discussion:
            //
            //    The routine also returns the corresponding Pruefer codes.
            //
            //  Formula:
            //
            //    There are N^(N-2) labeled trees on N nodes (Cayley's formula).
            //
            //    The number of trees in which node I has degree D(I) is the
            //    multinomial coefficient: ( N-2; D(1)-1, D(2)-1, ..., D(N)-1 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2013
            //
            //  Parameters:
            //
            //    Input, int NNODE, the number of nodes to be used in 
            //    the trees.
            //
            //    Input/output, int CODE[NNODE].  The first NNODE-2 entries 
            //    of CODE contain the Pruefer code for the given labeled tree.
            //
            //    Output, int ITREE[NNODE].  The first NNODE-1 entries 
            //    of ITREE describe the edges that go between the nodes.  Each pair
            //    (I, ITREE(I)) represents an edge.  Thus if ITREE(5) = 3,
            //    there is an edge from node 3 to node 5.
            //
            //    Input/output, ref int MORE.  On the first call only, the
            //    user is required to set MORE = .FALSE.  Then call the routine, and
            //    it will return information about the first tree
            //    as well as setting MORE to the value .TRUE.
            //    Keep calling to get another tree until MORE is .FALSE.
            //    on return, at which point there are no more trees.
            //
        {
            int i;

            if (more)
            {
                for (i = 0; i < nnode - 2; i++)
                {
                    code[i] = code[i] - 1;
                }
            }

            VecNextData vData = data.vData;
            
            Vec.vec_next(ref vData, nnode - 2, nnode, ref code, ref more);

            data.vData = vData;

            for (i = 0; i < nnode - 2; i++)
            {
                code[i] = code[i] + 1;
            }

            Pruefer.pruefer_to_tree_2(nnode, code, ref itree);

        }

        public static void tree_parent_to_arc(int nnode, int[] parent, ref int nedge, ref int[] inode,
                ref int[] jnode)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_PARENT_TO_ARC converts a tree from parent to arc representation.
            //
            //  Discussion:
            //
            //    Parent representation lists the parent node of each node.  For a
            //    tree of N nodes, one node has a parent of 0, representing a null link.
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
            //    Input, int NNODE, the number of nodes in the tree.
            //
            //    Input, int PARENT[NNODE], the parent node representation 
            //    of the tree.
            //
            //    Output, ref int NEDGE, the number of edges, normally NNODE-1.
            //
            //    Output, int INODE[NEDGE], JNODE[NEDGE], pairs of nodes
            //    that define the links.
            //
        {
            int i;

            nedge = 0;

            for (i = 1; i <= nnode; i++)
            {
                if (parent[i - 1] != 0)
                {
                    nedge = nedge + 1;
                    inode[nedge - 1] = i;
                    jnode[nedge - 1] = parent[i - 1];
                }
            }

        }

        public static int tree_rb_enum(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_RB_ENUM returns the number of rooted binary trees with N nodes.
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
            //    Input, int N, the number of nodes in the rooted 
            //    binary tree.  N should be odd.
            //
            //    Output, int TREE_RB_ENUM, the number of rooted binary trees 
            //    with N nodes.
            //
        {
            int[] c;
            int m;
            int num;

            if (n < 0)
            {
                num = 0;
            }
            else if (n == 0)
            {
                num = 1;
            }
            else if ((n % 2) == 0)
            {
                num = 0;
            }
            else
            {
                m = (n - 1) / 2;
                c = Catalan.catalan(m);
                num = c[m];
            }

            return num;
        }

        public static void tree_rb_lex_next(int n, int[] a, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_RB_LEX_NEXT generates rooted binary trees in lexicographic order.
            //
            //  Discussion:
            //
            //    The information definining the tree of N nodes is stored in a vector 
            //    of 0's and 1's, in preorder traversal form.  Essentially, the
            //    shape of the tree is traced out with a pencil that starts at the root,
            //    and stops at the very last null leaf.  The first time that a (non-null) 
            //    node is encountered, a 1 is added to the vector, and the left 
            //    descendant of the node is visited next.  When the path returns from
            //    the first descendant, the second descendant is visited.  When then path
            //    returns again, the path passes back up from the node to its parent.
            //    A null leaf is encountered only once, and causes a zero to be added to 
            //    the vector, and the path goes back up to the parent node.  
            //
            //    The lexicographic order is used to order the vectors of 1's and 0's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2013
            //
            //  Reference:
            //
            //    Frank Ruskey,
            //    Combinatorial Generation,
            //    To appear.
            //
            //  Parameters:
            //
            //    Input, int N, the number of nodes in the rooted binary
            //    tree.  N should be odd.
            //
            //    Input/output, int A[N], the preorder traversal form for
            //    the previous/next rooted binary tree.
            //
            //    Output, logical &MORE, is TRUE if the next rooted binary tree was
            //    returned on this call, or FALSE if there are no more rooted binary
            //    trees, and the output of the previous call was the last one.
            //
        {
            int i;
            int k;
            int p;
            int q;

            if (!more)
            {
                for (i = 1; i <= n - 2; i = i + 2)
                {
                    a[i - 1] = 1;
                }

                for (i = 2; i <= n - 1; i = i + 2)
                {
                    a[i - 1] = 0;
                }

                a[n - 1] = 0;
                more = true;
                return;
            }

            //
            //  Find the last 1 in A.
            //
            k = n;
            while (a[k - 1] == 0)
            {
                k = k - 1;
            }

            q = n - k - 1;
            //
            //  Find the last 0 preceding the last 1 in A.
            //  If there is none, then we are done, because 11...1100..00 
            //  is the final element.
            //
            for (;;)
            {
                if (k == 1)
                {
                    more = false;
                    return;
                }

                if (a[k - 1] == 0)
                {
                    break;
                }

                k = k - 1;
            }

            p = n - k - q - 1;

            a[k - 1] = 1;
            for (i = k + 1; i <= n - 2 * p + 1; i++)
            {
                a[i - 1] = 0;
            }

            for (i = n - 2 * p + 2; i <= n - 2; i = i + 2)
            {
                a[i - 1] = 1;
            }

            for (i = n - 2 * p + 3; i <= n - 1; i = i + 2)
            {
                a[i - 1] = 0;
            }

            a[n - 1] = 0;

            return;
        }

        public static int[] tree_rb_to_parent(int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_RB_TO_PARENT converts rooted binary tree to parent node representation.
            //
            //  Discussion:
            //
            //    Parent node representation of a tree assigns to each node a "parent" node,
            //    which represents the first link of the path between the node and the 
            //    root node.  The root node itself is assigned a parent of 0.
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
            //    Input, int N, the number of nodes in the tree.
            //
            //    Input, int A[N], the preorder traversal form for the
            //    rooted binary tree.
            //
            //    Output, int TREE_RB_TO_PARENT[N], the parent node representation 
            //    of the tree.
            //
        {
            int dad;
            int k;
            int node;
            int node_num;
            int[] parent;
            int[] use;

            parent = new int[n];
            use = new int[n];

            node = 0;
            node_num = 0;

            for (k = 1; k <= n; k++)
            {
                dad = node;
                node_num = node_num + 1;
                node = node_num;
                parent[node - 1] = dad;

                if (a[k - 1] == 1)
                {
                    use[node - 1] = 0;
                }
                else
                {
                    use[node - 1] = 2;

                    while (use[node - 1] == 2)
                    {
                        node = dad;
                        if (node == 0)
                        {
                            break;
                        }

                        use[node - 1] = use[node - 1] + 1;
                        dad = parent[node - 1];
                    }
                }
            }

            return parent;
        }

        public static void tree_rb_yule(ref int n, ref int seed, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_RB_YULE adds two nodes to a rooted binary tree using the Yule model.
            //
            //  Discussion:
            //
            //    The Yule model is a simulation of how an evolutionary family tree
            //    develops.  We start with a root node.  The internal nodes of the tree 
            //    are inactive and never change.  Each pendant or leaf node of the
            //    tree represents a biological family that can spontaneously "fission",
            //    developing two new distinct sub families.  In graphical terms, the node
            //    becomes internal, with two new leaf nodes depending from it.
            //
            //    The tree is stored in inorder traversal form.
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
            //    Input/output, ref int N, the number of nodes in the input
            //    tree.  On output, this number has been increased, usually by 2.
            //
            //    Input/output, ref int SEED, a seed for the random number
            //    generator.
            //
            //    Input/output, int A[*], the preorder traversal form 
            //    for the rooted binary tree.  The number of entries in A is N.
            //
        {
            int i;
            int ileaf;
            int j;
            int jleaf;
            int nleaf;

            if (n <= 0)
            {
                n = 1;
                a[0] = 0;
                return;
            }

            //
            //  Count the expected number of leaves, which are the 0 values.
            //
            nleaf = (n + 1) / 2;
            //
            //  Choose a random number between 1 and NLEAF.
            //
            ileaf = UniformRNG.i4_uniform_ab(1, nleaf, ref seed);
            //
            //  Locate leaf number ILEAF.
            //
            j = 0;
            jleaf = 0;
            for (i = 1; i <= n; i++)
            {
                if (a[i - 1] == 0)
                {
                    jleaf = jleaf + 1;
                }

                if (jleaf == ileaf)
                {
                    j = i;
                    break;
                }
            }

            //
            //  Replace '0' by '100'
            //
            for (i = n; j <= i; i--)
            {
                a[i + 1] = a[i - 1];
            }

            a[j - 1] = 1;
            a[j] = 0;

            n = n + 2;

            return;
        }

        public static int[] tree_rooted_code(int nnode, int[] parent)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ROOTED_CODE returns the code of a rooted tree.
            //
            //  Discussion:
            //
            //    This code for a rooted tree depends on the node ordering, so it's actually
            //    the code for a labeled rooted tree.  To eliminate the effects of node
            //    labeling, one could choose as the code for a tree the maximum of all
            //    the codes associated with the different possible labelings of the tree.
            //    There are more effective ways of arriving at this code than simply
            //    generating all possible codes and comparing them.  
            //
            //    For a tree with NNODES, the code is a list of 2*NNODE 0's and 1's,
            //    describing a traversal of the tree starting at an imaginary node 0,
            //    moving "down" to the root (a code entry of 1), and then moving
            //    "down" (1) or "up" (0) as the tree is traversed in a depth first
            //    manner.  The final move must be from the root up to the imaginary
            //    node 0.
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
            //    Input, int PARENT[NNODE], is the parent node of each node.
            //    The node with parent 0 is the root.
            //
            //    Output, int TREE_ROOTED_CODE[2*NNODE], the code for the tree.
            //
        {
            int[] code;
            int father;
            int i;
            int k = 0;
            int son;

            code = new int[2 * nnode];
            //
            //  Find the root.
            //
            father = 0;
            for (i = 1; i <= nnode; i++)
            {
                if (parent[i - 1] == 0)
                {
                    k = 1;
                    code[0] = 1;
                    father = i;
                    break;
                }
            }

            if (father == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ROOTED_CODE - Fatal error!");
                Console.WriteLine("  Could not find the root.");
                return null;
            }

            while (father != 0)
            {
                k = k + 1;
                code[k - 1] = 0;
                for (son = 1; son <= nnode; son++)
                {
                    if (parent[son - 1] == father)
                    {
                        code[k - 1] = 1;
                        father = son;
                        break;
                    }
                }

                if (code[k - 1] == 0)
                {
                    parent[father - 1] = -parent[father - 1];
                    father = -parent[father - 1];
                }
            }

            for (i = 0; i < nnode; i++)
            {
                parent[i] = -parent[i];
            }

            return code;
        }

        public static int tree_rooted_code_compare(int nnode, int npart, int[] code1, int[] code2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ROOTED_CODE_COMPARE compares a portion of the code for two rooted trees.
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
            //    Input, int NPART, the number of nodes for which the code
            //    has been determined.  This determines the portion of the codes to be
            //    compared.  We expect 0 <= NPART <= NNODE.
            //
            //    Input, int CODE1[2*NNODE], CODE2[2*NNODE], the two 
            //    rooted tree codes to be compared.
            //
            //    Output, int TREE_ROOTED_CODE_COMPARE, the result of the comparison.
            //    -1, CODE1 < CODE2,
            //     0, CODE1 = CODE2,
            //    +1, CODE1 > CODE2.
            //
        {
            int i;
            int ihi;
            int result;

            result = 0;

            if (npart <= 0)
            {
                return result;
            }

            ihi = 2 * nnode;
            if (npart < nnode)
            {
                ihi = 2 * npart;
            }

            for (i = 0; i < ihi; i++)
            {
                if (code1[i] < code2[i])
                {
                    result = -1;
                    return result;
                }
                else if (code2[i] < code1[i])
                {
                    result = +1;
                    return result;
                }
            }

            return result;
        }

        public static void tree_rooted_depth(int nnode, int[] parent, ref int depth, ref int[] depth_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ROOTED_DEPTH returns the depth of a rooted tree.
            //
            //  Discussion:
            //
            //    The depth of any node of a rooted tree is the number of edges in 
            //    the shortest path from the root to the node.
            //
            //    The depth of the rooted tree is the maximum of the depths
            //    of all the nodes.
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
            //    Input, int PARENT[NNODE], is the parent node of each node.
            //    The node with parent 0 is the root.
            //
            //    Output, ref int DEPTH, the depth of the tree.
            //
            //    Output, int DEPTH_NODE[NNODE], the depth of each node.
            //
        {
            int i;
            int j;
            int root;
            //
            //  Find the root.
            //
            root = -1;
            for (i = 1; i <= nnode; i++)
            {
                if (parent[i - 1] == 0)
                {
                    root = i;
                    break;
                }
            }

            if (root == -1)
            {
                Console.WriteLine("");
                Console.WriteLine("TREE_ROOTED_DEPTH - Fatal error!");
                Console.WriteLine("  Could not find the root.");
                return;
            }

            //
            //  Determine the depth of each node by moving towards the node.
            //  If you reach a node whose depth is already known, stop early.
            //
            for (i = 0; i < nnode; i++)
            {
                depth_node[i] = 0;
            }

            for (i = 1; i <= nnode; i++)
            {
                j = i;

                while (j != root)
                {
                    depth_node[i - 1] = depth_node[i - 1] + 1;
                    j = parent[j - 1];

                    if (0 < depth_node[j - 1])
                    {
                        depth_node[i - 1] = depth_node[i - 1] + depth_node[j - 1];
                        break;
                    }
                }
            }

            //
            //  Determine the maximum depth.
            //
            depth = typeMethods.i4vec_max(nnode, depth_node);

        }

        public static int[] tree_rooted_enum(int nnode)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ROOTED_ENUM counts the number of unlabeled rooted trees.
            //
            //  Example:
            //
            //    Input    Output
            //
            //      1         1
            //      2         1
            //      3         2
            //      4         4
            //      5         9
            //      6        20
            //      7        48
            //      8       115
            //      9       286
            //     10       719
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
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int NNODE, the number of nodes.
            //
            //    Output, int TREE_ROOTED_ENUM[NNODE].  NTREE(I) is the number of 
            //    rooted, unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
            //
        {
            int i;
            int id;
            int isum;
            int itd;
            int j;
            int nlast;
            int[] ntree;

            ntree = new int[nnode];

            ntree[0] = 1;

            for (nlast = 2; nlast <= nnode; nlast++)
            {
                isum = 0;

                for (id = 1; id <= nlast - 1; id++)
                {
                    i = nlast;
                    itd = ntree[id - 1] * id;

                    for (j = 1; j <= nlast - 1; j++)
                    {
                        i = i - id;

                        if (i <= 0)
                        {
                            break;
                        }

                        isum = isum + ntree[i - 1] * itd;
                    }
                }

                ntree[nlast - 1] = isum / (nlast - 1);
            }

            return ntree;
        }

        public static int[] tree_rooted_random(int nnode, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TREE_ROOTED_RANDOM selects a random unlabeled rooted tree.
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
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int NNODE, the number of nodes.
            //
            //    Input/output, ref int SEED, a seed for the random 
            //    number generator.
            //
            //    Output, int TREE_ROOTED_RANDOM[NNODE].  (I,ITREE(I)) is the I-th edge
            //    of the output tree for I = 2,NNODE.  ITREE(1)=0.
            //
        {
            int i;
            int id;
            int is1;
            int is2;
            int itd;
            int[] itree;
            int iz;
            int j;
            int l;
            int ll;
            int ls;
            int m;
            int[] ntree;
            int nval;
            double r;
            int[] stack;

            if (nnode <= 0)
            {
                Console.WriteLine("");
                ;
                Console.WriteLine("TREE_ROOTED_RANDOM - Fatal error!");
                Console.WriteLine("  NNODE = " + nnode + "");
                Console.WriteLine("  but NNODE must be at least 1.");
                return null;
            }

            itree = new int[nnode];
            stack = new int[2 * nnode];
            //
            //  Compute a table of the number of such trees for a given number of nodes.
            //
            ntree = tree_rooted_enum(nnode);
            //
            //  Now select one such tree at random.
            //
            l = 0;

            nval = nnode;
            is1 = 0;
            is2 = 0;

            for (;;)
            {
                while (2 < nval)
                {
                    r = UniformRNG.r8_uniform_01(ref seed);

                    iz = (int) ((nval - 1) * ntree[nval - 1] * r);

                    id = 0;

                    id = id + 1;
                    itd = id * ntree[id - 1];
                    m = nval;
                    j = 0;

                    for (;;)
                    {
                        j = j + 1;
                        m = m - id;

                        if (m < 1)
                        {
                            id = id + 1;
                            itd = id * ntree[id - 1];
                            m = nval;
                            j = 0;
                            continue;
                        }

                        iz = iz - ntree[m - 1] * itd;
                        if (iz < 0)
                        {
                            break;
                        }
                    }

                    is1 = is1 + 1;
                    stack[0 + (is1 - 1) * 2] = j;
                    stack[1 + (is1 - 1) * 2] = id;
                    nval = m;
                }

                itree[is2] = l;
                l = is2 + 1;
                is2 = is2 + nval;

                if (1 < nval)
                {
                    itree[is2 - 1] = is2 - 1;
                }

                for (;;)
                {
                    nval = stack[1 + (is1 - 1) * 2];

                    if (nval != 0)
                    {
                        stack[1 + (is1 - 1) * 2] = 0;
                        break;
                    }

                    j = stack[0 + (is1 - 1) * 2];
                    is1 = is1 - 1;
                    m = is2 - l + 1;
                    ll = itree[l - 1];
                    ls = l + (j - 1) * m - 1;

                    if (j != 1)
                    {
                        for (i = l; i <= ls; i++)
                        {
                            itree[i + m - 1] = itree[i - 1] + m;
                            if (((i - l) % m) == 0)
                            {
                                itree[i + m - 1] = ll;
                            }
                        }
                    }

                    is2 = ls + m;

                    if (is2 == nnode)
                    {
                        return itree;
                    }

                    l = ll;
                }
            }

            return itree;
        }
    }
}