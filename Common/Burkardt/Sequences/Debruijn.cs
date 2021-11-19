using System;
using Burkardt.Graph;

namespace Burkardt.Sequence;

public static class Debruijn
{
    public static void debruijn(int m, int n, ref int[] string_)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DEBRUIJN constructs a de Bruijn sequence.
        //
        //  Discussion:
        //
        //    Suppose we have an alphabet of M letters, and we are interested in
        //    all possible strings of length N.  If M = 2 and N = 3, then we are
        //    interested in the M**N strings:
        //
        //      000
        //      001
        //      010
        //      011
        //      100
        //      101
        //      110
        //      111
        //
        //    Now, instead of making a list like this, we prefer, if possible, to
        //    write a string of letters, such that every consecutive sequence of
        //    N letters is one of the strings, and every string occurs once, if
        //    we allow wraparound.
        //
        //    For the above example, a suitable sequence would be the 8 characters:
        //
        //      00011101(00...
        //
        //    where we have suggested the wraparound feature by repeating the first
        //    two characters at the end.
        //
        //    Such a sequence is called a de Bruijn sequence.  It can easily be
        //    constructed by considering a directed graph, whose nodes are all
        //    M**(N-1) strings of length N-1.  A node I has a directed edge to
        //    node J (labeled with character K) if the string at node J can
        //    be constructed by beheading the string at node I and adding character K.
        //
        //    In this setting, a de Bruijn sequence is simply an Eulerian circuit
        //    of the directed graph, with the edge labels being the entries of the
        //    sequence.  In general, there are many distinct de Bruijn sequences
        //    for the same parameter M and N.  This program will only find one
        //    of them.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 June 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of letters in the alphabet.
        //
        //    Input, int N, the number of letters in a codeword.
        //
        //    Output, int STRING[M**N], a deBruijn string.
        //
    {
        int i;
        bool success = false;
        //
        //  Construct the adjacency information.
        //
        int nnode = (int)Math.Pow(m, n - 1);
        int nedge = (int)Math.Pow(m, n);

        int[] inode = new int[nedge];
        int[] ivec = new int[n - 1];
        int[] jnode = new int[nedge];
        int[] jvec = new int[n - 1];
        int[] knode = new int[nedge];

        int iedge = 0;

        for (i = 1; i <= nnode; i++)
        {
            IndexNS.Index.index_unrank0(n - 1, m, i, ref ivec);

            int k;
            for (k = 1; k <= m; k++)
            {
                //
                //  Shift N-2 entries of IVEC down.
                //
                int j;
                for (j = 0; j < n - 2; j++)
                {
                    jvec[j] = ivec[j + 1];
                }

                jvec[n - 2] = k;

                j = IndexNS.Index.index_rank0(n - 1, m, jvec);

                inode[iedge] = i;
                jnode[iedge] = j;
                knode[iedge] = k;
                iedge += 1;
            }
        }

        //
        //  Determine a circuit.
        //
        int[] trail = new int[nedge];

        Digraph.digraph_arc_euler(nnode, nedge, inode, jnode, ref success, ref trail);
        //
        //  The string is constructed from the labels of the edges in the circuit.
        //
        for (i = 0; i < nedge; i++)
        {
            string_[i] = knode[trail[i] - 1];
        }
    }
}