using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Node
    {
        public static void node_merge(int dim_num, int node_num, double[] node_xy,
        double tolerance, ref int[] node_rep )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NODE_MERGE detects nodes that should be merged.
        //
        //  Discussion:
        //
        //    Two nodes "should" be merged if they are within TOLERANCE distance
        //    of each other.
        //
        //    With a tolerance of 0, only exactly equal nodes are counted.
        //
        //    With a positive tolerance, a pair of nodes inside a circle of
        //    radius TOLERANCE result in a count of 1 duplicate.
        //
        //    However, what do we do if nodes A, B and C are arranged in a line,!
        //    with A and B just within TOLERANCE of each other, and B and C just
        //    within tolerance of each other?  What we do here is make a choice
        //    that can be defended consistently.  A and B define an equivalence
        //    class because they are closer than TOLERANCE.  C is then added to
        //    this equivalence class, because it is within TOLERANCE of at least
        //    on thing in that equivalence class.
        //
        //    Thus, if 100 nodes are separated pairwise by slightly less
        //    than TOLERANCE, a total of 99 duplicates will be counted.
        //
        //    The program starts out by giving each node its own label.
        //    If it finds that two nodes should be merged, then the index of
        //    one node is used as the label for both.  This process continues
        //    until all nodes have been considered.  The number of unique nodes
        //    is the number of unique values in the output quantity NODE_REP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[DIM_NUM*NODE_NUM], the nodes.
        //
        //    Input, double TOLERANCE, the maximum distance between
        //    two nodes regarded as duplicate.
        //
        //    Output, int NODE_REP[NODE_NUM], the "representative" of each node.
        //    NODE_REP(NODE) is the index of a node which is within TOLERANCE of node
        //    NODE, or for which a chain of nodes can be found, all having the
        //    same representative, and all of which are pairwise closer than TOLERANCE.
        //
        {
            double dist;
            int i;
            int j;
            int node1;
            int node2;
            int rep;
            double[] rep_dist;

            rep_dist = new double[node_num];

            for (node1 = 0; node1 < node_num; node1++)
            {
                node_rep[node1] = node1;
            }

            for (node1 = 0; node1 < node_num; node1++)
            {
                for (j = 0; j < node_num; j++)
                {
                    rep_dist[j] = typeMethods.r8_huge();
                }

                for (node2 = 0; node2 < node_num; node2++)
                {
                    dist = 0.0;
                    for (i = 0; i < dim_num; i++)
                    {
                        dist = dist
                               + Math.Pow(node_xy[i + node1 * dim_num] - node_xy[i + node2 * dim_num], 2);
                    }

                    dist = Math.Sqrt(dist);

                    rep = node_rep[node2];

                    if (dist < rep_dist[rep])
                    {
                        rep_dist[rep] = dist;
                    }
                }

                for (node2 = 0; node2 < node_num; node2++)
                {
                    rep = node_rep[node2];
                    if (rep_dist[rep] <= tolerance)
                    {
                        node_rep[node2] = node1;
                    }
                }

            }
        }
    }
}